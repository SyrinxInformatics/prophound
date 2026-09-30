# SPDX-License-Identifier: Apache-2.0
#
# Copyright 2025 Jacob Alford <jalford0000@gmail.com>
#                and Benjamin Alford <ssolar5.dev@gmail.com>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""G5: recall retained on fragmented input, prophound versus geNomad and Jaeger.

docs/claims-and-evidence.md "G5 preregistration" fixes the design, the endpoint
and the bar; this script implements it and nothing is tuned here. Truth is
LAMBDA's curated clean subset (39 chromosomes, 149 elements). Every tool runs on
the complete chromosome and on identical shreds (shred_genomes.py --dist
lognormal; fragment ids `<sequence>__<offset>`), and every call is mapped back to
chromosome coordinates before scoring, so all conditions share one truth table.

Subcommands:
  collect-genomad  geNomad end-to-end output dirs -> one call table. `virus_summary`
                   rows are the calls: proviruses at their coordinates, whole-contig
                   viruses at full contig length (on fragments, an element that fills
                   its contig is reported as a virus, not a provirus). Also emits the
                   find-proviruses-only calls as tool `geNomad find-proviruses`.
  collect-jaeger   Jaeger predict -p output dirs -> one call table. Non-rejected
                   prophage regions [sstart, eend] (raw_start/raw_end where no att
                   site was found), plus contigs classified `phage` at full length,
                   which replace any region on the same contig.
  score            covered50 per tool x condition, retention (shredded / complete),
                   paired genome-bootstrap intervals for prophound minus each
                   comparator, and prophound at each comparator's per-genome call count.

Call table columns: tool, condition, assembly, sequence, start, end, score
(0-based half-open chromosome coordinates; score only where the tool gives one).

Example:
  python scripts/g5_score.py collect-genomad --root ../g5_genomad --out analysis/g5/genomad_calls.tsv
  python scripts/g5_score.py collect-jaeger --root analysis/g5/raw/jaeger --out analysis/g5/jaeger_calls.tsv
  python scripts/g5_score.py score --lambda_dir analysis/external_benchmarks/lambda \\
      --prophound complete=analysis/lambda/ens --prophound mag_s0=analysis/g5/prophound/mag_s0 ... \\
      --calls analysis/g5/genomad_calls.tsv --calls analysis/g5/jaeger_calls.tsv --out analysis/g5/results.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
SEP = "__"
B = 10_000
COLS = ["tool", "condition", "assembly", "sequence", "start", "end", "score"]
PRIMARY = "prophound default"
COMPARATORS = ["geNomad", "Jaeger 57341_1.5M"]


def split_frag(name: str) -> tuple[str, int]:
    """`CP002907__12345` -> (CP002907, 12345); an unshredded name has offset 0."""
    seq, _, off = str(name).partition(SEP)
    return seq, int(off) if off else 0


# ---- collectors ---------------------------------------------------------------

def collect_genomad(root: Path) -> pd.DataFrame:
    """root/<condition>/<assembly>/<assembly>_summary/<assembly>_virus_summary.tsv"""
    rows = []
    for vs in sorted(root.glob("*/*/*_summary/*_virus_summary.tsv")):
        asm, cond = vs.parent.parent.name, vs.parent.parent.parent.name
        t = pd.read_csv(vs, sep="\t")
        for r in t.itertuples():
            name = str(r.seq_name)
            contig = name.split("|provirus_")[0]
            seq, off = split_frag(contig)
            if "|provirus_" in name:
                s, e = (int(x) for x in str(r.coordinates).split("-"))
                s, e = s - 1, e
            else:
                s, e = 0, int(r.length)
            rows.append(("geNomad", cond, asm, seq, off + s, off + e, float(r.virus_score)))
        pv = vs.parent.parent / f"{asm}_find_proviruses" / f"{asm}_provirus.tsv"
        if pv.exists():
            for r in pd.read_csv(pv, sep="\t").itertuples():
                seq, off = split_frag(r.source_seq if hasattr(r, "source_seq") else str(r.seq_name).split("|provirus_")[0])
                rows.append(("geNomad find-proviruses", cond, asm, seq, off + int(r.start) - 1, off + int(r.end), np.nan))
    if not rows:
        raise SystemExit(f"no <condition>/<assembly>/*_summary/*_virus_summary.tsv under {root}")
    return pd.DataFrame(rows, columns=COLS)


def collect_jaeger(root: Path) -> pd.DataFrame:
    """root/<model label>/<condition>/<assembly>/<model dir>/{<asm>*.tsv, <asm>_prophages/prophages_jaeger.tsv}"""
    rows = []
    for asm_dir in sorted(root.glob("*/*/*/*")):
        if not asm_dir.is_dir():
            continue
        asm, cond, label = asm_dir.parent.name, asm_dir.parent.parent.name, asm_dir.parent.parent.parent.name
        tool = f"Jaeger {label}"
        summ = [p for p in asm_dir.glob(f"{asm}*.tsv") if "phages" not in p.name]
        if len(summ) != 1:
            raise SystemExit(f"expected one contig summary in {asm_dir}, found {[p.name for p in summ]}")
        cls = pd.read_csv(summ[0], sep="\t")
        whole = cls[cls.prediction == "phage"]
        whole_ids = set(whole.contig_id.astype(str))
        for r in whole.itertuples():
            seq, off = split_frag(r.contig_id)
            rows.append((tool, cond, asm, seq, off, off + int(r.length), np.nan))
        pps = list(asm_dir.glob("*prophages/prophages_jaeger.tsv"))   # `<asm>_prophages` or `<asm>_default_prophages`
        if len(pps) > 1:
            raise SystemExit(f"several prophage tables in {asm_dir}")
        pp = pps[0] if pps else asm_dir / "none"
        if pp.exists() and pp.stat().st_size:
            t = pd.read_csv(pp, sep="\t")
            if "reject" in t:
                t = t[~t.reject.astype(str).str.lower().eq("true")]
            for r in t.itertuples():
                if str(r.contig_id) in whole_ids:
                    continue
                s = r.sstart if pd.notna(r.sstart) else r.raw_start
                e = r.eend if pd.notna(r.eend) else r.raw_end
                seq, off = split_frag(r.contig_id)
                rows.append((tool, cond, asm, seq, off + int(s), off + int(e), float(r.phage_score)))
    if not rows:
        raise SystemExit(f"no Jaeger outputs under {root}")
    return pd.DataFrame(rows, columns=COLS)


def prophound_calls(cond: str, run_dir: str) -> pd.DataFrame:
    import mode_recheck as mr
    ens = mr.load_tracks(run_dir)
    seeds = [mr.load_tracks(run_dir, f"prob_ckpt{k}") for k in (1, 2, 3)]
    out = []
    for mode in mr.MODE_THRESHOLDS:
        reg = mr.shipped_calls(ens, seeds, mode)
        so = reg.contig_id.astype(str).map(split_frag)
        seq = [s.split(".")[0] for s, _ in so]
        off = np.array([o for _, o in so], dtype=int)
        out.append(pd.DataFrame({"tool": f"prophound {mode}", "condition": cond, "assembly": reg.genome_id.values,
                                 "sequence": seq, "start": reg.start.values.astype(int) + off,
                                 "end": reg.end.values.astype(int) + off, "score": reg.max_prob.values}))
    return pd.concat(out, ignore_index=True)


# ---- scoring ------------------------------------------------------------------

def cov_frac(truth: pd.DataFrame, calls: pd.DataFrame) -> np.ndarray:
    """Per element, the fraction of its bp under the union of calls on its sequence."""
    out = np.zeros(len(truth))
    by = {k: g.sort_values("start")[["start", "end"]].to_numpy() for k, g in calls.groupby(["assembly", "sequence"])}
    for i, r in enumerate(truth.itertuples()):
        c = by.get((r.assembly, r.sequence))
        if c is None:
            continue
        c = c[(c[:, 1] > r.start) & (c[:, 0] < r.end)]
        bp, cs, ce = 0, None, None
        for s, e in c:
            s, e = max(s, r.start), min(e, r.end)
            if ce is None or s > ce:
                bp += 0 if ce is None else ce - cs
                cs, ce = s, e
            else:
                ce = max(ce, e)
        bp += 0 if ce is None else ce - cs
        out[i] = bp / (r.end - r.start)
    return out


def kept_frac(truth: pd.DataFrame, shred_dir: Path) -> np.ndarray:
    """Per element, the fraction of its bp present in the shredded input at all."""
    frags = []
    for f in shred_dir.glob("*.fna"):
        name, n = None, 0
        for line in open(f):
            if line.startswith(">"):
                if name:
                    frags.append((f.stem, *split_frag(name), n))
                name, n = line[1:].split()[0], 0
            else:
                n += len(line.strip())
        if name:
            frags.append((f.stem, *split_frag(name), n))
    fr = pd.DataFrame(frags, columns=["assembly", "sequence", "start", "len"])
    fr["end"] = fr.start + fr.len
    return cov_frac(truth, fr)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    g = sub.add_parser("collect-genomad")
    g.add_argument("--root", required=True)
    g.add_argument("--out", required=True)
    j = sub.add_parser("collect-jaeger")
    j.add_argument("--root", required=True)
    j.add_argument("--out", required=True)
    s = sub.add_parser("score")
    s.add_argument("--lambda_dir", required=True, help="lambda_prepare.py output (truth, genomes)")
    s.add_argument("--prophound", action="append", default=[], metavar="COND=RUN_DIR")
    s.add_argument("--calls", action="append", default=[], help="collected call tables")
    s.add_argument("--shreds", default=None, help="dir of <condition>/ shred FASTAs, for the reachable ceiling")
    s.add_argument("--seed", type=int, default=20260930)
    s.add_argument("--out", required=True)
    a = ap.parse_args()

    if a.cmd in ("collect-genomad", "collect-jaeger"):
        df = (collect_genomad if a.cmd == "collect-genomad" else collect_jaeger)(Path(a.root))
        Path(a.out).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(a.out, sep="\t", index=False)
        print(df.groupby(["tool", "condition"]).assembly.agg(["nunique", "size"]).to_string())
        print(f"[done] wrote {a.out}")
        return

    D = Path(a.lambda_dir)
    G = pd.read_csv(D / "lambda_genomes.tsv", sep="\t")
    T = pd.read_csv(D / "lambda_truth.tsv", sep="\t")
    T["kb"] = (T.end - T.start) / 1000
    calls = [pd.read_csv(c, sep="\t") for c in a.calls]
    calls += [prophound_calls(*x.split("=", 1)) for x in a.prophound]
    C = pd.concat(calls, ignore_index=True)
    genomes = list(G.assembly)
    mb = G.set_index("assembly").genome_length / 1e6

    # per tool x condition: every one must cover all 39 genomes (a missing genome is a failed run, not a miss)
    res = {"n_genomes": len(G), "n_elements": len(T), "cells": {}}
    hit = {}                                   # (tool, condition) -> per-element covered50
    for (tool, cond), c in C.groupby(["tool", "condition"]):
        missing = set(genomes) - set(c.assembly)
        f = cov_frac(T, c)
        hit[(tool, cond)] = f >= 0.5
        res["cells"][f"{tool} | {cond}"] = {
            "covered": float((f > 0).mean()), "covered50": float((f >= 0.5).mean()), "k50": int((f >= 0.5).sum()),
            "calls_per_genome": len(c) / len(genomes), "calls_per_mb": len(c) / mb.sum(),
            "genomes_without_calls": len(missing)}
    if a.shreds:
        res["reachable_covered50"] = {p.name: float((kept_frac(T, p) >= 0.5).mean())
                                      for p in sorted(Path(a.shreds).iterdir()) if p.is_dir()}

    # conditions: `complete`, `mag` (s0-s2 pooled: elements x 3), `draft` (s0)
    def pooled(tool: str, prefix: str):
        keys = [k for k in hit if k[0] == tool and k[1].startswith(prefix)]
        return [hit[k] for k in sorted(keys)]

    tools = sorted({k[0] for k in hit})
    el_genome = T.assembly.to_numpy()
    idx_of = {g: np.flatnonzero(el_genome == g) for g in genomes}

    def per_genome(vs: list) -> np.ndarray:
        """found counts per genome, summed over replicate shreds -> (n_genomes,)"""
        return np.array([sum(v[idx_of[g]].sum() for v in vs) for g in genomes], float)

    den_el = np.array([len(idx_of[g]) for g in genomes], float)
    rng = np.random.default_rng(a.seed)
    boot_idx = rng.integers(0, len(genomes), size=(B, len(genomes)))
    summary = {}
    for regime in ("mag", "draft"):
        for tool in tools:
            comp, shr = pooled(tool, "complete"), pooled(tool, regime)
            if not comp or not shr:
                continue
            nc, ns = per_genome(comp), per_genome(shr)
            rep = len(shr)
            c50, s50 = nc.sum() / den_el.sum(), ns.sum() / (rep * den_el.sum())
            bc = nc[boot_idx].sum(1) / den_el[boot_idx].sum(1)
            bs = ns[boot_idx].sum(1) / (rep * den_el[boot_idx].sum(1))
            summary.setdefault(regime, {})[tool] = {
                "complete_covered50": c50, "shredded_covered50": s50, "retention": s50 / c50, "n_shreds": rep,
                "_boot_ret": np.divide(bs, bc, out=np.full_like(bs, np.nan), where=bc > 0), "_boot_abs": bs,
                "shredded_calls_per_genome": float(np.mean([res["cells"][f"{tool} | {k[1]}"]["calls_per_genome"]
                                                            for k in hit if k[0] == tool and k[1].startswith(regime)]))}

    # paired differences, prophound minus comparator, same genome resamples
    for regime, d in summary.items():
        for mode in [t for t in d if t.startswith("prophound")]:
            for other in [t for t in d if not t.startswith("prophound")]:
                dr = d[mode]["_boot_ret"] - d[other]["_boot_ret"]
                da = d[mode]["_boot_abs"] - d[other]["_boot_abs"]
                d[mode].setdefault("vs", {})[other] = {
                    "retention_diff": d[mode]["retention"] - d[other]["retention"],
                    "retention_diff_boot95": [float(np.nanpercentile(dr, 2.5)), float(np.nanpercentile(dr, 97.5))],
                    "abs_diff": d[mode]["shredded_covered50"] - d[other]["shredded_covered50"],
                    "abs_diff_boot95": [float(np.nanpercentile(da, 2.5)), float(np.nanpercentile(da, 97.5))]}
        for t in d.values():
            t["retention_boot95"] = [float(np.nanpercentile(t["_boot_ret"], 2.5)), float(np.nanpercentile(t["_boot_ret"], 97.5))]
            t["shredded_covered50_boot95"] = [float(np.nanpercentile(t["_boot_abs"], 2.5)),
                                              float(np.nanpercentile(t["_boot_abs"], 97.5))]
            del t["_boot_ret"], t["_boot_abs"]
    res["summary"] = summary

    # matched cost: prophound Exploratory calls ranked by score, cut per genome to the comparator's count
    mc = {}
    for (tool, cond), c in C.groupby(["tool", "condition"]):
        if tool.startswith("prophound"):
            continue
        p = C[(C.tool == "prophound exploratory") & (C.condition == cond)].sort_values("score", ascending=False)
        if p.empty:
            continue
        k = c.groupby("assembly").size()
        cut = pd.concat([p[p.assembly == g].head(int(k.get(g, 0))) for g in genomes])
        mc[f"{tool} | {cond}"] = {"competitor_covered50": res["cells"][f"{tool} | {cond}"]["covered50"],
                                  "prophound_at_matched_cost": float((cov_frac(T, cut) >= 0.5).mean())}
    res["matched_cost"] = mc

    # secondary subsets, pooled the same way
    subsets = {"not_in_corpus": T.assembly.isin(G[G.corpus_status == "not_in_corpus"].assembly).to_numpy(),
               "under_10kb": (T.kb < 10).to_numpy()}
    sec = {}
    for name, m in subsets.items():
        for tool in tools:
            for regime in ("complete", "mag", "draft"):
                vs = pooled(tool, regime)
                if vs:
                    sec.setdefault(name, {}).setdefault(tool, {})[regime] = float(np.mean([v[m].mean() for v in vs]))
        sec[name]["n_elements"] = int(m.sum())
    res["secondary"] = sec

    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(res, indent=2) + "\n")
    rows = [{"cell": k, **v} for k, v in res["cells"].items()]
    pd.DataFrame(rows).to_csv(out.with_suffix(".csv"), index=False)
    for regime, d in summary.items():
        print(f"[{regime}]")
        for tool, v in sorted(d.items(), key=lambda kv: -kv[1]["retention"]):
            print(f"  {tool:28s} complete {v['complete_covered50']:.3f} shredded {v['shredded_covered50']:.3f} "
                  f"retention {v['retention']:.3f} [{v['retention_boot95'][0]:.3f}, {v['retention_boot95'][1]:.3f}] "
                  f"calls/genome {v['shredded_calls_per_genome']:.1f}")
        if PRIMARY in d:
            for other, v in d[PRIMARY].get("vs", {}).items():
                print(f"  {PRIMARY} - {other}: retention {v['retention_diff']:+.3f} {v['retention_diff_boot95']}, "
                      f"absolute {v['abs_diff']:+.3f} {v['abs_diff_boot95']}")
    print(f"[done] wrote {out}")


if __name__ == "__main__":
    main()
