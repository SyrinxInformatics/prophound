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
"""Score prophound and every LAMBDA-published tool on the clean subset with one rule.

docs/claims-and-evidence.md "Release gate". The bars were fixed on 2026-09-28,
before any prophound output on LAMBDA existed, and nothing is tuned on LAMBDA.

Rule: covered50, meaning an element counts as found when >=50% of its bp lie
under the union of a tool's calls. Every tool is scored only on the genomes it
was run on (lambda_tools_run.tsv), so "not run" is never scored as "missed".
Intervals are 95% genome bootstrap (10,000 resamples).

  per tool      covered, covered50 (with interval), calls/genome, calls/Mb, genomes run
  G1 complement geNomad ∪ X on genomes where both ran: elements gained over geNomad alone
  G2 matched    prophound's calls ranked by max_prob and cut, per genome, to a competitor's
     cost       call count on that genome; covered50 against that competitor
  G3 recall     prophound's covered50 against the best published tool
  landing       prophound calls overlapping no curated element, split by the LAMBDA
                candidate class they overlap (LIKELY / POSSIBLE / WEAK PHAGE, GI, ...)
Subsets: all 39, not_in_corpus (34), test (5), and elements under 10 kb. With
--containment, also elements by 21-mer containment against training positives.

prophound calls are rebuilt from the run's tracks per shipped mode
(mode_recheck.shipped_calls), which needs windows.tsv written with --write-seed-scores.

Example:
  python scripts/lambda_score.py --data analysis/external_benchmarks/lambda \\
      --prophound_run analysis/lambda/ens --out analysis/external_benchmarks/lambda/scores.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

B = 10_000


def union_cov(truth: pd.DataFrame, calls: pd.DataFrame) -> np.ndarray:
    """Per element, the fraction of its bp under the union of calls on its sequence."""
    out = np.zeros(len(truth))
    by = {k: g.sort_values("start")[["start", "end"]].to_numpy() for k, g in calls.groupby("sequence")}
    for i, r in enumerate(truth.itertuples()):
        c = by.get(r.sequence)
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


def boot(num: np.ndarray, den: np.ndarray, rng) -> list:
    if len(num) == 0 or den.sum() == 0:
        return [float("nan")] * 2
    idx = rng.integers(0, len(num), size=(B, len(num)))
    r = num[idx].sum(1) / np.maximum(1, den[idx].sum(1))
    return [float(np.percentile(r, 2.5)), float(np.percentile(r, 97.5))]


def score(truth, calls, genomes, rng) -> dict:
    gs = list(genomes.assembly)
    t = truth[truth.assembly.isin(gs)].reset_index(drop=True)
    c = calls[calls.assembly.isin(gs)]
    f = union_cov(t, c)
    h50 = f >= 0.5
    num = np.array([h50[(t.assembly == g).to_numpy()].sum() for g in gs], float)
    den = np.array([(t.assembly == g).sum() for g in gs], float)
    mb = genomes.genome_length.sum() / 1e6
    return {"n_genomes": len(gs), "n_elements": int(len(t)), "covered": float((f > 0).mean()) if len(t) else None,
            "covered50": float(h50.mean()) if len(t) else None, "k50": int(h50.sum()),
            "covered50_boot95": boot(num, den, rng), "calls_per_genome": len(c) / max(1, len(gs)),
            "calls_per_mb": len(c) / max(1e-9, mb), "found": set(zip(t.assembly[h50], t.start[h50], t.end[h50]))}


def prophound_calls(run_dir: str, genomes: pd.DataFrame, genomes_dir: str | None = None) -> dict:
    """Calls per shipped mode. With `genomes_dir`, every stage `predict` runs since 0.1.0: edge
    snap and the N filter, refinement (from the run's fine_windows.tsv) and att boundaries."""
    import evalkit as ek
    ens = ek.load_tracks(run_dir)
    seeds = [ek.load_tracks(run_dir, f"prob_ckpt{k}") for k in (1, 2, 3)]
    missing = set(genomes.assembly) - set(ens)
    if missing:
        raise SystemExit(f"prophound run lacks {len(missing)} genomes, e.g. {sorted(missing)[:3]}")
    fine = None
    if genomes_dir is not None:
        fine = {g: pd.read_csv(Path(run_dir) / g / "fine_windows.tsv", sep="\t", float_precision="round_trip")
                for g in ens if (Path(run_dir) / g / "fine_windows.tsv").exists()}
        if set(fine) != set(ens):
            raise SystemExit(f"{len(set(ens) - set(fine))} genomes lack fine_windows.tsv: run with refinement on")
    out = {}
    for mode in ek.MODE_THRESHOLDS:
        reg = ek.shipped_calls(ens, seeds, mode, genomes_dir, fine_tracks=fine)
        reg = reg.assign(assembly=reg.genome_id, sequence=reg.contig_id.astype(str).str.split(".").str[0])
        out[f"prophound {mode}"] = reg[["assembly", "sequence", "start", "end", "max_prob"]]
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", required=True, help="lambda_prepare.py output dir")
    ap.add_argument("--prophound_run", default=None, help="run dir: <assembly>/windows.tsv with --write-seed-scores")
    ap.add_argument("--genomes_dir", default=None,
                    help="<assembly>.fna: rebuild with every 0.1.0 stage (snap, N filter, refinement, att)")
    ap.add_argument("--containment", default=None, help="eval_containment.py output with a `lambda` set")
    ap.add_argument("--seed", type=int, default=20260929)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    D = Path(a.data)
    G = pd.read_csv(D / "lambda_genomes.tsv", sep="\t")
    T = pd.read_csv(D / "lambda_truth.tsv", sep="\t")
    C = pd.read_csv(D / "lambda_calls.tsv", sep="\t")
    RAN = pd.read_csv(D / "lambda_tools_run.tsv", sep="\t")
    CAND = pd.read_csv(D / "lambda_candidates.tsv", sep="\t")
    T["kb"] = (T.end - T.start) / 1000

    tools = {t: C[C.tool == t] for t in sorted(C.tool.unique())}
    ran = {t: set(RAN[RAN.tool == t].assembly) for t in tools}
    if a.prophound_run:
        for k, v in prophound_calls(a.prophound_run, G, a.genomes_dir).items():
            tools[k], ran[k] = v, set(G.assembly)

    subsets = {"all": G, "not_in_corpus": G[G.corpus_status == "not_in_corpus"], "test": G[G.corpus_status == "test"]}
    truth_sets = {"all_sizes": T, "under_10kb": T[T.kb < 10]}
    if a.containment:
        ct = pd.read_csv(a.containment, sep="\t")
        ct = ct[ct.set == "lambda"].rename(columns={"accession": "assembly", "contig": "sequence"})
        ct["start"] = ct.start - 1
        T2 = T.merge(ct[["assembly", "sequence", "start", "end", "cont_train_positive"]], how="left")
        truth_sets["containment_lt_0.1"] = T2[T2.cont_train_positive < 0.1]
        truth_sets["containment_ge_0.5"] = T2[T2.cont_train_positive >= 0.5]

    res, rows = {"n_genomes": len(G), "n_elements": len(T), "per_tool": {}}, []
    for sname, g in subsets.items():
        for tname, tt in truth_sets.items():
            for tool, calls in tools.items():
                gg = g[g.assembly.isin(ran[tool])]
                sc = score(tt, calls, gg, rng)
                found = sc.pop("found")
                res["per_tool"].setdefault(f"{sname}/{tname}", {})[tool] = sc
                rows.append({"subset": sname, "truth": tname, "tool": tool, **{k: v for k, v in sc.items() if not isinstance(v, list)},
                             "boot_lo": sc["covered50_boot95"][0], "boot_hi": sc["covered50_boot95"][1]})
                if sname == "all" and tname == "all_sizes":
                    res.setdefault("_found", {})[tool] = found
    found = res.pop("_found")

    # G1: complement to geNomad, on genomes where both ran
    comp = {}
    base = "geNomad"
    for tool in tools:
        if tool == base:
            continue
        both = G[G.assembly.isin(ran[tool] & ran[base])]
        n_el = int(T.assembly.isin(both.assembly).sum())
        gen = {x for x in found[base] if x[0] in set(both.assembly)}
        oth = {x for x in found[tool] if x[0] in set(both.assembly)}
        comp[tool] = {"n_genomes": len(both), "n_elements": n_el, "genomad": len(gen), "union": len(gen | oth),
                      "gained_over_genomad": len(oth - gen),
                      "extra_calls_per_genome": len(tools[tool][tools[tool].assembly.isin(both.assembly)]) / max(1, len(both))}
    res["G1_complement_to_genomad"] = comp

    # G2: matched cost, prophound cut per genome to each competitor's call count
    if a.prophound_run:
        mc = {}
        pc = tools["prophound exploratory"].sort_values("max_prob", ascending=False)
        for tool in [t for t in tools if not t.startswith("prophound")]:
            gg = G[G.assembly.isin(ran[tool])]
            k = tools[tool].groupby("assembly").size()
            cut = pd.concat([pc[pc.assembly == g].head(int(k.get(g, 0))) for g in gg.assembly] or [pc.head(0)])
            mine, theirs = score(T, cut, gg, rng), score(T, tools[tool], gg, rng)
            mc[tool] = {"n_genomes": len(gg), "competitor_covered50": theirs["covered50"],
                        "prophound_at_matched_cost_covered50": mine["covered50"],
                        "calls_per_genome": theirs["calls_per_genome"]}
        res["G2_matched_cost"] = mc
        # where off-truth calls land
        land = {}
        for mode in ("exploratory", "default", "conservative"):
            p = tools[f"prophound {mode}"]
            off = []
            for r in p.itertuples():
                tr = T[(T.sequence == r.sequence) & (T.end > r.start) & (T.start < r.end)]
                if len(tr):
                    continue
                cd = CAND[(CAND.sequence == r.sequence) & (CAND.end > r.start) & (CAND.start < r.end)]
                off.append(cd.classification.iloc[0] if len(cd) else "no LAMBDA candidate")
            land[mode] = {"off_truth_calls": len(off), "by_class": pd.Series(off).value_counts().to_dict() if off else {}}
        res["off_truth_landing"] = land

    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(res, indent=2, default=str) + "\n")
    pd.DataFrame(rows).to_csv(out.with_suffix(".csv"), index=False)
    head = res["per_tool"]["all/all_sizes"]
    for tool, s in sorted(head.items(), key=lambda kv: -(kv[1]["covered50"] or 0)):
        print(f"{tool:24s} genomes {s['n_genomes']:2d} elements {s['n_elements']:3d} covered50 {s['covered50']:.3f} "
              f"[{s['covered50_boot95'][0]:.3f}, {s['covered50_boot95'][1]:.3f}] calls/genome {s['calls_per_genome']:.1f}",
              flush=True)
    print("[G1] geNomad + X, elements gained:", {k: v["gained_over_genomad"] for k, v in
                                                  sorted(comp.items(), key=lambda kv: -kv[1]["gained_over_genomad"])[:8]})
    print(f"[done] wrote {out}")


if __name__ == "__main__":
    main()
