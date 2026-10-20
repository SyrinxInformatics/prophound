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
"""Development sweep for the fragmented-input fix, on corpus shreds only.

docs/claims-and-evidence.md "G5b preregistration": the fix is built here, on the
42 size-stratified eval genomes (tuning) and checked on the 18 confirmation
genomes, with geNomad truth. LAMBDA and Philympics are never read. Only
prophound's own retention is measured: covered50 on shreds over covered50 on the
same genomes complete.

Every variant post-processes the shipped calls (mode_recheck.shipped_calls: the
mode's thresholds, coverage builder, merging, short-contig scales, seed vote)
with one or both of:
  edge snap  ES(d)     a call within d bp of its fragment's end is extended to that
                       end -- an element cut by a breakpoint runs off the fragment
  whole      WC(t, f)  a scanned fragment where >= f of its windows score >= t is
  contig               called end to end; in a voting mode, at least the mode's
                       min seed support of checkpoints must pass on their own
                       scores (prob_ckpt<k>), so WC obeys the same vote as calls
WC calls replace the fragment's other calls. Fragments under 3 kb are never
scanned, so WC cannot reach them either.

Selection rule, fixed before the data exist (same "never exceed" principle as
the short-contig calibration, predict_genome.SHORT_CONTIG_THRESHOLDS): on the 42,
the variant with the highest Default MAG retention whose MAG false calls per Mb
of element-free sequence do not exceed the shipped pipeline's on the complete
genomes, in that mode. The 18 then check it; nothing is chosen on them. Ties go
to the simpler variant (fewer active rules, then smaller d, then larger t and f).

Complete-genome regression check (reported, not part of selection): the
selected variant's calls on the 20 geNomad-negative test-split genomes
(--phagefree), against the shipped pipeline's, and on the `complete` condition
of both sets. By design neither rule should change much on intact chromosomes;
this measures it.

Inputs: tracks <run_root>/<condition>/<genome>/{windows.tsv, run_metadata.json}
written with --write-seed-scores (g5dev_prophound.slurm), and a fragment-length
table (condition, genome, contig_id, length) covering every contig of every
condition, so whole-contig calls get exact ends.

Example:
  python scripts/g5dev_eval.py --run_root analysis/g5dev/prophound --lengths analysis/g5dev/contig_lengths.tsv \\
      --set tune=scripts/eval_sizestrat_elements.csv --set check=scripts/eval_confirm18_elements.csv \\
      --phagefree analysis/phagefree --out analysis/g5dev/sweep.json
"""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import evalkit as ek
import numpy as np
import pandas as pd
from evalkit import cov_frac, split_frag

CONDS = ["complete", "mag_s0", "mag_s1", "mag_s2", "draft_s0"]
ES_GRID = [0, 2000, 5000, 10000]
WC_T = [None, 0.5, 0.8, "mode"]                 # None = WC off; "mode" = the mode's seed threshold
WC_F = [0.5, 0.7, 0.9]


def load(run_root: Path, cond: str, genomes: list) -> tuple[dict, list]:
    """Ensemble and per-seed tracks for one condition, restricted to `genomes`."""
    d = run_root / cond
    keep = lambda t: {g: v for g, v in t.items() if g in genomes}  # noqa: E731
    ens = keep(ek.load_tracks(str(d)))
    seeds = [keep(ek.load_tracks(str(d), f"prob_ckpt{k}")) for k in (1, 2, 3)]
    missing = set(genomes) - set(ens)
    if missing:
        raise SystemExit(f"{d}: no tracks for {len(missing)} genomes, e.g. {sorted(missing)[:3]}")
    return ens, seeds


def window_table(ens: dict, seeds: list) -> pd.DataFrame:
    """All windows of a condition with ensemble and per-seed scores side by side."""
    parts = []
    for g, (w, _) in ens.items():
        x = w[["genome_id", "contig_id", "start", "end", "prob"]].copy()
        for k, s in enumerate(seeds, 1):
            x[f"p{k}"] = s[g][0]["prob"].to_numpy()
        parts.append(x)
    return pd.concat(parts, ignore_index=True)


def edge_snap(reg: pd.DataFrame, lengths: dict, d: int) -> pd.DataFrame:
    if d <= 0 or reg.empty:
        return reg
    reg = reg.copy()
    L = reg.apply(lambda r: lengths[(r.genome_id, str(r.contig_id))], axis=1).to_numpy()
    reg.loc[reg.start.to_numpy() <= d, "start"] = 0
    snap = reg.end.to_numpy() >= L - d
    reg.loc[snap, "end"] = L[snap]
    return reg


def whole_contig(win: pd.DataFrame, lengths: dict, t: float, f: float, k: int) -> pd.DataFrame:
    """Fragments called end to end: >= f of windows >= t, by >= k checkpoints (k=1: the ensemble)."""
    if k <= 1:
        ok = win.groupby(["genome_id", "contig_id"]).prob.apply(lambda p: (p >= t).mean() >= f)
    else:
        per = [win.groupby(["genome_id", "contig_id"])[f"p{j}"].apply(lambda p: (p >= t).mean() >= f) for j in (1, 2, 3)]
        ok = sum(p.astype(int) for p in per) >= k
    hit = ok[ok].index
    mp = win.groupby(["genome_id", "contig_id"]).prob.max()
    return pd.DataFrame([{"genome_id": g, "contig_id": c, "start": 0, "end": lengths[(g, str(c))],
                          "max_prob": float(mp[(g, c)]), "wc": True} for g, c in hit],
                        columns=["genome_id", "contig_id", "start", "end", "max_prob", "wc"])


def to_genome(reg: pd.DataFrame) -> pd.DataFrame:
    """Fragment coordinates -> chromosome coordinates, in g5_score's call layout."""
    so = reg.contig_id.astype(str).map(split_frag)
    return pd.DataFrame({"assembly": reg.genome_id.to_numpy(), "sequence": [s for s, _ in so],
                         "start": reg.start.to_numpy().astype(int) + np.array([o for _, o in so], int),
                         "end": reg.end.to_numpy().astype(int) + np.array([o for _, o in so], int)})


def apply_variant(base: pd.DataFrame, win: pd.DataFrame, lengths: dict, mode: str, d: int, t, f,
                  cache: dict) -> tuple[pd.DataFrame, int]:
    """Shipped calls -> calls under variant (d, t, f); returns the calls and the number of WC calls."""
    reg = edge_snap(base, lengths, d)
    if t is None:
        return reg, 0
    k = ek.MODE_MIN_SEED_SUPPORT[mode]
    tt = ek.MODE_THRESHOLDS[mode] if t == "mode" else t
    if (tt, f, k) not in cache:
        cache[(tt, f, k)] = whole_contig(win, lengths, tt, f, k)
    wc = cache[(tt, f, k)]
    if len(wc):
        key = set(zip(wc.genome_id, wc.contig_id.astype(str)))
        reg = reg[[(g, str(c)) not in key for g, c in zip(reg.genome_id, reg.contig_id)]]
        reg = pd.concat([reg, wc], ignore_index=True)
    return reg, len(wc)


def phagefree_check(root: Path, summary_csv: str, modes: list, variants: list) -> dict:
    """Calls on the 20 geNomad-negative test-split genomes (complete), per mode and variant.

    The tracks keep each checkpoint in its own run dir (<root>/{ens,r1,r2,r3}/<genome>/,
    from confirm_eval.slurm), joined here row by row. Exact contig lengths are not
    available locally, so a contig's length is taken as its furthest window end, which
    is within one stride (1 kb) of the true end; only edge snapping reads it.
    """
    ens = ek.load_tracks(str(root / "ens"))
    seeds = []
    for r in ("r1", "r2", "r3"):
        t = ek.load_tracks(str(root / r))
        for g, (w, _) in t.items():
            if not (w[["contig_id", "start", "end"]].to_numpy() == ens[g][0][["contig_id", "start", "end"]].to_numpy()).all():
                raise SystemExit(f"{root}/{r}/{g}: windows do not line up with the ensemble run")
        seeds.append(t)
    lengths = {(g, str(c)): int(e) for g, (w, _) in ens.items() for c, e in w.groupby("contig_id").end.max().items()}
    bp = pd.read_csv(summary_csv).set_index("accession").total_bp.reindex(list(ens))
    if bp.isna().any():
        raise SystemExit(f"{summary_csv} lacks total_bp for {list(bp[bp.isna()].index)[:3]}")
    win = window_table(ens, seeds)
    out, cache = {}, {}
    for mode in modes:
        base = ek.shipped_calls(ens, seeds, mode)
        for d, t, f in variants:
            reg, n_wc = apply_variant(base, win, lengths, mode, d, t, f, cache)
            out[f"{mode}|es={d}|wc={t},{f}"] = {
                "calls_per_genome": len(reg) / len(ens), "calls_per_mb": len(reg) / (bp.sum() / 1e6),
                "frac_bp_called": float((reg.end - reg.start).sum() / bp.sum()), "wc_calls": n_wc}
    return {"n_genomes": len(ens), "cells": out}


def metrics(T: pd.DataFrame, calls: pd.DataFrame, genome_bp: pd.Series) -> dict:
    f = cov_frac(T, calls)
    hit = f >= 0.5
    on = np.zeros(len(calls), bool)
    for i, c in enumerate(calls.itertuples()):
        on[i] = ((T.assembly == c.assembly) & (T.sequence == c.sequence) & (T.end > c.start) & (T.start < c.end)).any()
    el_bp = (T.end - T.start).groupby(T.assembly).sum()
    free_mb = float((genome_bp - el_bp.reindex(genome_bp.index).fillna(0)).sum()) / 1e6
    per_stratum = {s: float(hit[(T.stratum == s).to_numpy()].mean()) for s in sorted(T.stratum.unique())}
    return {"covered50": float(hit.mean()), "k50": int(hit.sum()), "n": int(len(T)),
            "calls_per_genome": len(calls) / len(genome_bp), "false_per_mb": float((~on).sum()) / free_mb,
            "per_stratum": per_stratum, "_hit": hit}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run_root", required=True)
    ap.add_argument("--lengths", required=True, help="condition, genome, contig_id, length (tsv)")
    ap.add_argument("--set", action="append", required=True, metavar="NAME=ELEMENTS_CSV",
                    help="tune=... (selection) and check=... (confirmation only)")
    ap.add_argument("--modes", default="default,exploratory,conservative")
    ap.add_argument("--phagefree", default=None, help="dir with ens/ r1/ r2/ r3/ tracks of geNomad-negative genomes")
    ap.add_argument("--summary", default="summary.csv", help="accession,total_bp (for --phagefree)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    Lt = pd.read_csv(a.lengths, sep="\t", dtype={"contig_id": str})
    sets = dict(x.split("=", 1) for x in a.set)
    if "tune" not in sets:
        raise SystemExit("--set tune=... is required (the selection sample)")
    modes = a.modes.split(",")
    variants = [(d, t, f) for d, t in itertools.product(ES_GRID, WC_T) for f in (WC_F if t is not None else [None])]

    res = {"variants": [{"es": d, "wc_t": t, "wc_f": f} for d, t, f in variants], "sets": {}}
    for sname, path in sets.items():
        ev = pd.read_csv(path)
        T = ev.rename(columns={"accession": "assembly", "contig": "sequence"})[["assembly", "sequence", "start", "end", "stratum"]]
        T = T.reset_index(drop=True)
        genomes = sorted(T.assembly.unique())
        out = res["sets"][sname] = {"n_genomes": len(genomes), "n_elements": len(T), "cells": {}}
        hits = {}
        for cond in CONDS:
            ens, seeds = load(Path(a.run_root), cond, genomes)
            lc = Lt[(Lt.condition == cond) & Lt.genome.isin(genomes)]
            lengths = {(g, c): int(n) for g, c, n in zip(lc.genome, lc.contig_id, lc.length)}
            genome_bp = lc.groupby("genome").length.sum() if cond == "complete" else None
            if cond == "complete":
                res.setdefault("_bp", {})[sname] = genome_bp
            genome_bp = res["_bp"][sname]
            win = window_table(ens, seeds)
            wc_cache = {}
            for mode in modes:
                base = ek.shipped_calls(ens, seeds, mode)
                for d, t, f in variants:
                    reg, n_wc = apply_variant(base, win, lengths, mode, d, t, f, wc_cache)
                    m = metrics(T, to_genome(reg), genome_bp)
                    hits[(cond, mode, d, t, f)] = m.pop("_hit")
                    out["cells"][f"{cond}|{mode}|es={d}|wc={t},{f}"] = {**m, "wc_calls": n_wc}
        # retention: MAG pooled over seeds and draft, each over the same variant on complete genomes
        ret = {}
        for mode in modes:
            for d, t, f in variants:
                c = hits[("complete", mode, d, t, f)].mean()
                key = f"{mode}|es={d}|wc={t},{f}"
                mag = np.mean([hits[(s, mode, d, t, f)].mean() for s in ("mag_s0", "mag_s1", "mag_s2")])
                ret[key] = {"complete": float(c), "mag": float(mag), "draft": float(hits[("draft_s0", mode, d, t, f)].mean()),
                            "mag_retention": float(mag / c) if c else None,
                            "mag_false_per_mb": float(np.mean([out["cells"][f"{s}|{mode}|es={d}|wc={t},{f}"]["false_per_mb"]
                                                               for s in ("mag_s0", "mag_s1", "mag_s2")])),
                            "mag_calls_per_genome": float(np.mean([out["cells"][f"{s}|{mode}|es={d}|wc={t},{f}"]["calls_per_genome"]
                                                                   for s in ("mag_s0", "mag_s1", "mag_s2")]))}
        out["retention"] = ret
    res.pop("_bp", None)
    if a.phagefree:
        res["phagefree"] = phagefree_check(Path(a.phagefree), a.summary, modes, variants)

    # selection on `tune` only, per the rule in the docstring
    tune = res["sets"]["tune"]
    sel = {}
    for mode in modes:
        cap = tune["cells"][f"complete|{mode}|es=0|wc=None,None"]["false_per_mb"]
        ok = [(k, v) for k, v in tune["retention"].items() if k.startswith(mode + "|") and v["mag_false_per_mb"] <= cap]

        def simplicity(k):
            d = int(k.split("es=")[1].split("|")[0]); t, f = k.split("wc=")[1].split(",")
            tv = 0.0 if t == "None" else 1.0 if t == "mode" else float(t)       # "mode" (>= 0.95) is the strictest
            return ((d > 0) + (t != "None"), d, -tv, -(0 if f == "None" else float(f)))
        if not ok:
            sel[mode] = {"variant": None, "false_per_mb_cap": cap, "note": "no variant within the cap"}
            continue
        best = max(ok, key=lambda kv: (round(kv[1]["mag_retention"], 6), tuple(-x for x in simplicity(kv[0]))))
        sel[mode] = {"variant": best[0], "false_per_mb_cap": cap, **best[1],
                     "baseline": tune["retention"][f"{mode}|es=0|wc=None,None"],
                     "check": {s: res["sets"][s]["retention"][best[0]] | {
                         "baseline": res["sets"][s]["retention"][f"{mode}|es=0|wc=None,None"]}
                         for s in res["sets"] if s != "tune"}}
    res["selected"] = sel

    o = Path(a.out)
    o.parent.mkdir(parents=True, exist_ok=True)
    o.write_text(json.dumps(res, indent=2, default=str) + "\n")
    rows = [{"set": s, "variant": k, **{x: y for x, y in v.items()}} for s, d in res["sets"].items() for k, v in d["retention"].items()]
    pd.DataFrame(rows).to_csv(o.with_suffix(".csv"), index=False)
    for mode, v in sel.items():
        if v["variant"] is None:
            print(f"[{mode}] no variant within the false/Mb cap {v['false_per_mb_cap']:.2f}")
            continue
        b = v["baseline"]
        print(f"[{mode}] baseline MAG retention {b['mag_retention']:.3f} (calls/genome {b['mag_calls_per_genome']:.1f}, "
              f"false/Mb {b['mag_false_per_mb']:.2f}, cap {v['false_per_mb_cap']:.2f})")
        print(f"   selected {v['variant']}: retention {v['mag_retention']:.3f}, complete {v['complete']:.3f}, "
              f"MAG {v['mag']:.3f}, draft {v['draft']:.3f}, calls/genome {v['mag_calls_per_genome']:.1f}, "
              f"false/Mb {v['mag_false_per_mb']:.2f}")
        if "phagefree" in res:
            pf = res["phagefree"]["cells"]
            b0, b1 = pf[f"{mode}|es=0|wc=None,None"], pf[v["variant"]]
            print(f"   phage-free complete genomes (n={res['phagefree']['n_genomes']}): calls/genome "
                  f"{b1['calls_per_genome']:.2f} vs baseline {b0['calls_per_genome']:.2f}; WC calls {b1['wc_calls']}")
        for s, c in v["check"].items():
            print(f"   check on {s}: retention {c['mag_retention']:.3f} vs baseline {c['baseline']['mag_retention']:.3f}; "
                  f"false/Mb {c['mag_false_per_mb']:.2f} vs {c['baseline']['mag_false_per_mb']:.2f}")
    print(f"[done] wrote {o}")


if __name__ == "__main__":
    main()
