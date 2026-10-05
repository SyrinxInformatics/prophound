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
"""Head-to-head of the fragmented-input fix candidates, for the freeze decision.

Compares, in every mode and on both corpus sets (42 tune, 18 check), the shipped
pipeline against edge snap only (ES) and edge snap plus the whole-contig call
(ES+WC), applied as one fix to all modes (user, 2026-10-05). Same inputs and
same post-processing code as scripts/g5dev_eval.py, so its numbers match the sweep.

Reported per set x mode:
  - MAG covered50 (three shred seeds pooled), complete, draft; calls/genome; false/Mb
  - paired genome-bootstrap 95% interval of the MAG covered50 difference
    (ES - shipped, ES+WC - ES), genomes resampled, seeds pooled inside a genome
  - MAG covered50 by how many fragments an element is cut into (1, 2, 3+), the
    split where G5's post-hoc look found the deficit
  - WC accounting on MAG: WC calls, how many touch an element, elements gained by WC
  - what covered50 and false/Mb cannot see: bp called, bp called outside elements,
    median IoU of detected elements (best single call), and for ES the calls it
    moves and the bp it adds, split into element and non-element bp
Phage-free genomes are not repeated here; g5dev_eval.py's check covers them.

Example:
  python scripts/g5dev_compare.py --run_root analysis/g5dev/prophound --lengths analysis/g5dev/contig_lengths.tsv \\
      --set tune=scripts/eval_sizestrat_elements.csv --set check=scripts/eval_confirm18_elements.csv \\
      --out analysis/g5dev/compare.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mode_recheck as mr  # noqa: E402
from g5_score import cov_frac, split_frag  # noqa: E402
from g5dev_eval import apply_variant, load, to_genome, window_table  # noqa: E402

MAGS = ["mag_s0", "mag_s1", "mag_s2"]
CONDS = ["complete", *MAGS, "draft_s0"]
VARIANTS = {"shipped": (0, None, None), "ES": (10000, None, None), "ES+WC": (10000, 0.8, 0.7)}
MODES = ["default", "exploratory", "conservative"]
B = 10_000


def n_fragments(T: pd.DataFrame, lc: pd.DataFrame) -> np.ndarray:
    """Per element, the number of fragments of this shred that overlap it."""
    so = lc.contig_id.map(split_frag)
    fr = pd.DataFrame({"assembly": lc.genome.to_numpy(), "sequence": [s for s, _ in so],
                       "start": [o for _, o in so]})
    fr["end"] = fr.start + lc.length.to_numpy()
    by = {k: g[["start", "end"]].to_numpy() for k, g in fr.groupby(["assembly", "sequence"])}
    out = np.zeros(len(T), int)
    for i, r in enumerate(T.itertuples()):
        f = by.get((r.assembly, r.sequence))
        if f is not None:
            out[i] = int(((f[:, 1] > r.start) & (f[:, 0] < r.end)).sum())
    return out


def off_element_bp(calls: pd.DataFrame, T: pd.DataFrame) -> np.ndarray:
    """Per call, bp not inside any element on its sequence."""
    by = {k: g[["start", "end"]].to_numpy() for k, g in T.groupby(["assembly", "sequence"])}
    out = np.zeros(len(calls), np.int64)
    for i, c in enumerate(calls.itertuples()):
        e = by.get((c.assembly, c.sequence))
        inside = 0
        if e is not None:
            iv = sorted((max(s, c.start), min(t, c.end)) for s, t in e if t > c.start and s < c.end)
            ce = c.start
            for s, t in iv:
                s = max(s, ce)
                if t > s:
                    inside += t - s
                    ce = t
        out[i] = (c.end - c.start) - inside
    return out


def best_iou(T: pd.DataFrame, calls: pd.DataFrame) -> np.ndarray:
    """Per element, the IoU of its best-overlapping single call (0 if none)."""
    by = {k: g[["start", "end"]].to_numpy() for k, g in calls.groupby(["assembly", "sequence"])}
    out = np.zeros(len(T))
    for i, r in enumerate(T.itertuples()):
        c = by.get((r.assembly, r.sequence))
        if c is None:
            continue
        inter = np.clip(np.minimum(c[:, 1], r.end) - np.maximum(c[:, 0], r.start), 0, None)
        union = np.maximum(c[:, 1], r.end) - np.minimum(c[:, 0], r.start)
        out[i] = float((inter / union).max()) if len(c) else 0.0
    return out


def boot_diff(hits_a: np.ndarray, hits_b: np.ndarray, genome: np.ndarray, rng) -> tuple[float, float, float]:
    """Paired genome bootstrap of mean(b) - mean(a); rows are element x seed."""
    gs = np.unique(genome)
    idx = [np.flatnonzero(genome == g) for g in gs]
    a_s = np.array([hits_a[i].sum() for i in idx], float)
    b_s = np.array([hits_b[i].sum() for i in idx], float)
    n = np.array([len(i) for i in idx], float)
    draw = rng.integers(0, len(gs), (B, len(gs)))
    d = (b_s[draw].sum(1) - a_s[draw].sum(1)) / n[draw].sum(1)
    return float((b_s.sum() - a_s.sum()) / n.sum()), float(np.quantile(d, 0.025)), float(np.quantile(d, 0.975))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run_root", required=True)
    ap.add_argument("--lengths", required=True)
    ap.add_argument("--set", action="append", required=True, metavar="NAME=ELEMENTS_CSV")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rng = np.random.default_rng(0)
    Lt = pd.read_csv(a.lengths, sep="\t", dtype={"contig_id": str})
    res = {}
    for sname, path in (x.split("=", 1) for x in a.set):
        T = pd.read_csv(path).rename(columns={"accession": "assembly", "contig": "sequence"})
        T = T[["assembly", "sequence", "start", "end", "stratum"]].reset_index(drop=True)
        genomes = sorted(T.assembly.unique())
        bp = Lt[(Lt.condition == "complete") & Lt.genome.isin(genomes)].groupby("genome").length.sum()
        el_bp = (T.end - T.start).groupby(T.assembly).sum()
        free_mb = float((bp - el_bp.reindex(bp.index).fillna(0)).sum()) / 1e6
        hits, cells, nfrag, wc_acc = {}, {}, {}, {}
        for cond in CONDS:
            ens, seeds = load(Path(a.run_root), cond, genomes)
            lc = Lt[(Lt.condition == cond) & Lt.genome.isin(genomes)]
            lengths = {(g, c): int(n) for g, c, n in zip(lc.genome, lc.contig_id, lc.length)}
            if cond in MAGS:
                nfrag[cond] = n_fragments(T, lc)
            win = window_table(ens, seeds)
            cache = {}
            for mode in MODES:
                base = mr.shipped_calls(ens, seeds, mode)
                for vname, (d, t, f) in VARIANTS.items():
                    reg, n_wc = apply_variant(base, win, lengths, mode, d, t, f, cache)
                    calls = to_genome(reg)
                    h = cov_frac(T, calls) >= 0.5
                    on = np.array([((T.assembly == c.assembly) & (T.sequence == c.sequence) & (T.end > c.start)
                                    & (T.start < c.end)).any() for c in calls.itertuples()], bool)
                    hits[(cond, mode, vname)] = h
                    off = off_element_bp(calls, T)
                    clen = (calls.end - calls.start).to_numpy()
                    iou = best_iou(T, calls)
                    cells[(cond, mode, vname)] = {"covered50": float(h.mean()), "calls_per_genome": len(calls) / len(genomes),
                                                  "false_per_mb": float((~on).sum()) / free_mb, "wc_calls": n_wc,
                                                  "bp_called": int(clen.sum()), "bp_off_element": int(off.sum()),
                                                  "iou": iou[h]}
                    if vname == "ES":
                        b = to_genome(base)
                        moved = ((b.start.to_numpy() != calls.start.to_numpy()) | (b.end.to_numpy() != calls.end.to_numpy()))
                        boff = off_element_bp(b, T)
                        added = clen - (b.end - b.start).to_numpy()
                        cells[(cond, mode, vname)]["snap"] = {
                            "calls_moved": int(moved.sum()), "calls": len(calls),
                            "moved_false": int((moved & ~on).sum()),
                            "bp_added": int(added.sum()), "bp_added_off_element": int((off - boff).sum()),
                            "median_bp_added_per_moved": float(np.median(added[moved])) if moved.any() else 0.0}
                    if vname == "ES+WC" and cond in MAGS:
                        is_wc = reg["wc"].eq(True).to_numpy() if "wc" in reg else np.zeros(len(reg), bool)
                        wc_acc.setdefault(mode, {"wc_calls": 0, "wc_on_element": 0})
                        wc_acc[mode]["wc_calls"] += int(is_wc.sum())
                        wc_acc[mode]["wc_on_element"] += int(on[is_wc].sum())
        out = res[sname] = {"n_genomes": len(genomes), "n_elements": len(T), "modes": {}}
        genome_rows = np.concatenate([T.assembly.to_numpy()] * len(MAGS))
        nf = np.concatenate([nfrag[c] for c in MAGS])
        for mode in MODES:
            m = out["modes"][mode] = {}
            mag = {v: np.concatenate([hits[(c, mode, v)] for c in MAGS]) for v in VARIANTS}
            for v in VARIANTS:
                m[v] = {"complete": cells[("complete", mode, v)]["covered50"],
                        "mag": float(mag[v].mean()),
                        "draft": cells[("draft_s0", mode, v)]["covered50"],
                        "mag_calls_per_genome": float(np.mean([cells[(c, mode, v)]["calls_per_genome"] for c in MAGS])),
                        "mag_false_per_mb": float(np.mean([cells[(c, mode, v)]["false_per_mb"] for c in MAGS])),
                        "complete_false_per_mb": cells[("complete", mode, v)]["false_per_mb"],
                        "mag_bp_called": int(sum(cells[(c, mode, v)]["bp_called"] for c in MAGS)),
                        "mag_bp_off_element": int(sum(cells[(c, mode, v)]["bp_off_element"] for c in MAGS)),
                        "mag_median_iou_detected": float(np.median(np.concatenate([cells[(c, mode, v)]["iou"] for c in MAGS]))),
                        "complete_median_iou_detected": float(np.median(cells[("complete", mode, v)]["iou"])),
                        "complete_bp_called": cells[("complete", mode, v)]["bp_called"],
                        "by_fragments": {k: float(mag[v][sel].mean()) for k, sel in
                                         (("1", nf == 1), ("2", nf == 2), ("3+", nf >= 3)) if sel.any()}}
            m["n_by_fragments"] = {"1": int((nf == 1).sum()), "2": int((nf == 2).sum()), "3+": int((nf >= 3).sum()),
                                   "0 (lost)": int((nf == 0).sum())}
            m["diff_ES_vs_shipped"] = boot_diff(mag["shipped"], mag["ES"], genome_rows, rng)
            m["diff_ESWC_vs_ES"] = boot_diff(mag["ES"], mag["ES+WC"], genome_rows, rng)
            m["elements_gained_by_WC"] = int((mag["ES+WC"] & ~mag["ES"]).sum())
            m["elements_lost_by_WC"] = int((mag["ES"] & ~mag["ES+WC"]).sum())
            m["wc_accounting"] = wc_acc.get(mode, {})
            m["snap"] = {c: cells[(c, mode, "ES")]["snap"] for c in CONDS}
    Path(a.out).write_text(json.dumps(res, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o)))

    for sname, s in res.items():
        print(f"\n=== {sname}: {s['n_genomes']} genomes, {s['n_elements']} elements (x3 MAG seeds)")
        for mode, m in s["modes"].items():
            print(f"[{mode}]  elements by fragment count: {m['n_by_fragments']}")
            for v in VARIANTS:
                r = m[v]
                bf = " ".join(f"{k}:{x:.3f}" for k, x in r["by_fragments"].items())
                print(f"  {v:8s} complete {r['complete']:.3f}  MAG {r['mag']:.3f}  draft {r['draft']:.3f}  "
                      f"calls/g {r['mag_calls_per_genome']:5.1f}  false/Mb {r['mag_false_per_mb']:.2f} "
                      f"(complete {r['complete_false_per_mb']:.2f})  by frags [{bf}]")
            for k in ("diff_ES_vs_shipped", "diff_ESWC_vs_ES"):
                d, lo, hi = m[k]
                print(f"  {k:20s} {d:+.3f} ({lo:+.3f} to {hi:+.3f})")
            for v in VARIANTS:
                r = m[v]
                print(f"  {v:8s} MAG bp called {r['mag_bp_called']/1e6:7.2f} Mb, off-element {r['mag_bp_off_element']/1e6:6.2f} Mb; "
                      f"median IoU detected MAG {r['mag_median_iou_detected']:.3f}, complete {r['complete_median_iou_detected']:.3f}")
            for c, sn in m["snap"].items():
                print(f"  ES snap {c:9s} moved {sn['calls_moved']}/{sn['calls']} calls ({sn['moved_false']} false); "
                      f"+{sn['bp_added']/1e3:.0f} kb, of which off-element {sn['bp_added_off_element']/1e3:.0f} kb; "
                      f"median +{sn['median_bp_added_per_moved']/1e3:.1f} kb per moved call")
            print(f"  WC: gained {m['elements_gained_by_WC']}, lost {m['elements_lost_by_WC']} element-seeds; "
                  f"{m['wc_accounting']}")
    print(f"\n[done] wrote {a.out}")


if __name__ == "__main__":
    main()
