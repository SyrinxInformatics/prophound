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
"""Offline sweep of the fine refinement pass across all four downstream problems.

docs/downstream-postprocessing.md, docs/adaptive-multiscale-detection.md. Each
run directory holds `<genome>/windows.tsv` and `<genome>/fine_windows.tsv` from
`predict_genome --refine-window-len W` (the score tracks); pass B's parameters
are cheap, so they are swept here without re-scoring. Coarse calls are rebuilt
from windows.tsv using the run's own metadata, so the baseline row is the
shipped operating point.

For every (refine_threshold, edge_alpha, split_min_gap_bp, split_max_prob)
variant, overall and per size stratum:
  - recall:         covered and matched recall   (nothing dropped)
  - conflation:     conflated calls, absorbed elements
  - fragmentation:  elements hit by >1 call      (gap-filling not undone)
  - truncation:     median undershoot bp          (edges inside the element)
  - overextension:  median overshoot bp           (edges outside the element)
  - boundary:       median IoU; calls/genome as the cost
A variant is only interesting if it improves one of these without paying for it
in another; the table is built so that trade is visible, not averaged away.

Example:
  python scripts/refine_sweep.py --elements scripts/eval_sizestrat_elements.csv \\
      --run r1_3k=analysis/multiscale/r1_3000 --out analysis/multiscale/sweep.json
"""

from __future__ import annotations

import argparse
import itertools
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # data_refactored/ layout: src/ is a sibling
from compare_variants import score_variant  # noqa: E402

try:  # repo layout
    from mge_detect.models.prophage.src.multiscale import refine_regions
    from mge_detect.models.prophage.src.predict_genome import build_region_table_coverage
except ImportError:  # HPCC data_refactored/ layout
    from src.multiscale import refine_regions
    from src.predict_genome import build_region_table_coverage


def edge_errors(ev: pd.DataFrame, reg: pd.DataFrame) -> dict:
    """Median overshoot / undershoot (bp) of the best-IoU call per covered element."""
    over, under = [], []
    for r in ev.itertuples():
        g = reg[(reg.genome_id == r.accession) & (reg.contig_id == r.contig)
                & (reg.end > r.start) & (reg.start < r.end)]
        if g.empty:
            continue
        inter = np.minimum(g.end, r.end) - np.maximum(g.start, r.start)
        union = np.maximum(g.end, r.end) - np.minimum(g.start, r.start)
        c = g.iloc[int(np.argmax((inter / union).to_numpy()))]
        over.append(max(0, r.start - c.start) + max(0, c.end - r.end))
        under.append(max(0, c.start - r.start) + max(0, r.end - c.end))
    return {"median_overshoot_bp": float(np.median(over)) if over else None,
            "median_undershoot_bp": float(np.median(under)) if under else None}


KEEP = ("covered_recall", "covered50_recall", "matched_recall", "median_iou", "conflated_calls",
        "absorbed_elements", "fragmented_elements", "calls_per_genome", "n_elements")


def covered50(ev: pd.DataFrame, reg: pd.DataFrame) -> float:
    """Fraction of elements with >= 50% of their bp under some call.

    Covered recall (>= 1 bp) credits grazes -- a neighbouring call's overshoot
    clipping an element's edge. Both losses 817370's refinement showed were
    grazes (docs/downstream-postprocessing.md "Per-element dig"), so every
    variant is also scored on this graze-proof version.
    """
    ok = 0
    for r in ev.itertuples():
        g = reg[(reg.genome_id == r.accession) & (reg.contig_id == r.contig)
                & (reg.end > r.start) & (reg.start < r.end)].sort_values("start")
        bp, cur_s, cur_e = 0, None, None
        for c in g.itertuples():  # union of overlapping calls, clipped to the element
            s, e = max(int(c.start), r.start), min(int(c.end), r.end)
            if cur_e is None or s > cur_e:
                bp += 0 if cur_e is None else cur_e - cur_s
                cur_s, cur_e = s, e
            else:
                cur_e = max(cur_e, e)
        bp += 0 if cur_e is None else cur_e - cur_s
        ok += bp >= 0.5 * (r.end - r.start)
    return ok / max(1, len(ev))


def score(ev: pd.DataFrame, reg: pd.DataFrame) -> dict:
    m = score_variant(ev, reg)
    m["covered50_recall"] = covered50(ev, reg)
    out = {k: m[k] for k in KEEP}
    out.update(edge_errors(ev, reg))
    out["by_stratum"] = {}
    for st, sub in ev.groupby("stratum"):
        sm = score_variant(sub, reg)
        sm["covered50_recall"] = covered50(sub, reg)
        row = {k: sm[k] for k in KEEP if k != "calls_per_genome"}
        row.update(edge_errors(sub, reg))
        out["by_stratum"][str(st)] = row
    return out


def load_run(d: str):
    """(coarse regions, fine windows) per genome, rebuilt from the run's metadata."""
    coarse, fine, meta0 = [], [], None
    for gdir in sorted(Path(d).iterdir()):
        w, f, m = gdir / "windows.tsv", gdir / "fine_windows.tsv", gdir / "run_metadata.json"
        if not (w.exists() and m.exists()):
            continue
        meta = json.loads(m.read_text())
        meta0 = meta0 or meta
        wdf = pd.read_csv(w, sep="\t")
        fm = meta.get("fragment_merge", {})
        _, reg = build_region_table_coverage(
            wdf, seed_threshold=meta["seed_threshold"], extend_threshold=meta["extend_threshold"],
            min_coverage_frac=meta["coverage_frac"], max_gap_bp=fm.get("max_gap_bp", 0),
            bridge_prob=fm.get("bridge_prob", 0.0))
        coarse.append(reg)
        if f.exists():
            fine.append(pd.read_csv(f, sep="\t", usecols=["genome_id", "contig_id", "start", "end", "prob"]))
    if not coarse:
        raise SystemExit(f"no <genome>/windows.tsv + run_metadata.json under {d}")
    fine_df = pd.concat(fine, ignore_index=True) if fine else pd.DataFrame(
        columns=["genome_id", "contig_id", "start", "end", "prob"])
    return pd.concat(coarse, ignore_index=True), fine_df, meta0, len(coarse), len(fine)


def fine_regions_by_genome(fine: pd.DataFrame, thr: float, frac: float) -> dict:
    """genome -> (fine windows, fine regions) at one threshold; reused across variants."""
    out = {}
    for g, fg in fine.groupby("genome_id", sort=False):
        out[g] = build_region_table_coverage(fg.assign(fasta_path="x"), seed_threshold=thr,
                                             extend_threshold=thr, min_coverage_frac=frac, max_gap_bp=0)
    return out


def refine_all(coarse: pd.DataFrame, fine_by_genome: dict, alpha: float, gap: int,
               maxp: float, ramp=None) -> pd.DataFrame:
    out = []
    for g, cg in coarse.groupby("genome_id", sort=False):
        if g not in fine_by_genome:
            out.append(cg)
            continue
        fg, fr = fine_by_genome[g]
        out.append(refine_regions(cg, fg, fr, edge_alpha=alpha, split_min_gap_bp=gap,
                                  split_max_prob=maxp, alpha_ramp_bp=ramp))
    return pd.concat(out, ignore_index=True)


def floats(s: str):
    return [float(x) for x in s.split(",") if x]


def ramps(s: str):
    """'none,10000:60000' -> [None, (10000, 60000)]"""
    out = []
    for x in s.split(","):
        if x == "none":
            out.append(None)
        elif x:
            lo, hi = x.split(":")
            out.append((int(lo), int(hi)))
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--elements", required=True)
    ap.add_argument("--run", action="append", required=True, metavar="LABEL=DIR")
    ap.add_argument("--refine_thresholds", default="0.95,0.99,0.995")
    ap.add_argument("--alphas", default="1,0.75,0.5,0.25,0")
    ap.add_argument("--split_gaps", default="1000,2000,4000")
    ap.add_argument("--split_max_probs", default="0.05,0.1,0.2")
    ap.add_argument("--ramps", default="none",
                    help="alpha ramps over fine-core length, e.g. none,10000:60000 (alpha = --alphas value at HI)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    ev = pd.read_csv(args.elements)
    res = {"n_elements": int(len(ev)), "runs": {}}
    flat = []
    for spec in args.run:
        label, _, d = spec.partition("=")
        coarse, fine, meta, n_g, n_f = load_run(d)
        if n_f < n_g:
            print(f"[warn] {label}: fine_windows.tsv in {n_f}/{n_g} genomes", flush=True)
        frac = float(meta["coverage_frac"])
        run = {"n_genomes": n_g, "fine_window_len": (meta.get("refine") or {}).get("window_len"),
               "baseline": score(ev, coarse), "variants": {}}
        b = run["baseline"]
        print(f"[{label}] baseline covered {b['covered_recall']:.3f} matched {b['matched_recall']:.3f} "
              f"absorbed {b['absorbed_elements']} over {b['median_overshoot_bp']} under "
              f"{b['median_undershoot_bp']}", flush=True)
        flat.append({"run": label, "variant": "baseline", **{k: b[k] for k in KEEP},
                     "overshoot": b["median_overshoot_bp"], "undershoot": b["median_undershoot_bp"]})
        fine_cache = {thr: fine_regions_by_genome(fine, thr, frac) for thr in floats(args.refine_thresholds)}
        for thr, alpha, gap, maxp, ramp in itertools.product(
                floats(args.refine_thresholds), floats(args.alphas),
                floats(args.split_gaps), floats(args.split_max_probs), ramps(args.ramps)):
            key = f"t{thr:g}_a{alpha:g}_g{int(gap)}_p{maxp:g}" + (f"_r{ramp[0]}-{ramp[1]}" if ramp else "")
            reg = refine_all(coarse, fine_cache[thr], alpha, int(gap), maxp, ramp)
            m = score(ev, reg)
            m["refine_status"] = {k: int(v) for k, v in reg.get(
                "refine_status", pd.Series(dtype=str)).value_counts().items()}
            run["variants"][key] = m
            flat.append({"run": label, "variant": key, **{k: m[k] for k in KEEP},
                         "overshoot": m["median_overshoot_bp"], "undershoot": m["median_undershoot_bp"]})
        res["runs"][label] = run
        print(f"[{label}] {len(run['variants'])} variants done", flush=True)

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    Path(args.out).write_text(json.dumps(res, indent=2) + "\n")
    pd.DataFrame(flat).to_csv(Path(args.out).with_suffix(".csv"), index=False)
    print(f"[done] wrote {args.out} and {Path(args.out).with_suffix('.csv')}", flush=True)


if __name__ == "__main__":
    main()
