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
"""Edges the fine pass leaves on the coarse edge ("stuck"), and what could place them
(docs/att-boundaries.md "Unrefined edges").

A refined model edge is stuck when the fine run at the refinement threshold reaches the
coarse edge, so refinement (which only moves edges inward) leaves it there. Variants,
applied to stuck edges only; every other edge and every att edge is unchanged:

  base          as shipped (refinement 0.6, then att)
  cascade       retry the shipped refinement at stricter fine thresholds (0.8, 0.9, 0.95,
                0.99, in that order) and take the first that moves the edge
  cascade_<t>   only one stricter threshold t
  outward_<b>   follow the 0.6 fine run past the coarse edge (into the scored margin) and
                move the edge a fraction b of the way to the run's end

Ranges: `ranges` reports, for the base calls, the coverage of the shipped (pooled) ranges
and of ranges calibrated separately for stuck and moved edges on --calibrate_from.

Input: Exploratory refinement-on runs (<genome>/windows.tsv with prob_ckpt1..3 and
fine_windows.tsv), as for boundary_eval.py.

Example:
  python scripts/edge_followups.py --runs ../../../arc_data/att_runs/tune42 \\
      --genomes_dir ../../../arc_data/raw/genomes --elements scripts/eval_sizestrat_elements.csv \\
      --out analysis/att/tune42_edges
"""

from __future__ import annotations

import argparse
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import boundary_eval as be
import evalkit as ek
import numpy as np
import pandas as pd

try:
    from prophound.models.prophage.src.att import Genome, apply_boundaries, edge_ranges
except ImportError:
    from src.att import Genome, apply_boundaries, edge_ranges

CASCADE = [0.8, 0.9, 0.95, 0.99]
OUTWARD = [0.5, 1.0]
MODES = ["exploratory", "default", "conservative"]
KEY = ["contig_id", "coarse_start", "coarse_end"]


def build(job: dict, mode: str, fine: pd.DataFrame, w: pd.DataFrame, lengths: dict, fasta: str, thr: float):
    """The shipped calls at one fine threshold (before att)."""
    wl = job["wl"]
    rules = ek.CallRules.for_mode(mode, wl, 3)
    span = (w.end - w.start).to_numpy()
    per = [w[f"prob_ckpt{k}"].to_numpy() for k in (1, 2, 3)]
    ww = w[["genome_id", "contig_id", "start", "end", "prob"]]
    full = span == wl
    coarse = ek.build_pass(ww[full].reset_index(drop=True), [x[full] for x in per], rules.seed_threshold,
                           rules.extend_threshold, rules, max_gap_bp=rules.max_gap_bp, bridge_prob=rules.bridge_prob)
    shorts = []
    for scale in sorted(set(span[span < wl])):
        m = span == scale
        t = rules.short_thresholds[int(scale)]
        shorts.append(ek.build_pass(ww[m].reset_index(drop=True), [x[m] for x in per], t, t, rules, max_gap_bp=0,
                                    kind_suffix=f"_short{int(scale)}"))
    calls = ek.refine_pass(coarse.calls, fine, lengths, job["genome"], fasta, mode,
                           extra=["--refine-threshold", str(thr)])
    allc = ek.concat_nonempty([calls.assign(short=False), *(s.calls.assign(short=True) for s in shorts)]) \
        if len(calls) or shorts else calls
    filt = ek.concat_nonempty([coarse.filtered, *(s.filtered for s in shorts)])
    per_seed = [[pp.per_seed[k] for pp in (coarse, *shorts) if pp.per_seed] for k in range(3)]
    return ek.finalize(allc, filt, per_seed, lengths, lambda df: ek.call_n_fraction(df, fasta), rules,
                       ek.DEFAULT_EDGE_SNAP_FRAC).calls.assign(genome_id=job["genome"])


def stuck_mask(c: pd.DataFrame, side: str) -> np.ndarray:
    edge = c.contig_edge.to_numpy()
    at_end = np.isin(edge, ["start", "both"] if side == "start" else ["end", "both"])
    return ((c.boundary == "model") & c.refine_status.isin(["refined", "split"])
            & (c[side] == c[f"coarse_{side}"])).to_numpy() & ~at_end


def genome_job(job: dict) -> list:
    d, fasta = Path(job["dir"]), job["fasta"]
    w = pd.read_csv(d / "windows.tsv", sep="\t", float_precision="round_trip")
    fine = pd.read_csv(d / "fine_windows.tsv", sep="\t", float_precision="round_trip")
    job["wl"] = int(json.loads((d / "run_metadata.json").read_text())["window_len"])
    lengths = dict(ek.iter_fasta_lengths(fasta))
    genome = Genome.from_fasta(fasta)
    # fine runs at the shipped threshold, unclipped, for the outward variants
    _, runs = ek.build_region_table_coverage(fine.assign(fasta_path=fasta), seed_threshold=0.6, extend_threshold=0.6,
                                             min_coverage_frac=ek.DEFAULT_COVERAGE_FRAC, max_gap_bp=0)
    out = []
    for mode in job["modes"]:
        base = apply_boundaries(build(job, mode, fine, w, lengths, fasta, 0.6), genome, lengths)
        strict = {t: build(job, mode, fine, w, lengths, fasta, t) for t in CASCADE}
        strict = {t: c[c.n_pieces == 1].drop_duplicates(KEY).set_index(KEY) for t, c in strict.items()}
        variants = {"base": base}
        for name, order in [("cascade", CASCADE), *((f"cascade_{t}", [t]) for t in CASCADE)]:
            v = base.copy()
            for side in ("start", "end"):
                m = stuck_mask(v, side) & (v.n_pieces == 1).to_numpy()
                for i in np.flatnonzero(m):
                    k = tuple(v.iloc[i][KEY])
                    for t in order:
                        if k in strict[t].index and strict[t].loc[k, side] != v.iloc[i][side]:
                            v.iloc[i, v.columns.get_loc(side)] = int(strict[t].loc[k, side])
                            break
            variants[name] = v
        for b in OUTWARD:
            v = base.copy()
            for side in ("start", "end"):
                for i in np.flatnonzero(stuck_mask(v, side)):
                    r = v.iloc[i]
                    e = int(r[side])
                    on = runs[(runs.contig_id == r.contig_id) & (runs.start < e + 1) & (runs.end > e - 1)]
                    if on.empty:
                        continue
                    far = int(on.start.min()) if side == "start" else int(on.end.max())
                    if (side == "start" and far < e) or (side == "end" and far > e):
                        v.iloc[i, v.columns.get_loc(side)] = int(round(e + b * (far - e)))
            variants[f"outward_{b:g}"] = v
        for name, v in variants.items():
            v = v.assign(region_len_bp=v.end - v.start)
            v = edge_ranges(v.drop(columns=["start_lo", "start_hi", "end_lo", "end_hi"]), lengths)
            v["stuck_start"], v["stuck_end"] = stuck_mask(base, "start"), stuck_mask(base, "end")
            out.append((mode, name, v))
    return out


def edge_table(calls: pd.DataFrame, ev: pd.DataFrame) -> pd.DataFrame:
    """One row per real edge of a call on exactly one element: inward-positive error and range hit."""
    rows = []
    for c in calls.itertuples():
        t = ev[(ev.accession == c.genome_id) & (ev.contig == c.contig_id) & (ev.end > c.start) & (ev.start < c.end)]
        if len(t) != 1 or c.boundary != "model":
            continue
        y = t.iloc[0]
        for side, true, stuck in (("start", int(y.start), c.stuck_start), ("end", int(y.end), c.stuck_end)):
            if c.contig_edge in ((side, "both")):
                continue
            val = int(getattr(c, side))
            dev = (true - val) if side == "start" else (val - true)
            rows.append({"side": side, "stuck": bool(stuck), "L": int(c.end - c.start), "dev": dev,
                         "hit": int(getattr(c, f"{side}_lo")) <= true <= int(getattr(c, f"{side}_hi"))})
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", required=True)
    ap.add_argument("--genomes_dir", required=True)
    ap.add_argument("--elements", required=True)
    ap.add_argument("--modes", nargs="+", default=MODES)
    ap.add_argument("--calibrate_from", help="a previous --out prefix: its <prefix>.edges.tsv calibrates the ranges")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    dirs = sorted(d for d in Path(a.runs).iterdir() if (d / "fine_windows.tsv").exists())
    jobs = [{"genome": d.name, "dir": str(d), "fasta": str(Path(a.genomes_dir) / f"{d.name}.fna"), "modes": a.modes}
            for d in dirs]
    by = {}
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        for res in ex.map(genome_job, jobs):
            for mode, name, c in res:
                by.setdefault((mode, name), []).append(c)
    truth = pd.read_csv(a.elements)
    ev = truth[truth.accession.isin([d.name for d in dirs])]
    rows, edges = [], []
    for (mode, name), parts in by.items():
        reg = ek.concat_nonempty(parts)
        e = edge_table(reg, ev)
        st = e[e.stuck]
        rows.append({"mode": mode, "variant": name, **be.metrics(ev, reg, len(dirs)),
                     "stuck_edges": int(len(st)), "stuck_median_dev_bp": float(st.dev.median()) if len(st) else np.nan,
                     "stuck_undershoot": float((st.dev < 0).mean()) if len(st) else np.nan,
                     "stuck_range_cover": float(st.hit.mean()) if len(st) else np.nan,
                     "all_range_cover": float(e.hit.mean()) if len(e) else np.nan})
        if name == "base":
            edges.append(e.assign(mode=mode))
    edges = pd.concat(edges, ignore_index=True)
    edges.to_csv(f"{a.out}.edges.tsv", sep="\t", index=False)
    res = pd.DataFrame(rows)
    res.to_csv(f"{a.out}.csv", index=False)
    cols = ["mode", "variant", "covered50", "covered", "matched", "median_iou", "Mb_called", "Mb_off", "stuck_edges",
            "stuck_median_dev_bp", "stuck_undershoot", "stuck_range_cover", "all_range_cover"]
    print(res[cols].round(3).to_string(index=False))
    # option 2: ranges calibrated per class (stuck / moved) on --calibrate_from, or on this sample
    cal = pd.read_csv(f"{a.calibrate_from}.edges.tsv", sep="\t") if a.calibrate_from else edges
    q = {s: ((g.dev / g.L).quantile(0.05), (g.dev / g.L).quantile(0.95)) for s, g in cal.groupby("stuck")}
    rel = edges.dev / edges.L
    lo = edges.stuck.map(lambda s: q[s][0])
    hi = edges.stuck.map(lambda s: q[s][1])
    two = (rel >= lo) & (rel <= hi)
    fracs = {str(k): tuple(round(x, 3) for x in v) for k, v in q.items()}
    print("\nranges calibrated per class on", a.calibrate_from or "this sample", fracs)
    for s, g in edges.assign(two=two.to_numpy(), w=((hi - lo) * edges.L).to_numpy()).groupby("stuck"):
        print(f"  {'stuck' if s else 'moved'} edges n={len(g)}: pooled ranges cover {g.hit.mean():.3f}; "
              f"per-class ranges cover {g.two.mean():.3f}, median width {g.w.median() / 1e3:.1f} kb")


if __name__ == "__main__":
    main()
