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
"""Fine-scale refinement of coarse region calls (pass B).

docs/cli-and-packaging.md "What a run does".
The coarse pass (10 kb windows, coverage builder, gap-fill merge) handles
detection and fragmentation. This pass re-scores each coarse call +/- a margin
at a finer window size and addresses the other three downstream problems:

  - conflation: a coarse call is SPLIT where the fine track has a gap of at
    least `split_min_gap_bp` whose windows look like background (mean prob
    below `split_max_prob`). The background test is what keeps splitting from
    undoing gap-filling, which bridged coarse gaps whose mean was >= bridge_prob.
  - truncation / overextension: each piece's edge is placed between its fine
    core (known to truncate: the fine scan fires on the element's core) and its
    outer bound (the coarse edge, known to overextend by ~half a coarse window,
    or the midpoint of a split gap):  edge = fine + alpha * (outer - fine).
    alpha=1 reproduces the coarse call exactly; alpha=0 is the fine core. So a
    refined piece is never wider than its coarse call nor narrower than its
    fine core -- neither failure mode can exceed its current bound.
    alpha may ramp with the piece's fine-core length (`alpha_ramp_bp`): the fine
    core suits small elements and the coarse edge large ones, so one global
    value cannot serve both.
  - recall: a coarse call with no fine support is kept unchanged and flagged,
    never dropped.

Pure logic only -- no model, no FASTA -- so it is unit-testable on synthetic
score tracks. predict_genome.py plans the fine windows, scores them, runs the
coverage builder over them and hands the result here.
"""

from __future__ import annotations

import math
from typing import Dict, List, Optional, Tuple

import pandas as pd

REFINE_COLUMNS = [
    "genome_id", "contig_id", "start", "end", "n_windows", "region_len_bp",
    "mean_prob", "max_prob", "boundary_kind",
    "refine_status", "coarse_start", "coarse_end", "n_pieces",
]


def plan_neighborhoods(
    regions_df: pd.DataFrame,
    contig_lengths: Dict[str, int],
    margin_bp: int,
) -> Dict[str, List[Tuple[int, int]]]:
    """Per contig, the union of [start - margin, end + margin) over its calls."""
    out: Dict[str, List[Tuple[int, int]]] = {}
    for contig_id, grp in regions_df.groupby("contig_id", sort=False):
        length = int(contig_lengths[contig_id])
        ivs = sorted(
            (max(0, int(s) - margin_bp), min(length, int(e) + margin_bp))
            for s, e in zip(grp["start"], grp["end"])
        )
        merged: List[Tuple[int, int]] = []
        for s, e in ivs:
            if merged and s <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], e))
            else:
                merged.append((s, e))
        out[contig_id] = merged
    return out


def plan_fine_windows(
    neighborhoods: Dict[str, List[Tuple[int, int]]],
    contig_lengths: Dict[str, int],
    genome_id: str,
    fasta_path: str,
    window_len: int,
    stride: int,
) -> pd.DataFrame:
    """Fine windows tiling each neighborhood on one global grid per contig.

    Starts are multiples of `stride`, so windows from different neighborhoods
    never disagree about bin alignment and the coverage builder's absolute
    `start // stride` binning holds. Every window lies fully inside its contig.
    A neighborhood shorter than one window gets the single grid window nearest
    its centre; a contig shorter than one window gets none.
    """
    rows = []
    for contig_id, ivs in neighborhoods.items():
        length = int(contig_lengths[contig_id])
        max_start = ((length - window_len) // stride) * stride
        if length < window_len:
            continue
        seen = set()
        for s, e in ivs:
            first = (s // stride) * stride
            last = min(e - window_len, max_start)
            if last < first:
                centre = (s + e) // 2 - window_len // 2
                first = last = min(max_start, max(0, (centre // stride) * stride))
            for start in range(first, last + 1, stride):
                if start not in seen:
                    seen.add(start)
                    rows.append({"genome_id": genome_id, "contig_id": contig_id, "start": int(start),
                                 "end": int(start + window_len), "fasta_path": fasta_path})
    cols = ["genome_id", "contig_id", "start", "end", "fasta_path"]
    return pd.DataFrame(rows, columns=cols)


def _group_runs(runs: List[Tuple[int, int]], fine_df: pd.DataFrame,
                split_min_gap_bp: int, split_max_prob: float) -> List[Tuple[int, int]]:
    """Merge fine runs unless the gap between them is long AND background."""
    if not runs:
        return []
    pieces = [list(runs[0])]
    centers = (fine_df["start"] + fine_df["end"]) // 2
    for s, e in runs[1:]:
        gap_s, gap_e = pieces[-1][1], s
        in_gap = fine_df[(centers >= gap_s) & (centers < gap_e)]
        background = (not in_gap.empty) and float(in_gap["prob"].mean()) < split_max_prob
        if gap_e - gap_s >= split_min_gap_bp and background:
            pieces.append([s, e])
        else:
            pieces[-1][1] = max(pieces[-1][1], e)
    return [(a, b) for a, b in pieces]


def ramp_alpha(core_len: int, alpha_max: float, ramp_bp: Optional[Tuple[int, int]]) -> float:
    """alpha for a piece whose fine core spans `core_len` bp.

    Without a ramp, `alpha_max` everywhere (the constant --edge-alpha). With
    ramp (lo, hi): 0 at or below lo, `alpha_max` at or above hi, linear in
    log(length) between -- element sizes span two orders of magnitude, so the
    ramp is set in log space.
    """
    if ramp_bp is None:
        return alpha_max
    lo, hi = ramp_bp
    if core_len <= lo:
        return 0.0
    if core_len >= hi:
        return alpha_max
    return alpha_max * math.log(core_len / lo) / math.log(hi / lo)


def refine_regions(
    coarse_df: pd.DataFrame,
    fine_windows_df: pd.DataFrame,
    fine_regions_df: pd.DataFrame,
    edge_alpha: float,
    split_min_gap_bp: int,
    split_max_prob: float,
    kind_suffix: str = "",
    alpha_ramp_bp: Optional[Tuple[int, int]] = None,
) -> pd.DataFrame:
    """Split and re-edge each coarse call using fine runs. See module docstring."""
    if not 0.0 <= edge_alpha <= 1.0:
        raise ValueError("edge_alpha must be in [0, 1]")
    if alpha_ramp_bp is not None and not 0 < alpha_ramp_bp[0] < alpha_ramp_bp[1]:
        raise ValueError("alpha_ramp_bp must be (lo, hi) with 0 < lo < hi")
    rows = []
    for r in coarse_df.itertuples(index=False):
        cs, ce = int(r.start), int(r.end)
        c_fine = fine_windows_df[(fine_windows_df["contig_id"] == r.contig_id)
                                 & (fine_windows_df["end"] > cs) & (fine_windows_df["start"] < ce)]
        runs = sorted(
            (max(cs, int(f.start)), min(ce, int(f.end)))
            for f in fine_regions_df[(fine_regions_df["contig_id"] == r.contig_id)
                                     & (fine_regions_df["end"] > cs)
                                     & (fine_regions_df["start"] < ce)].itertuples(index=False)
        )
        runs = [(s, e) for s, e in runs if e > s]
        base = {"genome_id": r.genome_id, "contig_id": r.contig_id,
                "coarse_start": cs, "coarse_end": ce}
        if not runs:
            rows.append({**base, "start": cs, "end": ce, "n_windows": int(r.n_windows),
                         "region_len_bp": ce - cs, "mean_prob": float(r.mean_prob),
                         "max_prob": float(r.max_prob), "boundary_kind": f"{r.boundary_kind}{kind_suffix}",
                         "refine_status": "no_fine_support", "n_pieces": 1})
            continue
        pieces = _group_runs(runs, c_fine, split_min_gap_bp, split_max_prob)
        n = len(pieces)
        for i, (fs, fe) in enumerate(pieces):
            lo = cs if i == 0 else (pieces[i - 1][1] + fs) // 2
            hi = ce if i == n - 1 else (fe + pieces[i + 1][0]) // 2
            alpha = ramp_alpha(fe - fs, edge_alpha, alpha_ramp_bp)
            start = int(round(fs - alpha * (fs - lo)))
            end = int(round(fe + alpha * (hi - fe)))
            inside = c_fine[(c_fine["end"] > start) & (c_fine["start"] < end)]
            rows.append({**base, "start": start, "end": end, "n_windows": int(len(inside)),
                         "region_len_bp": end - start,
                         "mean_prob": float(inside["prob"].mean()) if len(inside) else float(r.mean_prob),
                         "max_prob": float(inside["prob"].max()) if len(inside) else float(r.max_prob),
                         "boundary_kind": f"{r.boundary_kind}{kind_suffix}",
                         "refine_status": "split" if n > 1 else "refined", "n_pieces": n})
    return pd.DataFrame(rows, columns=REFINE_COLUMNS)
