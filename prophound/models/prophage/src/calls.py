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
"""From window scores to calls: the call pipeline shared by `predict` and the offline scorers.

No model and no FASTA here, only score tracks, so `scripts/mode_recheck.py` rebuilds
calls with exactly the code the CLI runs. Each scoring pass (the coarse grid, then each
short-contig scale) goes through `build_pass`; `finalize` then snaps edges, measures N
content and seed support, and applies the N filter. Calls that fail the seed vote or
the N filter are kept aside with a `filter` reason instead of being deleted.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Mapping, Tuple

import numpy as np
import pandas as pd

from .options import (
    DEFAULT_BRIDGE_PROB,
    DEFAULT_COVERAGE_FRAC,
    DEFAULT_EDGE_SNAP_FRAC,
    DEFAULT_MAX_GAP_BP,
    MAX_CALL_N_FRAC,
    MODE_MIN_SEED_SUPPORT,
    MODE_THRESHOLDS,
    SHORT_CONTIG_THRESHOLDS,
)

REGION_COLUMNS = ["genome_id", "contig_id", "start", "end", "n_windows",
                  "region_len_bp", "mean_prob", "max_prob", "boundary_kind"]
WINDOW_COORDS = ["genome_id", "contig_id", "start", "end"]
FILTER_SEED_VOTE = "seed_vote"
FILTER_N_GAP = "n_gap"


def _coverage_boundary_kind(frac: float, max_gap_bp: int, bridge_prob: float) -> str:
    """`boundary_kind` for the coverage builder: its settings, so a regions.tsv says how it was built."""
    kind = f"coverage_frac_{frac:g}"
    if max_gap_bp > 0:
        kind += f"_gap{max_gap_bp:d}"
        if bridge_prob > 0:
            kind += f"_br{bridge_prob:g}"
    return kind


def _close_gaps(
    inside: np.ndarray,
    bin_mean: np.ndarray,
    covered: np.ndarray,
    max_gap_bins: int,
    bridge_prob: float,
) -> np.ndarray:
    """Merge inside-runs separated by a short gap that still looks element-like.

    A morphological closing of the per-bin `inside` mask: a run of outside bins is filled
    if it is at most `max_gap_bins` long and every bin in it is covered with a mean window
    probability of at least `bridge_prob`. This repairs large elements whose score dips
    partway through. `bridge_prob` is the guard against bridging two
    separate neighbouring elements; at 0 the merge is distance-only. Gaps at a contig end
    are never filled, and an uncovered bin never passes a `bridge_prob > 0` gate.
    """
    if max_gap_bins <= 0:
        return inside

    out = inside.copy()
    n = out.shape[0]
    idx = 0
    # Skip a leading outside run: nothing to its left to merge with.
    while idx < n and not out[idx]:
        idx += 1
    while idx < n:
        while idx < n and out[idx]:
            idx += 1
        if idx >= n:
            break
        gap_start = idx
        while idx < n and not out[idx]:
            idx += 1
        if idx >= n:
            break  # trailing outside run, no right-hand neighbour
        gap_end = idx  # exclusive
        if gap_end - gap_start <= max_gap_bins:
            if bridge_prob <= 0:
                out[gap_start:gap_end] = True
            else:
                span = slice(gap_start, gap_end)
                if covered[span].all() and (bin_mean[span] >= bridge_prob).all():
                    out[span] = True
    return out


def _bin_grid(contig_df: pd.DataFrame) -> Tuple[int, int, int, int]:
    """Stride-sized bins over one contig's windows (sorted by start): stride (the
    commonest start-to-start step), contig end, bin count, bins covered per window."""
    stride_values = contig_df["start"].astype(int).diff().dropna()
    window_len = int(contig_df["end"].iloc[0] - contig_df["start"].iloc[0])
    stride = int(stride_values.mode().iloc[0]) if not stride_values.empty else window_len
    contig_end = int(contig_df["end"].max())
    n_bins = int(math.ceil(contig_end / stride))
    return stride, contig_end, n_bins, max(1, int(math.ceil(window_len / stride)))


def _runs(mask: np.ndarray) -> List[Tuple[int, int]]:
    """(first, last) bin index of each run of True in `mask`."""
    out, i, n = [], 0, len(mask)
    while i < n:
        if not mask[i]:
            i += 1
            continue
        j = i
        while j + 1 < n and mask[j + 1]:
            j += 1
        out.append((i, j))
        i = j + 1
    return out


def build_region_table_coverage(
    windows_df: pd.DataFrame,
    seed_threshold: float,
    extend_threshold: float,
    min_coverage_frac: float = DEFAULT_COVERAGE_FRAC,
    max_gap_bp: int = DEFAULT_MAX_GAP_BP,
    bridge_prob: float = DEFAULT_BRIDGE_PROB,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Turn per-window scores into regions, per contig, on a grid of stride-sized bins.

    A bin is inside if at least `min_coverage_frac` of the windows covering it score at or
    above `extend_threshold` (0 approaches the union of high windows, 1 the intersection).
    The model's score is spiky inside real elements, so requiring the mean of the covering
    windows instead drops sparse-signal elements and trims every call by about half a
    window per side. With `max_gap_bp` > 0, runs separated
    by a short gap are merged first (`_close_gaps`). A run of inside bins becomes a region if
    any bin in it is covered by a window at or above `seed_threshold`.

    Returns the windows, sorted by start within each contig, with `region_pred_label` (1 if
    the window overlaps a region), and the regions (`REGION_COLUMNS`). Region edges are bin
    edges.
    """
    if not 0.0 <= min_coverage_frac <= 1.0:
        raise ValueError("min_coverage_frac must be in [0, 1]")
    if max_gap_bp < 0:
        raise ValueError("max_gap_bp must be >= 0")
    if not 0.0 <= bridge_prob <= 1.0:
        raise ValueError("bridge_prob must be in [0, 1]")

    annotated_parts: List[pd.DataFrame] = []
    region_rows: List[Dict[str, object]] = []

    for (genome_id, contig_id), contig_df in windows_df.groupby(
        ["genome_id", "contig_id"],
        sort=False,
    ):
        contig_df = contig_df.sort_values(["start", "end"]).reset_index(drop=True).copy()
        contig_df["region_pred_label"] = 0

        stride, contig_end, n_bins, bins_per_window = _bin_grid(contig_df)

        high_count = np.zeros(n_bins, dtype=np.int32)
        cover_count = np.zeros(n_bins, dtype=np.int32)
        prob_sum = np.zeros(n_bins, dtype=np.float64)
        seed_bin_mask = np.zeros(n_bins, dtype=bool)
        starts = contig_df["start"].astype(int).to_numpy()
        probs = contig_df["prob"].astype(float).to_numpy()
        for start, prob in zip(starts, probs):
            start_bin = int(start // stride)
            end_bin = min(n_bins, start_bin + bins_per_window)
            cover_count[start_bin:end_bin] += 1
            prob_sum[start_bin:end_bin] += prob
            if prob >= extend_threshold:
                high_count[start_bin:end_bin] += 1
            if prob >= seed_threshold:
                seed_bin_mask[start_bin:end_bin] = True

        covered = cover_count > 0
        frac = np.zeros(n_bins, dtype=np.float64)
        frac[covered] = high_count[covered] / cover_count[covered]
        inside = covered & (high_count > 0) & (frac >= min_coverage_frac)

        if max_gap_bp > 0:
            bin_mean = np.zeros(n_bins, dtype=np.float64)
            bin_mean[covered] = prob_sum[covered] / cover_count[covered]
            inside = _close_gaps(
                inside,
                bin_mean=bin_mean,
                covered=covered,
                max_gap_bins=int(max_gap_bp // stride),
                bridge_prob=bridge_prob,
            )

        for bin_index, right_bin in _runs(inside):
            if seed_bin_mask[bin_index : right_bin + 1].any():
                region_start = bin_index * stride
                region_end = min(contig_end, (right_bin + 1) * stride)
                overlap_mask = (contig_df["end"].astype(int) > region_start) & (
                    contig_df["start"].astype(int) < region_end
                )
                contig_df.loc[overlap_mask, "region_pred_label"] = 1
                overlapping = contig_df[overlap_mask]
                region_rows.append(
                    {
                        "genome_id": genome_id,
                        "contig_id": contig_id,
                        "start": int(region_start),
                        "end": int(region_end),
                        "n_windows": int(len(overlapping)),
                        "region_len_bp": int(region_end - region_start),
                        "mean_prob": float(overlapping["prob"].astype(float).mean()),
                        "max_prob": float(overlapping["prob"].astype(float).max()),
                        "boundary_kind": _coverage_boundary_kind(
                            min_coverage_frac, max_gap_bp, bridge_prob
                        ),
                    }
                )

        annotated_parts.append(contig_df)

    annotated_windows_df = pd.concat(annotated_parts, ignore_index=True)
    regions_df = pd.DataFrame(region_rows, columns=REGION_COLUMNS)
    return annotated_windows_df, regions_df


def seed_support(regions_df: pd.DataFrame, per_seed_regions: List[pd.DataFrame]) -> List[int]:
    """Per call, how many checkpoints make an overlapping call on their own.

    With the `max` ensemble any one seed can create a
    call, and on the 42-genome sample calls only one of three seeds makes were
    99% off-target (152 of 154). Reported so users can rank calls; no call is
    dropped here.
    """
    out = []
    for r in regions_df.itertuples(index=False):
        n = 0
        for sr in per_seed_regions:
            if len(sr) and ((sr["contig_id"] == r.contig_id) & (sr["end"].astype(int) > int(r.start))
                            & (sr["start"].astype(int) < int(r.end))).any():
                n += 1
        out.append(n)
    return out


def vote_extent(regions_df: pd.DataFrame) -> pd.DataFrame:
    """The extent the seed vote was taken on: `coarse_start`/`coarse_end` for refined
    calls and split pieces, the call itself otherwise. Support counted on a piece alone
    can read below the vote its parent call passed."""
    if "coarse_start" not in regions_df.columns or regions_df.empty:
        return regions_df
    has = regions_df["coarse_start"].notna().to_numpy()
    return regions_df.assign(start=np.where(has, regions_df["coarse_start"], regions_df["start"]).astype(np.int64),
                             end=np.where(has, regions_df["coarse_end"], regions_df["end"]).astype(np.int64))


def apply_seed_vote(regions_df: pd.DataFrame, per_seed_regions: List[pd.DataFrame],
                    min_support: int) -> Tuple[pd.DataFrame, int]:
    """Keep calls that at least `min_support` checkpoints make on their own.

    The call-level seed vote (docs/cli-and-packaging.md "Modes"): extent comes from the
    ensemble (`max`) call, the vote only decides whether it is kept. Returns the
    kept calls and the number dropped. min_support <= 1 keeps everything.
    """
    if min_support <= 1 or regions_df.empty:
        return regions_df, 0
    sup = np.asarray(seed_support(regions_df, per_seed_regions))
    keep = sup >= min_support
    return regions_df[keep].reset_index(drop=True), int((~keep).sum())


def snap_contig_edges(regions_df: pd.DataFrame, contig_lengths: Dict[str, int],
                      snap_bp: int) -> Tuple[pd.DataFrame, int]:
    """Move call edges within `snap_bp` of a contig end onto that end, and flag
    every edge that sits on a contig end afterwards.

    `contig_edge` is one of none/start/end/both whether or not the edge moved:
    either way the reported boundary is the end of the sequence, not a detected
    element edge, and the element may continue on another contig. Returns the
    calls and how many moved.
    """
    if regions_df.empty:
        return regions_df.assign(contig_edge=pd.Series(dtype=str)), 0
    out = regions_df.copy()
    L = out["contig_id"].astype(str).map(contig_lengths).to_numpy(dtype=np.int64)
    start, end = out["start"].to_numpy(dtype=np.int64), out["end"].to_numpy(dtype=np.int64)
    new_start = np.where(start <= snap_bp, 0, start)
    new_end = np.where(end >= L - snap_bp, L, end)
    moved = int(((new_start != start) | (new_end != end)).sum())
    out["start"], out["end"] = new_start, new_end
    if "region_len_bp" in out.columns:
        out["region_len_bp"] = new_end - new_start
    at_s, at_e = new_start == 0, new_end == L
    out["contig_edge"] = np.select([at_s & at_e, at_s, at_e], ["both", "start", "end"], "none")
    return out, moved


def label_windows(windows_df: pd.DataFrame, regions_df: pd.DataFrame) -> pd.DataFrame:
    """Set region_pred_label to 1 on windows overlapping any call, else 0."""
    windows_df = windows_df.copy()
    windows_df["region_pred_label"] = 0
    for r in regions_df.itertuples(index=False):
        m = ((windows_df["contig_id"] == r.contig_id)
             & (windows_df["end"].astype(int) > int(r.start)) & (windows_df["start"].astype(int) < int(r.end)))
        windows_df.loc[m, "region_pred_label"] = 1
    return windows_df


def concat_nonempty(frames: List[pd.DataFrame]) -> pd.DataFrame:
    """Concatenate, skipping empty frames (pandas would let them set column dtypes); the
    first frame when all are empty."""
    nonempty = [f for f in frames if len(f)]
    if not nonempty:
        return frames[0]
    return pd.concat(nonempty, ignore_index=True) if len(nonempty) > 1 else nonempty[0].reset_index(drop=True)


@dataclass(frozen=True)
class CallRules:
    """Every setting that turns score tracks into calls.

    `min_seed_support` is the vote actually applied: already capped at the number of
    checkpoints, and 1 with a single checkpoint.
    """

    seed_threshold: float
    extend_threshold: float
    coverage_frac: float = DEFAULT_COVERAGE_FRAC
    max_gap_bp: int = DEFAULT_MAX_GAP_BP
    bridge_prob: float = DEFAULT_BRIDGE_PROB
    min_seed_support: int = 1
    short_thresholds: Mapping[int, float] = field(default_factory=dict)
    snap_bp: int = 0
    max_call_n_frac: float = MAX_CALL_N_FRAC

    @classmethod
    def for_mode(cls, mode: str, window_len: int, n_checkpoints: int) -> CallRules:
        """The shipped defaults for `mode` (what `predict --mode` runs with no overrides)."""
        seed = MODE_THRESHOLDS[mode]
        return cls(seed_threshold=seed, extend_threshold=min(seed, 0.99),
                   min_seed_support=min(MODE_MIN_SEED_SUPPORT[mode], n_checkpoints) if n_checkpoints > 1 else 1,
                   short_thresholds=dict(SHORT_CONTIG_THRESHOLDS[mode]),
                   snap_bp=int(DEFAULT_EDGE_SNAP_FRAC * window_len))


@dataclass
class PassCalls:
    """One scoring pass turned into calls."""

    windows: pd.DataFrame          # the pass's windows, sorted per contig, with region_pred_label
    calls: pd.DataFrame            # calls that passed the vote
    filtered: pd.DataFrame         # calls that failed it, with a `filter` column
    per_seed: List[pd.DataFrame]   # each checkpoint's own calls (empty with one checkpoint)


def build_pass(
    windows: pd.DataFrame,
    per_checkpoint: List[np.ndarray],
    seed_threshold: float,
    extend_threshold: float,
    rules: CallRules,
    max_gap_bp: int,
    bridge_prob: float = 0.0,
    kind_suffix: str = "",
) -> PassCalls:
    """Calls from one pass: the coverage builder on the ensemble track (`windows.prob`),
    the same builder on each checkpoint's track (`per_checkpoint`, aligned to `windows`'
    rows), and the call-level seed vote. The vote keeps a call's ensemble extent and
    only decides whether it passes. `kind_suffix` is appended to `boundary_kind`
    (the short-contig passes add `_short<scale>`).
    """
    annotated, calls = build_region_table_coverage(
        windows, seed_threshold=seed_threshold, extend_threshold=extend_threshold,
        min_coverage_frac=rules.coverage_frac, max_gap_bp=max_gap_bp, bridge_prob=bridge_prob,
    )
    per_seed = []
    if len(per_checkpoint) > 1:
        coords = windows[WINDOW_COORDS]
        per_seed = [build_region_table_coverage(
            coords.assign(prob=p), seed_threshold=seed_threshold, extend_threshold=extend_threshold,
            min_coverage_frac=rules.coverage_frac, max_gap_bp=max_gap_bp, bridge_prob=bridge_prob,
        )[1] for p in per_checkpoint]
    if kind_suffix:
        calls["boundary_kind"] = calls["boundary_kind"].astype(str) + kind_suffix
    filtered = calls.iloc[0:0].assign(filter=pd.Series(dtype=str))
    if rules.min_seed_support > 1 and len(calls):
        keep = np.asarray(seed_support(calls, per_seed)) >= rules.min_seed_support
        filtered = calls[~keep].reset_index(drop=True).assign(filter=FILTER_SEED_VOTE)
        calls = calls[keep].reset_index(drop=True)
    return PassCalls(annotated, calls, filtered, per_seed)


@dataclass
class FinalCalls:
    calls: pd.DataFrame       # what regions.tsv reports
    filtered: pd.DataFrame    # every candidate a filter removed, with the reason(s) in `filter`
    meta: Dict[str, dict]     # edge_snap, n_gap_filter and (with several checkpoints) seed_support


def finalize(
    calls: pd.DataFrame,
    filtered: pd.DataFrame,
    per_seed: List[List[pd.DataFrame]],
    contig_lengths: Dict[str, int],
    n_fraction: Callable[[pd.DataFrame], np.ndarray],
    rules: CallRules,
    edge_snap_frac: float,
    requested_min_seed_support: int | None = None,
) -> FinalCalls:
    """Snap edges at contig ends, measure N content, apply the N filter, count seed support.

    `per_seed[k]` holds checkpoint k's calls from every pass; with one checkpoint it is
    empty and no `seed_support` column is added. Filtered candidates get the same columns
    (snapped edges, `n_frac`, `seed_support`) so they can be judged like any call.
    Snapping comes first, so `n_frac` describes the reported coordinates. Support is
    counted on the extent the vote saw (`vote_extent`), so a split piece reports its
    parent call's support.
    """
    calls, n_snapped = snap_contig_edges(calls, contig_lengths, rules.snap_bp)
    edge_meta = {"frac": edge_snap_frac, "snap_bp": rules.snap_bp, "calls_moved": n_snapped,
                 "calls_on_contig_edge": int((calls["contig_edge"] != "none").sum()) if len(calls) else 0}

    n_frac = n_fraction(calls)
    gappy = n_frac >= rules.max_call_n_frac
    calls = calls.assign(n_frac=np.round(n_frac, 4))
    n_dropped = int(gappy.sum())
    if n_dropped:
        filtered = concat_nonempty([filtered, calls[gappy].assign(filter=FILTER_N_GAP)])
        calls = calls[~gappy].reset_index(drop=True)
    gap_meta = {"max_call_n_frac": rules.max_call_n_frac, "calls_dropped": n_dropped,
                "calls_with_any_n": int((calls.n_frac > 0).sum())}

    voted_out = filtered["filter"] == FILTER_SEED_VOTE
    if voted_out.any():
        snapped, _ = snap_contig_edges(filtered[voted_out].drop(columns="filter"), contig_lengths, rules.snap_bp)
        snapped = snapped.assign(n_frac=np.round(n_fraction(snapped), 4))
        snapped["filter"] = np.where(snapped.n_frac >= rules.max_call_n_frac,
                                     f"{FILTER_SEED_VOTE};{FILTER_N_GAP}", FILTER_SEED_VOTE)
        filtered = concat_nonempty([snapped, filtered[~voted_out]])

    meta: Dict[str, dict] = {"edge_snap": edge_meta, "n_gap_filter": gap_meta}
    if per_seed:
        seed_regs = [pd.concat([x for x in r if len(x)], ignore_index=True) if any(len(x) for x in r)
                     else pd.DataFrame(columns=["contig_id", "start", "end"]) for r in per_seed]
        calls = calls.copy()
        calls["seed_support"] = seed_support(vote_extent(calls), seed_regs) if len(calls) else []
        if len(filtered):
            filtered = filtered.assign(seed_support=seed_support(vote_extent(filtered), seed_regs))
        support: Dict[str, object] = {
            "n_checkpoints": len(per_seed),
            "min_seed_support": rules.min_seed_support,
            "requested_min_seed_support": (rules.min_seed_support if requested_min_seed_support is None
                                           else requested_min_seed_support),
            "calls_dropped_by_vote": int((filtered["filter"].str.startswith(FILTER_SEED_VOTE)).sum())
            if len(filtered) else 0,
            "calls_by_support": {str(k): int(v) for k, v in
                                 calls["seed_support"].value_counts().sort_index().items()},
        }
        meta["seed_support"] = support
    return FinalCalls(calls, filtered, meta)
