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
    """`boundary_kind` string for the coverage builder, encoding the merge
    settings when fragment merging is on (so a regions.tsv is self-describing).
    Unchanged from the pre-merge format when it is off."""
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
    """Merge inside-runs separated by a short, optionally still-elevated gap.

    The fragmentation failure mode (docs/gap-filling.md, measured 2026-09-20 in
    docs/model-quality-baseline.md): 60 kb+ elements are matched 86-100% of the
    time but at median IoU 0.218 and median predicted length 22 kb, because the
    score track dips below `extend_threshold` partway through and the builder
    closes the region there. This is a morphological closing on the per-bin
    `inside` mask: a run of outside bins is filled if it is at most
    `max_gap_bins` long AND every covered bin in it has mean covering-window
    probability >= `bridge_prob`.

    `bridge_prob` is the anti-overextension guard (docs/anti-overextension.md):
    at 0 the merge is distance-only and will happily bridge two genuinely
    separate nearby elements; raising it requires the intervening sequence to
    still look somewhat element-like, well below `extend_threshold`. Gaps at a
    contig edge (leading/trailing) are never filled -- there is no second run to
    merge with. Uncovered bins inside a gap cannot be judged, so they never
    satisfy a `bridge_prob > 0` gate.
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


def build_region_table_coverage(
    windows_df: pd.DataFrame,
    seed_threshold: float,
    extend_threshold: float,
    min_coverage_frac: float = DEFAULT_COVERAGE_FRAC,
    max_gap_bp: int = DEFAULT_MAX_GAP_BP,
    bridge_prob: float = DEFAULT_BRIDGE_PROB,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Coverage-fraction region builder (2026-09-17, docs/multi-model-architecture.md).

    The removed `strict` builder marked a stride-sized bin as inside a region only if
    the MEAN of the ~window_len/stride windows covering it is >= `extend_threshold`
    -- effectively the intersection of high-scoring windows. Because the model's
    score track inside real elements is spiky (islands of high windows amid
    background), that drops sparse-signal elements outright and trims every
    region by ~half a window per side. This builder instead marks a bin as
    inside if at least `min_coverage_frac` of its covering windows individually
    score >= `extend_threshold` (frac -> 0 is the union of high windows; frac ~1
    recovers intersection behavior). Runs of inside bins that contain at least
    one bin covered by a window >= `seed_threshold` become regions.

    On the 18-genome eval set this lifts region recall 76.9% -> 88.5% (frac 0.3)
    / 94.2% (frac <= 0.1) and median IoU 0.56 -> 0.75-0.77, at the cost of
    2-4x more calls per genome. frac 0.3 is where start/end errors are
    symmetric; the union brackets true edges by ~4.5 kb a side. Same output
    contract as the strict builder had; `region_pred_label` marks windows that
    overlap a called region.

    `max_gap_bp` / `bridge_prob` (both default off) additionally merge calls
    separated by a short gap, to repair elements fragmented by an internal dip
    in the score track -- see `_close_gaps` and docs/gap-filling.md. With
    `max_gap_bp=0` this function's output is identical to its pre-2026-09-21
    behavior.
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

        stride_values = contig_df["start"].astype(int).diff().dropna()
        stride = int(stride_values.mode().iloc[0]) if not stride_values.empty else int(
            contig_df["end"].iloc[0] - contig_df["start"].iloc[0]
        )
        window_len = int(contig_df["end"].iloc[0] - contig_df["start"].iloc[0])
        contig_end = int(contig_df["end"].max())
        n_bins = int(math.ceil(contig_end / stride))
        bins_per_window = max(1, int(math.ceil(window_len / stride)))

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

        bin_index = 0
        while bin_index < n_bins:
            if not inside[bin_index]:
                bin_index += 1
                continue
            right_bin = bin_index
            while right_bin + 1 < n_bins and inside[right_bin + 1]:
                right_bin += 1
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
            bin_index = right_bin + 1

        annotated_parts.append(contig_df)

    annotated_windows_df = pd.concat(annotated_parts, ignore_index=True)
    regions_df = pd.DataFrame(
        region_rows,
        columns=[
            "genome_id", "contig_id", "start", "end", "n_windows",
            "region_len_bp", "mean_prob", "max_prob", "boundary_kind",
        ],
    )
    return annotated_windows_df, regions_df


def seed_support(regions_df: pd.DataFrame, per_seed_regions: List[pd.DataFrame]) -> List[int]:
    """Per call, how many checkpoints make an overlapping call on their own.

    docs/error-analysis.md: with the `max` ensemble any one seed can create a
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


def apply_seed_vote(regions_df: pd.DataFrame, per_seed_regions: List[pd.DataFrame],
                    min_support: int) -> Tuple[pd.DataFrame, int]:
    """Keep calls that at least `min_support` checkpoints make on their own.

    The call-level vote of docs/error-analysis.md: extent comes from the
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
    Snapping comes first, so `n_frac` describes the reported coordinates.
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

    meta = {"edge_snap": edge_meta, "n_gap_filter": gap_meta}
    if per_seed:
        seed_regs = [pd.concat([x for x in r if len(x)], ignore_index=True) if any(len(x) for x in r)
                     else pd.DataFrame(columns=["contig_id", "start", "end"]) for r in per_seed]
        calls = calls.copy()
        calls["seed_support"] = seed_support(calls, seed_regs) if len(calls) else []
        if len(filtered):
            filtered = filtered.assign(seed_support=seed_support(filtered, seed_regs))
        meta["seed_support"] = {
            "n_checkpoints": len(per_seed),
            "min_seed_support": rules.min_seed_support,
            "requested_min_seed_support": (rules.min_seed_support if requested_min_seed_support is None
                                           else requested_min_seed_support),
            "calls_dropped_by_vote": int((filtered["filter"].str.startswith(FILTER_SEED_VOTE)).sum())
            if len(filtered) else 0,
            "calls_by_support": {str(k): int(v) for k, v in
                                 calls["seed_support"].value_counts().sort_index().items()},
        }
    return FinalCalls(calls, filtered, meta)
