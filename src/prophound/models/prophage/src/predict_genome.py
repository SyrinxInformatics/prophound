from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd
import torch
from pyfaidx import Fasta
from torch.utils.data import DataLoader, Dataset
from tqdm import tqdm

from .io_checks import package_version, sha256, validate_fasta, write_bed, write_gff3, write_region_fasta
from .multiscale import plan_fine_windows, plan_neighborhoods, refine_regions
from .train_prophage import (
    DataConfig,
    ProphageTrainConfig,
    build_model,
    collate_batch,
    load_cfg,
    onehot_encode,
    resolve_device,
    stable_sigmoid,
)


# Fraction of the way from local background to a region's plateau score at
# which the true (geNomad-called) boundary tends to sit -- see
# `build_region_table_multi_tile`. Calibrated 2026-09-08 from 10 boundary
# crossings (5 regions x 2 edges) in one genome (E. coli K-12 MG1655):
# mean 0.79, median 0.77, std 0.11. A small sample from a single genome --
# treat as a reasonable starting point, not a settled constant. Recalibrate
# against more genomes before relying on it heavily.
RAMP_MIDPOINT_FRACTION = 0.77

# Seed threshold and call-level seed vote per mode. Re-derived 2026-09-29 (user
# decision) from the vote frontier (docs/error-analysis.md "Seed-vote
# confirmation"): an ensemble call is kept only if at least MODE_MIN_SEED_SUPPORT
# checkpoints make an overlapping call of their own at the same thresholds.
# Ensemble covered50 / calls per genome, 42-genome tuning sample and the 18
# (docs/research-integrity.md C1: neither is held out):
#   exploratory 0.90, keep all (flagged) -> 0.994 / 35.6 and 1.000 / 39.6
#   default     0.95, 3 of 3            -> 0.972 / 14.1 and 1.000 / 13.3
#   conservative 0.995, 3 of 3          -> 0.944 /  8.5 and 0.962 /  6.4
# Default was 0.98 / 3 of 3 for part of 2026-09-29; moved to 0.95 the same day
# (user) because 0.98 lost small and geNomad-negative elements that 0.95 keeps
# (docs/error-analysis.md "Cost on evidence that isn't geNomad agreement").
# Exploratory then moved 0.95 -> 0.90 (user) so it keeps a reach beyond Default:
# recall saturates at 0.90 (lower thresholds add calls, not elements).
# Before 2026-09-29: 0.95 / 0.995 / 0.999 with no vote (reproduce with
# --threshold and --min-seed-support 1).
MODE_THRESHOLDS = {
    "exploratory": 0.90,
    "default": 0.95,
    "conservative": 0.995,
}
MODE_MIN_SEED_SUPPORT = {
    "exploratory": 1,
    "default": 3,
    "conservative": 3,
}


AMBIGUITY_POLICIES = ("n-channel", "uniform", "impute")
# Default flipped to impute 2026-09-24 (user decision): n-channel silently zeroes
# scores on ambiguity codes; impute is in-distribution and a no-op on clean
# sequence (816446: covered recall identical to 3 dp on the 42-genome sample).
DEFAULT_AMBIGUITY_POLICY = "impute"
# Contigs shorter than --window-len (docs/adaptive-multiscale-detection.md Part 4).
# 'scan' (default since 2026-09-28, user decision) tiles each with the largest
# scale that fits; 'skip' leaves them unscanned. Either way they are reported.
SHORT_CONTIG_MODES = ("skip", "scan")
DEFAULT_SHORT_CONTIGS = "scan"
SHORT_CONTIG_SCALES = (5000, 3000)   # below ~3 kb a window cannot judge (spike-in floor)
SHORT_CONTIG_STRIDE = 500
# Seed/extend threshold per mode and short-contig scale. Calibrated 2026-09-28
# (818235, scripts/short_contig_sweep.py): the 42-genome size-stratified sample
# shredded into 27,719 fragments of 3-10 kb; each threshold is the one whose false
# calls per Mb of element-free sequence match the intact-genome 10 kb pipeline in
# the same mode (shipped merging). Ensemble, short vs 10 kb reference (false/Mb):
# exploratory 0.95 -> 6.55 vs 6.04; default 0.99 -> 2.64 vs 2.52; conservative
# 0.998 -> 1.13 vs 1.33 (0.995 overshoots at 1.84). Single seeds agree. The sweep
# applied one threshold to both scales, so they share a value; a
# fragment-length-normalised split puts 5 kb ~13% above the reference and 3 kb
# below it, inside the reference's own seed spread.
# Re-derived 2026-09-29 under the seed vote (822053, short_contig_vote_sweep.py;
# rule "never exceed", user): the lowest threshold whose shredded false calls/Mb,
# with the mode's vote, do not exceed the intact 10 kb pipeline's in that mode.
#   default 0.95/3-of-3: shredded 2.06 vs reference 2.45 (0.90 gives 3.12)
#   conservative 0.995/3-of-3: shredded 0.84 at 0.99 vs 1.09 (0.985 gives 1.11)
#   exploratory 0.90/keep all: shredded 6.55 at 0.95 vs reference 7.60 (0.90 gives 9.68)
SHORT_CONTIG_THRESHOLDS = {
    "exploratory": {5000: 0.95, 3000: 0.95},
    "default": {5000: 0.95, 3000: 0.95},
    "conservative": {5000: 0.99, 3000: 0.99},
}
REGION_COLUMNS = ["genome_id", "contig_id", "start", "end", "n_windows",
                  "region_len_bp", "mean_prob", "max_prob", "boundary_kind"]


def encode_window(sequence: str, policy: str = DEFAULT_AMBIGUITY_POLICY) -> torch.Tensor:
    """One-hot encode a window, choosing how IUPAC ambiguity codes are represented.

    `n-channel` is the training-time encoding: every non-ACGT base sets the
    dedicated N channel. It was the default until 2026-09-24 and reproduces
    every number published before then; `impute` is now the default.

    Measured 2026-09-22: under `n-channel`, two contiguous ambiguity codes take
    a confident window from p=1.000 to exactly 0.000, and a real 37 kb prophage
    in EDL933 scores ~0 because of 2-7 stray IUPAC bases per window. Both
    alternatives below restore p~1.000 there and are exact no-ops on
    ambiguity-free sequence. See docs/scientific-validity.md "1a".

    `impute` replaces each ambiguous base with a concrete A/C/G/T, chosen
    deterministically from the window's own content so repeated runs agree.
    This is the only option that keeps the input *in distribution*: the model
    was trained on one-hot vectors, and the N channel is effectively untrained
    (ambiguity codes appear in ~0.2% of training windows, equally in both
    classes, and the first-layer N-channel weights stay ~2.1x larger than the
    ACGT ones across all three seeds -- i.e. near initialization). Under max
    pooling a single such activation can dominate the pooled embedding, which
    is why one or two bases suffice.

    `uniform` instead spreads each ambiguous base as 0.25 across A/C/G/T. It is
    deterministic and invents no specific base, but fractional inputs are
    themselves out of distribution, so prefer `impute` unless measurement says
    otherwise.
    """
    x = onehot_encode(sequence)
    if policy == "n-channel":
        return x
    if policy not in AMBIGUITY_POLICIES:
        raise ValueError(f"unknown ambiguity policy {policy!r}; expected one of {AMBIGUITY_POLICIES}")
    ambiguous = x[4] > 0
    if not bool(ambiguous.any()):
        return x
    x[4][ambiguous] = 0.0
    if policy == "uniform":
        x[:4, ambiguous] = 0.25
        return x
    # impute: deterministic per window -- seed from the unambiguous content so
    # the same window always yields the same substitution, across runs and hosts
    idx = torch.nonzero(ambiguous, as_tuple=False).flatten()
    seed = (int(x[:4].sum().item()) * 1000003 + int(idx.sum().item())) & 0x7FFFFFFF
    gen = torch.Generator().manual_seed(seed)
    picks = torch.randint(0, 4, (idx.numel(),), generator=gen)
    x[picks, idx] = 1.0
    return x


class GenomeWindowDataset(Dataset):
    """Window dataset for direct inference on a FASTA file."""

    def __init__(self, windows_df: pd.DataFrame, window_len: int, ambiguity_policy: str = DEFAULT_AMBIGUITY_POLICY):
        self.windows_df = windows_df.reset_index(drop=True)
        self.window_len = window_len
        self.ambiguity_policy = ambiguity_policy
        fasta_paths = self.windows_df["fasta_path"].drop_duplicates().tolist()
        if len(fasta_paths) != 1:
            raise ValueError("GenomeWindowDataset expects exactly one FASTA path.")
        self.fasta_path = fasta_paths[0]
        self.fasta = Fasta(self.fasta_path, as_raw=True, sequence_always_upper=True)

    def __len__(self) -> int:
        return len(self.windows_df)

    def __getitem__(self, index: int) -> Tuple[torch.Tensor, torch.Tensor]:
        row = self.windows_df.iloc[index]
        contig_id = row["contig_id"]
        start = int(row["start"])
        end = int(row["end"])
        sequence = self.fasta[contig_id][start:end]
        if len(sequence) < self.window_len:
            # Unreachable via `build_windows`, which emits only full-length
            # windows -- but padding with N would silently zero the score
            # (ambiguity codes saturate the model; see
            # docs/scientific-validity.md "1a"), so fail rather than pad.
            raise RuntimeError(
                f"short window {contig_id}:{start}-{end} ({len(sequence)} < {self.window_len} bp); "
                "padding with N would silently zero its score"
            )
        elif len(sequence) > self.window_len:
            sequence = sequence[: self.window_len]
        return encode_window(sequence, self.ambiguity_policy), torch.tensor(0.0, dtype=torch.float32)

    def close(self) -> None:
        try:
            self.fasta.close()
        except Exception:
            pass

    def __del__(self):
        self.close()


@torch.no_grad()
def predict_probs(
    model: torch.nn.Module,
    loader: DataLoader,
    device: torch.device,
) -> np.ndarray:
    model.eval()
    all_logits: List[np.ndarray] = []
    for batch_x, _ in tqdm(loader, desc="predict", leave=False):
        batch_x = batch_x.to(device, non_blocking=True)
        output = model(batch_x)
        # A checkpoint trained with model.predict_boundary=True returns a
        # (classification_logit, boundary_raw) tuple; the boundary head isn't
        # wired into inference yet (see docs/boundary-refinement.md), so only
        # the classification output is used here for now.
        logits = (output[0] if isinstance(output, tuple) else output).view(-1)
        all_logits.append(logits.detach().cpu().numpy())
    if not all_logits:
        return np.array([], dtype=np.float32)
    return stable_sigmoid(np.concatenate(all_logits))


def iter_fasta_lengths(fasta_path: str) -> Iterable[Tuple[str, int]]:
    fasta = Fasta(fasta_path, as_raw=True, sequence_always_upper=True)
    try:
        for contig_id in fasta.keys():
            yield contig_id, len(fasta[contig_id])
    finally:
        fasta.close()


def build_windows(
    fasta_path: str,
    genome_id: str,
    data_config: DataConfig,
    offset: int = 0,
) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    for contig_id, contig_length in iter_fasta_lengths(fasta_path):
        if contig_length < data_config.window_len:
            continue
        for start in range(
            offset,
            contig_length - data_config.window_len + 1,
            data_config.stride,
        ):
            rows.append(
                {
                    "genome_id": genome_id,
                    "contig_id": contig_id,
                    "start": int(start),
                    "end": int(start + data_config.window_len),
                    "fasta_path": fasta_path,
                }
            )
    return pd.DataFrame(rows)


def build_short_contig_windows(
    fasta_path: str,
    genome_id: str,
    window_len: int,
    scales: Tuple[int, ...] = SHORT_CONTIG_SCALES,
    stride: int = SHORT_CONTIG_STRIDE,
) -> Tuple[pd.DataFrame, Dict[str, object]]:
    """Windows for contigs shorter than `window_len`, plus a report of them all.

    `build_windows` skips these contigs. Each one at least as long as the
    smallest scale is tiled with the largest scale that fits (column `scale`);
    shorter ones stay unscanned. The report always covers every short contig, so
    the caller can warn even when scanning is off.
    """
    rows: List[Dict[str, object]] = []
    report = {"n_short": 0, "bp_short": 0, "n_scannable": 0, "n_too_short": 0, "bp_too_short": 0}
    for contig_id, length in iter_fasta_lengths(fasta_path):
        if length >= window_len:
            continue
        report["n_short"] += 1
        report["bp_short"] += int(length)
        fitting = [sc for sc in scales if sc <= length and sc < window_len]
        if not fitting:
            report["n_too_short"] += 1
            report["bp_too_short"] += int(length)
            continue
        report["n_scannable"] += 1
        scale = max(fitting)
        for start in range(0, length - scale + 1, stride):
            rows.append({"genome_id": genome_id, "contig_id": contig_id, "start": int(start),
                         "end": int(start + scale), "fasta_path": fasta_path, "scale": int(scale)})
    cols = ["genome_id", "contig_id", "start", "end", "fasta_path", "scale"]
    return pd.DataFrame(rows, columns=cols), report


def build_multi_tile_windows(
    fasta_path: str,
    genome_id: str,
    data_config: DataConfig,
    num_tiles: int,
    grid_offset: int = 0,
) -> pd.DataFrame:
    """Generate windows from `num_tiles` phase-shifted grids over the same stride.

    A single grid (num_tiles=1, i.e. `build_windows`) can only place a region
    boundary as precisely as `data_config.stride`. Scoring several offset
    copies of the grid and merging them (see `build_region_table_multi_tile`)
    lets boundaries land at a finer effective resolution of
    `stride // num_tiles`, at the cost of ~`num_tiles`x more model calls.

    This is opt-in (see POST_V1_BACKLOG.md's "randomized windowing" idea) and
    must never change output when num_tiles <= 1.
    """
    if num_tiles <= 1:
        return build_windows(fasta_path, genome_id, data_config, offset=grid_offset)

    tile_stride = max(1, data_config.stride // num_tiles)
    tiles = []
    for tile_index in range(num_tiles):
        tile_df = build_windows(
            fasta_path, genome_id, data_config, offset=grid_offset + tile_index * tile_stride
        )
        tile_df["tile"] = tile_index
        tiles.append(tile_df)
    return pd.concat(tiles, ignore_index=True)


def build_hysteresis_windows(
    prob_values: List[float],
    seed_threshold: float,
    extend_threshold: float,
) -> List[int]:
    if not prob_values:
        return []

    positive_values = [0] * len(prob_values)
    for seed_index, seed_prob in enumerate(prob_values):
        if seed_prob < seed_threshold or positive_values[seed_index] == 1:
            continue

        left_index = seed_index
        while left_index - 1 >= 0 and prob_values[left_index - 1] >= extend_threshold:
            left_index -= 1

        right_index = seed_index
        while (
            right_index + 1 < len(prob_values)
            and prob_values[right_index + 1] >= extend_threshold
        ):
            right_index += 1

        for region_index in range(left_index, right_index + 1):
            positive_values[region_index] = 1

    return positive_values


def build_region_table(
    windows_df: pd.DataFrame,
    seed_threshold: float,
    extend_threshold: float,
) -> pd.DataFrame:
    annotated_parts: List[pd.DataFrame] = []
    region_rows: List[Dict[str, object]] = []

    for (genome_id, contig_id), contig_df in windows_df.groupby(
        ["genome_id", "contig_id"],
        sort=False,
    ):
        contig_df = contig_df.sort_values(["start", "end"]).reset_index(drop=True).copy()
        contig_df["region_pred_label"] = build_hysteresis_windows(
            contig_df["prob"].astype(float).to_list(),
            seed_threshold=seed_threshold,
            extend_threshold=extend_threshold,
        )
        annotated_parts.append(contig_df)

        stride_values = contig_df["start"].astype(int).diff().dropna()
        stride = int(stride_values.mode().iloc[0]) if not stride_values.empty else int(
            contig_df["end"].iloc[0] - contig_df["start"].iloc[0]
        )
        window_len = int(contig_df["end"].iloc[0] - contig_df["start"].iloc[0])
        contig_end = int(contig_df["end"].max())
        n_bins = int(math.ceil(contig_end / stride))
        bins_per_window = max(1, int(math.ceil(window_len / stride)))

        bin_score_sum = np.zeros(n_bins, dtype=np.float64)
        bin_score_count = np.zeros(n_bins, dtype=np.int32)
        seed_bin_mask = np.zeros(n_bins, dtype=bool)

        for row in contig_df.itertuples(index=False):
            start_bin = int(row.start // stride)
            end_bin = min(n_bins, start_bin + bins_per_window)
            bin_score_sum[start_bin:end_bin] += float(row.prob)
            bin_score_count[start_bin:end_bin] += 1
            if float(row.prob) >= seed_threshold:
                seed_bin_mask[start_bin:end_bin] = True

        bin_scores = np.zeros(n_bins, dtype=np.float64)
        covered_bins = bin_score_count > 0
        bin_scores[covered_bins] = bin_score_sum[covered_bins] / bin_score_count[covered_bins]
        extend_bin_mask = covered_bins & (bin_scores >= extend_threshold)

        region_bin_ranges: List[Tuple[int, int]] = []
        used_bins = np.zeros(n_bins, dtype=bool)
        for bin_index in range(n_bins):
            if used_bins[bin_index] or not seed_bin_mask[bin_index] or not extend_bin_mask[bin_index]:
                continue

            left_bin = bin_index
            while left_bin - 1 >= 0 and extend_bin_mask[left_bin - 1]:
                left_bin -= 1

            right_bin = bin_index
            while right_bin + 1 < n_bins and extend_bin_mask[right_bin + 1]:
                right_bin += 1

            used_bins[left_bin : right_bin + 1] = True
            region_bin_ranges.append((left_bin, right_bin))

        previous_end = 0
        for left_bin, right_bin in region_bin_ranges:
            refined_start = left_bin * stride
            refined_end = min(contig_end, (right_bin + 1) * stride)

            if left_bin > 0:
                left_score = bin_scores[left_bin - 1]
                right_score = bin_scores[left_bin]
                if right_score > left_score and left_score < extend_threshold <= right_score:
                    frac = (extend_threshold - left_score) / (right_score - left_score)
                    refined_start = int(round(((left_bin - 1) + frac) * stride))

            if right_bin + 1 < n_bins:
                left_score = bin_scores[right_bin]
                right_score = bin_scores[right_bin + 1]
                if left_score > right_score and right_score < extend_threshold <= left_score:
                    frac = (left_score - extend_threshold) / (left_score - right_score)
                    refined_end = int(round((right_bin + 1 - frac) * stride))

            refined_start = max(previous_end, max(0, refined_start))
            refined_end = max(refined_start + 1, min(contig_end, refined_end))
            previous_end = refined_end

            overlapping_windows = contig_df[
                (contig_df["end"].astype(int) > refined_start)
                & (contig_df["start"].astype(int) < refined_end)
                & (contig_df["region_pred_label"].astype(int) == 1)
            ].copy()
            if overlapping_windows.empty:
                overlapping_windows = contig_df[
                    (contig_df["start"].astype(int) // stride <= right_bin)
                    & (contig_df["end"].astype(int) // stride >= left_bin + 1)
                    & (contig_df["region_pred_label"].astype(int) == 1)
                ].copy()
            if overlapping_windows.empty:
                continue

            region_rows.append(
                {
                    "genome_id": genome_id,
                    "contig_id": contig_id,
                    "start": int(refined_start),
                    "end": int(refined_end),
                    "n_windows": int(len(overlapping_windows)),
                    "region_len_bp": int(refined_end - refined_start),
                    "mean_prob": float(overlapping_windows["prob"].astype(float).mean()),
                    "max_prob": float(overlapping_windows["prob"].astype(float).max()),
                    "boundary_kind": "non_overlapping_base_estimate",
                }
            )

    annotated_windows_df = pd.concat(annotated_parts, ignore_index=True)
    return annotated_windows_df, pd.DataFrame(region_rows)


REGION_BUILDERS = ("strict", "coverage")
ENSEMBLE_MODES = ("max", "mean")
DEFAULT_ENSEMBLE = "max"
# Defaults changed 2026-09-21 to the measured operating point for the ">=95%
# recall on known regions" contract -- docs/model-quality-baseline.md
# "Operating point". `strict` remains fully supported as the conservative,
# database-agreement path; `coverage` + merging is the discovery path.
DEFAULT_REGION_BUILDER = "coverage"
DEFAULT_COVERAGE_FRAC = 0.1
# Fragment merging (docs/gap-filling.md). 0 = off, which is the shipped default
# until a held-out sweep picks a value -- see `_close_gaps`.
DEFAULT_MAX_GAP_BP = 50000
DEFAULT_BRIDGE_PROB = 0.1
# Fine-scale refinement, pass B (multiscale.py; docs/adaptive-multiscale-detection.md).
# Off by default: 0 = no refinement, and with it off output is byte-identical.
DEFAULT_REFINE_WINDOW_LEN = 0
DEFAULT_REFINE_STRIDE = 250
# Fine-pass threshold for the opt-in "tight boundaries" setting (5 kb, alpha 0.6,
# ramp 15000 60000). Swept 2026-09-28 on the 42-genome tracks at 5 kb (ensemble):
# 0.95 gives median IoU 0.785 at covered50 0.955 (vs 0.713 / 0.961 unrefined); the
# previous fallback, the mode's extend threshold (0.99), gave 0.793 / 0.916. 0.6-0.8
# is recall-neutral (+0.04-0.06 IoU) and is being confirmed on the 18-genome
# sample. Default mode only; docs/downstream-postprocessing.md.
DEFAULT_REFINE_THRESHOLD = 0.95
DEFAULT_EDGE_ALPHA = 1.0        # 1 = keep coarse edges; calibrate before lowering
# Split gates as swept for the tight-boundaries setting (sweep_ramp, 2026-09-24):
# 4 kb / 0.05 adds ~0.3-0.9 fragmented elements; 0.2 added ~4.7. Were 2000 / 0.1,
# which no ramp sweep measured. Below DEFAULT_BRIDGE_PROB, so splitting cannot undo
# a merge; 0.05-0.1 is a dead band where neither acts.
DEFAULT_SPLIT_MIN_GAP_BP = 4000
DEFAULT_SPLIT_MAX_PROB = 0.05


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

    `build_region_table` marks a stride-sized bin as inside a region only if
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
    contract as `build_region_table`; `region_pred_label` marks windows that
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


def _walk_to_ramp_crossing(
    bin_scores: np.ndarray,
    covered_bins: np.ndarray,
    from_bin: int,
    direction: int,
    crossing_threshold: float,
    max_steps: int,
    n_bins: int,
) -> int:
    """Walk outward (direction=-1 left, +1 right) from `from_bin` while the
    score stays at or above `crossing_threshold`, and return the last such bin
    index.

    Used to locate where a region's score ramp crosses a calibrated fraction
    of the way between its plateau and local background
    (`RAMP_MIDPOINT_FRACTION`), rather than a single fixed threshold -- see
    `build_region_table_multi_tile` docstring for why this tracks true
    boundaries much better than a strict high-confidence threshold does.
    """
    bin_index = from_bin
    steps = 0
    while steps < max_steps:
        next_index = bin_index + direction
        if next_index < 0 or next_index >= n_bins:
            break
        if not covered_bins[next_index] or bin_scores[next_index] < crossing_threshold:
            break
        bin_index = next_index
        steps += 1
    return bin_index


def build_region_table_multi_tile(
    windows_df: pd.DataFrame,
    bin_size: int,
    seed_threshold: float,
    extend_threshold: float,
    window_len: int,
    stride: int,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Hybrid region construction: seed on the trusted single grid, refine
    boundaries on the finer multi-tile score track.

    Two earlier approaches were tried and both failed at scale (see
    docs/boundary-refinement.md and docs/model-quality-baseline.md):
    - Seeding directly from the *mean* multi-tile bin score under-detects:
      more tiles means more overlapping windows averaged into each bin
      (~10 per tile), which dilutes a single well-aligned strong window's
      score toward its worse-aligned neighbors. Measured region recall
      dropped from 76.9% (single-grid) to 67.3% (mean multi-tile) on an
      18-genome eval.
    - Seeding from the *max* multi-tile bin score over-detects badly: with
      ~4x more overlapping windows contributing to each bin, the chance that
      at least one of them exceeds a strict threshold by pure noise goes up
      accordingly (classic multiple-comparisons peak inflation). This
      hallucinated a dozen spurious regions in a genome with exactly one true
      region in a quick 3-genome check.

    So seeding and boundary refinement are deliberately decoupled here:
    - WHICH regions exist is decided using only the tile-0 windows (`tile`
      column == 0), replicating `build_region_table`'s seeding rule exactly:
      a bin is a seed if ANY single covering window clears `seed_threshold`,
      and extends while the bin MEAN clears `extend_threshold`. This carries
      over the frozen v1 detection behavior with no new false positives and
      no new dilution -- seeding is not where the finer tiling should help
      anyway. (Until 2026-09-15 this tested the bin mean against
      `seed_threshold`, which is stricter than the frozen path and dropped
      regions -- see docs/boundary-refinement.md "Discrepancy resolved".)
    - HOW PRECISELY each already-seeded region's boundaries are drawn uses
      the full multi-tile mean-based fine-grained bin track (`bin_size` =
      `stride / num_tiles`), walking outward from the seed to where the score
      crosses a calibrated fraction (`RAMP_MIDPOINT_FRACTION`) of the way from
      local background up to the region's plateau, rather than where it
      crosses `extend_threshold`.

    Why a ramp-crossing fraction instead of `extend_threshold` at all: the
    model's score does not jump sharply at a true biological boundary -- a
    window straddling the edge sees a mix of positive/negative sequence and
    scores roughly in proportion, so the score ramps smoothly over about one
    `window_len` on each side. Requiring the score to stay above a strict
    threshold like 0.99 before extending only claims the saturated middle of
    that ramp, landing the boundary roughly `window_len / 2` inside the true
    edge. A naive 50%-of-the-way midpoint overcorrects the other way (every
    region ends up bracketing the true call instead of sitting inside it).
    `RAMP_MIDPOINT_FRACTION` (0.77) was calibrated from 10 boundary crossings
    in one genome (E. coli K-12) -- a small sample, not a proven constant.

    Also returns a per-bin diagnostic track (`bins_df`) with a `bin_stability`
    heuristic (1 - stdev of per-window scores landing in that bin): closer to
    1 means tiles agreed on that position's score, closer to 0 means they
    disagreed. This is a heuristic, not a calibrated confidence interval.
    """
    region_rows: List[Dict[str, object]] = []
    bin_frames: List[pd.DataFrame] = []

    for (genome_id, contig_id), contig_df in windows_df.groupby(
        ["genome_id", "contig_id"], sort=False
    ):
        contig_df = contig_df.sort_values(["start", "end"]).reset_index(drop=True)
        contig_end = int(contig_df["end"].max())
        n_bins = int(math.ceil(contig_end / bin_size))

        # --- Fine-grained mean-based track from ALL tiles (refinement only) ---
        bin_score_sum = np.zeros(n_bins, dtype=np.float64)
        bin_score_sq_sum = np.zeros(n_bins, dtype=np.float64)
        bin_score_count = np.zeros(n_bins, dtype=np.int32)

        for row in contig_df.itertuples(index=False):
            start_bin = int(row.start // bin_size)
            end_bin = min(n_bins, int(math.ceil(row.end / bin_size)))
            prob = float(row.prob)
            bin_score_sum[start_bin:end_bin] += prob
            bin_score_sq_sum[start_bin:end_bin] += prob * prob
            bin_score_count[start_bin:end_bin] += 1

        covered_bins = bin_score_count > 0
        bin_scores = np.zeros(n_bins, dtype=np.float64)
        bin_scores[covered_bins] = bin_score_sum[covered_bins] / bin_score_count[covered_bins]
        bin_variance = np.zeros(n_bins, dtype=np.float64)
        bin_variance[covered_bins] = np.clip(
            bin_score_sq_sum[covered_bins] / bin_score_count[covered_bins]
            - bin_scores[covered_bins] ** 2,
            0.0,
            None,
        )
        bin_stability = np.full(n_bins, np.nan, dtype=np.float64)
        bin_stability[covered_bins] = 1.0 - np.sqrt(bin_variance[covered_bins])

        # --- Coarse single-grid (tile 0) track: decides WHICH regions exist ---
        tile0_df = contig_df[contig_df["tile"] == 0] if "tile" in contig_df.columns else contig_df
        coarse_n_bins = int(math.ceil(contig_end / stride))
        coarse_sum = np.zeros(coarse_n_bins, dtype=np.float64)
        coarse_count = np.zeros(coarse_n_bins, dtype=np.int32)
        # Seeding must mirror `build_region_table` exactly: a bin is a seed
        # candidate if ANY single window covering it clears `seed_threshold`
        # (per-window max), while extension uses the bin MEAN against
        # `extend_threshold`. An earlier version tested the bin mean against
        # `seed_threshold` too, which is strictly stricter and silently
        # dropped regions the frozen single-tile path emits (one 0.996 window
        # among 0.99x neighbours seeds there but never averaged to 0.995).
        coarse_seed_mask = np.zeros(coarse_n_bins, dtype=bool)
        for row in tile0_df.itertuples(index=False):
            cb_start = int(row.start // stride)
            cb_end = min(coarse_n_bins, int(math.ceil(row.end / stride)))
            prob = float(row.prob)
            coarse_sum[cb_start:cb_end] += prob
            coarse_count[cb_start:cb_end] += 1
            if prob >= seed_threshold:
                coarse_seed_mask[cb_start:cb_end] = True
        coarse_covered = coarse_count > 0
        coarse_scores = np.zeros(coarse_n_bins, dtype=np.float64)
        coarse_scores[coarse_covered] = coarse_sum[coarse_covered] / coarse_count[coarse_covered]

        coarse_extend_mask = coarse_covered & (coarse_scores >= extend_threshold)

        seed_region_ranges_bp: List[Tuple[int, int]] = []
        used_coarse = np.zeros(coarse_n_bins, dtype=bool)
        for coarse_index in range(coarse_n_bins):
            if (
                used_coarse[coarse_index]
                or not coarse_seed_mask[coarse_index]
                or not coarse_extend_mask[coarse_index]
            ):
                continue
            left = coarse_index
            while left - 1 >= 0 and coarse_extend_mask[left - 1]:
                left -= 1
            right = coarse_index
            while right + 1 < coarse_n_bins and coarse_extend_mask[right + 1]:
                right += 1
            used_coarse[left : right + 1] = True
            seed_region_ranges_bp.append((left * stride, min(contig_end, (right + 1) * stride)))

        # Background for the ramp-crossing fraction: fine bins outside every
        # seeded region (not "not seeded" at the fine level, since seeding now
        # happens on the coarse grid).
        in_any_seed_region = np.zeros(n_bins, dtype=bool)
        for seed_start_bp, seed_end_bp in seed_region_ranges_bp:
            in_any_seed_region[int(seed_start_bp // bin_size) : int(math.ceil(seed_end_bp / bin_size))] = True
        background_candidates = bin_scores[covered_bins & ~in_any_seed_region]
        background_level = (
            float(np.median(background_candidates)) if background_candidates.size else 0.0
        )
        # Cap the ramp walk generously past one window's width, since the
        # score ramp around a true edge is expected to span roughly window_len.
        max_walk_steps = max(1, int(math.ceil((2 * window_len) / bin_size)))

        previous_end = 0
        for seed_start_bp, seed_end_bp in seed_region_ranges_bp:
            left_bin = int(seed_start_bp // bin_size)
            right_bin = min(n_bins - 1, int(math.ceil(seed_end_bp / bin_size)) - 1)

            plateau = float(bin_scores[left_bin : right_bin + 1].max())
            crossing_threshold = background_level + RAMP_MIDPOINT_FRACTION * (
                plateau - background_level
            )
            # Guard against a degenerate plateau/background (e.g. a very weak
            # region): never target below the extend_threshold's own strict
            # floor, so this can only push boundaries outward, not inward.
            crossing_threshold = (
                min(crossing_threshold, extend_threshold)
                if plateau > background_level
                else extend_threshold
            )

            left_edge_bin = _walk_to_ramp_crossing(
                bin_scores, covered_bins, left_bin, -1, crossing_threshold, max_walk_steps, n_bins
            )
            right_edge_bin = _walk_to_ramp_crossing(
                bin_scores, covered_bins, right_bin, 1, crossing_threshold, max_walk_steps, n_bins
            )

            refined_start = left_edge_bin * bin_size
            refined_end = min(contig_end, (right_edge_bin + 1) * bin_size)

            if left_edge_bin > 0 and covered_bins[left_edge_bin - 1]:
                left_score = bin_scores[left_edge_bin - 1]
                right_score = bin_scores[left_edge_bin]
                if right_score > left_score and left_score < crossing_threshold <= right_score:
                    frac = (crossing_threshold - left_score) / (right_score - left_score)
                    refined_start = int(round(((left_edge_bin - 1) + frac) * bin_size))

            if right_edge_bin + 1 < n_bins and covered_bins[right_edge_bin + 1]:
                left_score = bin_scores[right_edge_bin]
                right_score = bin_scores[right_edge_bin + 1]
                if left_score > right_score and right_score < crossing_threshold <= left_score:
                    frac = (left_score - crossing_threshold) / (left_score - right_score)
                    refined_end = int(round((right_edge_bin + 1 - frac) * bin_size))

            refined_start = max(previous_end, max(0, refined_start))
            refined_end = max(refined_start + 1, min(contig_end, refined_end))
            previous_end = refined_end

            region_bin_slice = slice(left_edge_bin, right_edge_bin + 1)
            region_bin_covered = covered_bins[region_bin_slice]
            region_stability = bin_stability[region_bin_slice][region_bin_covered]

            region_rows.append(
                {
                    "genome_id": genome_id,
                    "contig_id": contig_id,
                    "start": int(refined_start),
                    "end": int(refined_end),
                    "region_len_bp": int(refined_end - refined_start),
                    "n_windows": int(
                        (
                            (contig_df["end"].astype(int) > refined_start)
                            & (contig_df["start"].astype(int) < refined_end)
                        ).sum()
                    ),
                    "mean_prob": float(bin_scores[region_bin_slice].mean()),
                    "max_prob": float(bin_scores[region_bin_slice].max()),
                    "mean_stability": float(region_stability.mean())
                    if region_stability.size
                    else float("nan"),
                    "boundary_kind": "multi_tile_midpoint_estimate",
                }
            )

        bin_frames.append(
            pd.DataFrame(
                {
                    "genome_id": genome_id,
                    "contig_id": contig_id,
                    "bin_start": np.arange(n_bins) * bin_size,
                    "bin_end": np.minimum((np.arange(n_bins) + 1) * bin_size, contig_end),
                    "bin_score": bin_scores,
                    "bin_stability": bin_stability,
                    "bin_n_windows": bin_score_count,
                }
            )
        )

    regions_df = pd.DataFrame(region_rows)
    bins_df = pd.concat(bin_frames, ignore_index=True) if bin_frames else pd.DataFrame()
    return bins_df, regions_df


def resolve_thresholds(
    mode: str,
    threshold: float | None,
    seed_threshold: float | None,
    extend_threshold: float | None,
) -> Tuple[float, float, float]:
    mode_threshold = MODE_THRESHOLDS[mode]
    score_threshold = threshold if threshold is not None else mode_threshold
    region_seed_threshold = seed_threshold if seed_threshold is not None else score_threshold
    region_extend_threshold = (
        extend_threshold if extend_threshold is not None else min(region_seed_threshold, 0.99)
    )
    return score_threshold, region_seed_threshold, region_extend_threshold


def write_metadata(
    cfg: ProphageTrainConfig,
    checkpoint_path,
    fasta_path: str,
    mode: str,
    threshold: float,
    seed_threshold: float,
    extend_threshold: float,
    device: str,
    windows_df: pd.DataFrame,
    regions_df: pd.DataFrame,
    out_path: Path,
    num_tiles: int = 1,
    ambiguity_policy: str = DEFAULT_AMBIGUITY_POLICY,
    grid_offset: int = 0,
    region_builder: str = DEFAULT_REGION_BUILDER,
    coverage_frac: float = DEFAULT_COVERAGE_FRAC,
    ensemble_mode: str = "single",
    max_gap_bp: int = DEFAULT_MAX_GAP_BP,
    bridge_prob: float = DEFAULT_BRIDGE_PROB,
    trained_window_len: int | None = None,
    trained_stride: int | None = None,
    refine: dict | None = None,
    short_contigs: dict | None = None,
    seed_support: dict | None = None,
    extra: dict | None = None,
) -> None:
    metadata = {
        "prophound_version": package_version(),
        "run_name": cfg.run_name,
        "checkpoint": checkpoint_path,
        "ensemble": ensemble_mode,
        "input_fasta": fasta_path,
        "mode": mode,
        "threshold": threshold,
        "seed_threshold": seed_threshold,
        "extend_threshold": extend_threshold,
        "device": device,
        "window_len": cfg.data.window_len,
        "stride": cfg.data.stride,
        "trained_window_len": trained_window_len if trained_window_len is not None else cfg.data.window_len,
        "trained_stride": trained_stride if trained_stride is not None else cfg.data.stride,
        "num_tiles": num_tiles,
        "ambiguity_policy": ambiguity_policy,
        "grid_offset": grid_offset,
        "n_windows": int(len(windows_df)),
        "n_regions": int(len(regions_df)),
    }
    if extra:
        metadata.update(extra)
    if short_contigs:
        metadata["short_contigs"] = short_contigs
    if seed_support:
        metadata["seed_support"] = seed_support
    if num_tiles > 1:
        metadata["region_builder"] = "multi_tile_hysteresis"
        metadata["boundary_refinement"] = "multi_tile_ramp_midpoint_crossing"
        metadata["bin_size"] = max(1, cfg.data.stride // num_tiles)
    elif region_builder == "coverage":
        metadata["region_builder"] = "coverage_fraction"
        metadata["boundary_refinement"] = "none_bin_edges"
        metadata["coverage_frac"] = coverage_frac
        if max_gap_bp > 0:
            metadata["fragment_merge"] = {"max_gap_bp": max_gap_bp, "bridge_prob": bridge_prob}
        if refine:
            metadata["refine"] = refine
    else:
        metadata["region_builder"] = "strict_hysteresis_stride_support"
        metadata["boundary_refinement"] = "interpolated_from_stride_support_track"
    out_path.write_text(json.dumps(metadata, indent=2) + "\n")


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run frozen MGE-Explorer inference on a genome FASTA."
    )
    parser.add_argument("--config", required=True)
    parser.add_argument(
        "--checkpoint",
        required=True,
        nargs="+",
        help=(
            "Trained checkpoint(s). Passing more than one runs a seed ensemble "
            "over them (see --ensemble); recommended, since it removes the "
            "dependence of recall on which replicate was packaged."
        ),
    )
    parser.add_argument(
        "--ensemble",
        choices=ENSEMBLE_MODES,
        default=DEFAULT_ENSEMBLE,
        help=(
            "How to combine multiple --checkpoint score tracks. 'max' "
            "(default) preserves borderline detections; 'mean' dilutes them "
            "below the seed threshold and measurably costs recall."
        ),
    )
    parser.add_argument("--input_fasta", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument(
        "--mode",
        choices=sorted(MODE_THRESHOLDS.keys()),
        default="default",
    )
    parser.add_argument("--threshold", type=float, default=None)
    parser.add_argument("--seed-threshold", type=float, default=None)
    parser.add_argument("--extend-threshold", type=float, default=None)
    parser.add_argument(
        "--min-seed-support",
        type=int,
        default=None,
        help=(
            "Keep an ensemble call only if at least this many checkpoints make an "
            "overlapping call on their own (default: the mode's; exploratory 1, "
            "default 3, conservative 3). Needs several --checkpoint files; 1 keeps "
            "every call. See docs/error-analysis.md."
        ),
    )
    parser.add_argument(
        "--write-region-fasta",
        action="store_true",
        help="Also write the called regions' sequences to regions.fna.",
    )
    parser.add_argument(
        "--write-seed-scores",
        action="store_true",
        help=(
            "With several --checkpoint files, also write each checkpoint's score to "
            "windows.tsv (prob_ckpt1, prob_ckpt2, ...), so seed votes can be re-applied "
            "offline. Off by default."
        ),
    )
    parser.add_argument("--device", default="auto")
    parser.add_argument("--batch-size", type=int, default=128)
    parser.add_argument("--num-workers", type=int, default=0)
    parser.add_argument(
        "--ambiguity-policy",
        choices=AMBIGUITY_POLICIES,
        default=DEFAULT_AMBIGUITY_POLICY,
        help=(
            "How IUPAC ambiguity codes are encoded. 'impute' (default) replaces each "
            "with a deterministic concrete base; a no-op on ambiguity-free sequence. "
            "'n-channel' is the training-time encoding and reproduces numbers "
            "published before 2026-09-24, but two contiguous codes can zero a "
            "window's score. 'uniform' spreads each as 0.25 across A/C/G/T. See "
            "docs/scientific-validity.md '1a'."
        ),
    )
    parser.add_argument(
        "--grid-offset",
        type=int,
        default=0,
        help="Shift the inference window grid by this many bp (default 0).",
    )
    parser.add_argument(
        "--num-tiles",
        type=int,
        default=1,
        help=(
            "EXPERIMENTAL, not recommended: score this many phase-shifted "
            "window tilings and merge them for finer region boundaries "
            "(default 1 = frozen v1 behavior, single grid). See "
            "docs/boundary-refinement.md for known issues before using this "
            "for anything but further investigation."
        ),
    )
    add_region_builder_args(parser)
    add_multiscale_args(parser)
    return parser


def add_multiscale_args(parser: argparse.ArgumentParser) -> None:
    """Window-size overrides and the fine refinement pass; mirrored in `prophound`."""
    parser.add_argument(
        "--window-len",
        type=int,
        default=None,
        help=(
            "Inference window length in bp (default: the checkpoint's trained "
            "length, 10000). The model's max pooling accepts any length, so this "
            "trades resolution for context. EXPERIMENTAL away from the trained "
            "length: mode thresholds were calibrated at it, and negative scores "
            "drift upward as windows shrink."
        ),
    )
    parser.add_argument(
        "--stride",
        type=int,
        default=None,
        help="Inference window stride in bp (default: the checkpoint's trained stride, 1000).",
    )
    parser.add_argument(
        "--short-contigs",
        choices=SHORT_CONTIG_MODES,
        default=DEFAULT_SHORT_CONTIGS,
        help=(
            "Contigs shorter than --window-len: 'scan' (default) tiles each with a 5 kb or "
            "3 kb window, whichever fits (contigs under 3 kb are never scanned); 'skip' "
            "leaves them unscanned. Either way a warning reports them. Reduced-scale "
            "thresholds are calibrated per mode to the 10 kb false-positive load "
            "(0.95 / 0.95 / 0.99, with the mode's seed vote)."
        ),
    )
    parser.add_argument(
        "--short-contig-threshold",
        type=float,
        default=None,
        help="Override the calibrated seed/extend threshold for --short-contigs scan.",
    )
    parser.add_argument(
        "--refine-window-len",
        type=int,
        default=DEFAULT_REFINE_WINDOW_LEN,
        help=(
            "EXPERIMENTAL fine-scale pass (default 0 = off). Re-score each call "
            "+/- --refine-margin at this window length, then split calls whose fine "
            "track shows a background gap (conflated elements) and place edges "
            "between the fine core and the coarse edge (--edge-alpha). Never drops "
            "a call. Requires --region-builder coverage. See "
            "docs/adaptive-multiscale-detection.md."
        ),
    )
    parser.add_argument("--refine-stride", type=int, default=DEFAULT_REFINE_STRIDE,
                        help="Fine-pass stride in bp. Default %(default)s.")
    parser.add_argument(
        "--refine-threshold",
        type=float,
        default=None,
        help=(
            "Fine-pass score threshold (default 0.95, calibrated at 5 kb for the "
            "tight-boundaries setting; 0.6-0.8 keeps recall unchanged)."
        ),
    )
    parser.add_argument("--refine-margin", type=int, default=None,
                        help="bp re-scored either side of each call (default: --window-len).")
    parser.add_argument(
        "--edge-alpha",
        type=float,
        default=DEFAULT_EDGE_ALPHA,
        help=(
            "Edge placement in [0, 1]: edge = fine_core + alpha * (coarse_edge - "
            "fine_core). 1 keeps coarse edges (overextend ~half a window); 0 uses "
            "the fine core (truncates). Default %(default)s."
        ),
    )
    parser.add_argument(
        "--edge-alpha-ramp",
        type=int,
        nargs=2,
        metavar=("LO_BP", "HI_BP"),
        default=None,
        help=(
            "Make --edge-alpha size-dependent: 0 for pieces whose fine core is at "
            "most LO_BP, --edge-alpha at HI_BP and above, linear in log length "
            "between. Small elements are best served by the fine core and large "
            "ones by edges near the coarse call. Default off (constant alpha)."
        ),
    )
    parser.add_argument(
        "--split-min-gap-bp",
        type=int,
        default=DEFAULT_SPLIT_MIN_GAP_BP,
        help="Minimum fine-track gap (bp) that can split a call. Default %(default)s.",
    )
    parser.add_argument(
        "--split-max-prob",
        type=float,
        default=DEFAULT_SPLIT_MAX_PROB,
        help=(
            "Split only across a gap whose fine windows average below this -- the "
            "mirror of --bridge-prob, so splitting cannot undo gap-filling. "
            "Default %(default)s."
        ),
    )


def add_region_builder_args(parser: argparse.ArgumentParser) -> None:
    """`--region-builder` / `--coverage-frac`, shared with `prophound`'s CLI."""
    parser.add_argument(
        "--region-builder",
        choices=REGION_BUILDERS,
        default=DEFAULT_REGION_BUILDER,
        help=(
            "How per-window scores become regions. 'coverage' (default since "
            "2026-09-21) marks a bin inside if at least --coverage-frac of its "
            "covering windows individually score >= the extend threshold; with "
            "--max-gap-bp merging it is the measured operating point for >=95% "
            "recall on known regions, and it is the discovery path. 'strict' is "
            "the frozen v1 builder: a bin is inside only if the MEAN of all "
            "windows covering it >= the extend threshold, which truncates "
            "regions and drops sparse-signal elements -- kept as the "
            "conservative, database-agreement comparison. See "
            "docs/model-quality-baseline.md 'Operating point'. 'coverage' is "
            "not combinable with --num-tiles > 1."
        ),
    )
    parser.add_argument(
        "--coverage-frac",
        type=float,
        default=DEFAULT_COVERAGE_FRAC,
        help=(
            "For --region-builder coverage: minimum fraction of covering windows "
            "that must be high for a bin to count as inside (0 -> union of high "
            "windows, 1 -> intersection). Default %(default)s."
        ),
    )
    parser.add_argument(
        "--max-gap-bp",
        type=int,
        default=DEFAULT_MAX_GAP_BP,
        help=(
            "For --region-builder coverage: merge two calls separated by at "
            "most this many bp into one, repairing elements fragmented by an "
            "internal dip in the score track (worst for 60 kb+ elements, where "
            "it lifts median IoU from ~0.38 to ~0.85). Default %(default)s; 0 "
            "disables merging. Always pair with --bridge-prob -- an ungated "
            "merge bridges genuinely separate neighbouring elements. Above "
            "~50000 the gate binds before the distance budget does. See "
            "docs/gap-filling.md."
        ),
    )
    parser.add_argument(
        "--bridge-prob",
        type=float,
        default=DEFAULT_BRIDGE_PROB,
        help=(
            "For --max-gap-bp: only merge across a gap whose every bin has mean "
            "covering-window probability >= this (0 = distance-only). The "
            "anti-overextension guard -- docs/anti-overextension.md. "
            "Default %(default)s."
        ),
    )


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


def label_windows(windows_df: pd.DataFrame, regions_df: pd.DataFrame) -> pd.DataFrame:
    """Set region_pred_label to 1 on windows overlapping any call, else 0."""
    windows_df = windows_df.copy()
    windows_df["region_pred_label"] = 0
    for r in regions_df.itertuples(index=False):
        m = ((windows_df["contig_id"] == r.contig_id)
             & (windows_df["end"].astype(int) > int(r.start)) & (windows_df["start"].astype(int) < int(r.end)))
        windows_df.loc[m, "region_pred_label"] = 1
    return windows_df


def score_windows(
    model: torch.nn.Module,
    checkpoints: List[str],
    windows_df: pd.DataFrame,
    window_len: int,
    ambiguity_policy: str,
    ensemble: str,
    device: torch.device,
    batch_size: int,
    num_workers: int,
    per_checkpoint_out: List[np.ndarray] | None = None,
) -> np.ndarray:
    """Score `windows_df` with every checkpoint and combine per `ensemble`.

    Shared by the coarse pass and the fine refinement pass so both use the same
    ensemble rule. If `per_checkpoint_out` is a list, each checkpoint's scores
    are appended to it (for seed support).
    """
    dataset = GenomeWindowDataset(windows_df, window_len=window_len, ambiguity_policy=ambiguity_policy)
    loader = DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=False,
        num_workers=num_workers,
        collate_fn=collate_batch,
        pin_memory=(device.type == "cuda"),
    )
    per_checkpoint = []
    for ckpt_path in checkpoints:
        checkpoint = torch.load(ckpt_path, map_location=device, weights_only=False)
        model.load_state_dict(checkpoint["model"])
        per_checkpoint.append(np.asarray(predict_probs(model, loader, device), dtype=np.float64))
    dataset.close()
    if per_checkpoint_out is not None:
        per_checkpoint_out.extend(per_checkpoint)

    if len(per_checkpoint) == 1:
        return per_checkpoint[0]
    stacked = np.stack(per_checkpoint)
    # max, not mean: averaging dilutes borderline seeds below the seed
    # threshold -- the same effect that sank multi-tile averaging
    # (docs/boundary-refinement.md) and mean-ensembling in the 2026-09-17
    # checks. Measured 2026-09-21: max reaches .978/.994 covered recall
    # (Default/Exploratory) against .950/.978 for mean.
    return stacked.max(axis=0) if ensemble == "max" else stacked.mean(axis=0)


def run_predict(args: argparse.Namespace) -> None:
    """Run inference end to end from a parsed argument namespace.

    Shared by this module's own CLI and by `prophound`'s `predict` subcommand,
    so both stay in sync with a single implementation.
    """
    cfg = load_cfg(args.config)
    device = resolve_device(args.device)
    threshold, seed_threshold, extend_threshold = resolve_thresholds(
        mode=args.mode,
        threshold=args.threshold,
        seed_threshold=args.seed_threshold,
        extend_threshold=args.extend_threshold,
    )

    trained_window_len, trained_stride = cfg.data.window_len, cfg.data.stride
    window_len_arg = getattr(args, "window_len", None)
    stride_arg = getattr(args, "stride", None)
    if window_len_arg:
        cfg.data.window_len = int(window_len_arg)
    if stride_arg:
        cfg.data.stride = int(stride_arg)
    if cfg.data.window_len <= 0 or cfg.data.stride <= 0:
        raise SystemExit("--window-len and --stride must be positive")
    if cfg.data.window_len != trained_window_len:
        # adaptive/global max pooling accepts any length, but score calibration
        # does not transfer: negatives drift up as windows shrink (median 0.003
        # at 10 kb -> 0.065 at 3 kb; docs/adaptive-multiscale-detection.md).
        print(
            f"[warn] --window-len {cfg.data.window_len} differs from the trained "
            f"{trained_window_len}; the mode thresholds were calibrated at "
            f"{trained_window_len} bp and are not validated at this length.",
            file=sys.stderr,
        )

    num_tiles = max(1, getattr(args, "num_tiles", 1))
    region_builder = getattr(args, "region_builder", DEFAULT_REGION_BUILDER)
    coverage_frac = float(getattr(args, "coverage_frac", DEFAULT_COVERAGE_FRAC))
    max_gap_bp = int(getattr(args, "max_gap_bp", DEFAULT_MAX_GAP_BP))
    bridge_prob = float(getattr(args, "bridge_prob", DEFAULT_BRIDGE_PROB))
    if region_builder != "strict" and num_tiles > 1:
        raise SystemExit("--region-builder coverage cannot be combined with --num-tiles > 1")
    refine_window_len = int(getattr(args, "refine_window_len", DEFAULT_REFINE_WINDOW_LEN) or 0)
    refine_stride = int(getattr(args, "refine_stride", DEFAULT_REFINE_STRIDE))
    edge_alpha = float(getattr(args, "edge_alpha", DEFAULT_EDGE_ALPHA))
    alpha_ramp = getattr(args, "edge_alpha_ramp", None)
    alpha_ramp = tuple(int(x) for x in alpha_ramp) if alpha_ramp else None
    if alpha_ramp and not 0 < alpha_ramp[0] < alpha_ramp[1]:
        raise SystemExit("--edge-alpha-ramp needs 0 < LO_BP < HI_BP")
    split_min_gap_bp = int(getattr(args, "split_min_gap_bp", DEFAULT_SPLIT_MIN_GAP_BP))
    split_max_prob = float(getattr(args, "split_max_prob", DEFAULT_SPLIT_MAX_PROB))
    refine_threshold = getattr(args, "refine_threshold", None)
    refine_threshold = DEFAULT_REFINE_THRESHOLD if refine_threshold is None else float(refine_threshold)
    refine_margin = getattr(args, "refine_margin", None)
    refine_margin = cfg.data.window_len if refine_margin is None else int(refine_margin)
    if refine_window_len:
        if region_builder != "coverage" or num_tiles > 1:
            raise SystemExit("--refine-window-len requires --region-builder coverage and --num-tiles 1")
        if refine_window_len >= cfg.data.window_len:
            raise SystemExit("--refine-window-len must be smaller than --window-len")
        if refine_stride <= 0 or not 0.0 <= edge_alpha <= 1.0 or refine_margin < 0:
            raise SystemExit("need --refine-stride > 0, 0 <= --edge-alpha <= 1, --refine-margin >= 0")

    input_path = Path(args.input_fasta)
    input_summary = validate_fasta(str(input_path))   # fails early with a clear message
    genome_id = input_path.stem
    windows_df = build_multi_tile_windows(
        str(input_path),
        genome_id=genome_id,
        data_config=cfg.data,
        num_tiles=num_tiles,
        grid_offset=max(0, int(getattr(args, "grid_offset", 0) or 0)),
    )
    short_contigs = getattr(args, "short_contigs", DEFAULT_SHORT_CONTIGS) or DEFAULT_SHORT_CONTIGS
    if short_contigs not in SHORT_CONTIG_MODES:
        raise SystemExit(f"--short-contigs must be one of {SHORT_CONTIG_MODES}")
    short_windows_df, short_report = build_short_contig_windows(
        str(input_path), genome_id, window_len=cfg.data.window_len,
    )
    if short_report["n_short"]:
        # never silent: this was the failure mode (docs/cli-and-packaging.md "Input robustness")
        verb = "scanned at reduced scale" if short_contigs == "scan" else "NOT scanned (--short-contigs skip)"
        print(
            f"[warn] {short_report['n_short']} contig(s) shorter than {cfg.data.window_len} bp "
            f"({short_report['bp_short']} bp): {short_report['n_scannable']} {verb}; "
            f"{short_report['n_too_short']} under {min(SHORT_CONTIG_SCALES)} bp never scanned "
            f"({short_report['bp_too_short']} bp).",
            file=sys.stderr,
        )
    if windows_df.empty and (short_contigs == "skip" or short_windows_df.empty):
        raise RuntimeError("No inference windows were generated from the input FASTA.")

    ambiguity_policy = getattr(args, "ambiguity_policy", DEFAULT_AMBIGUITY_POLICY) or DEFAULT_AMBIGUITY_POLICY
    checkpoints = args.checkpoint if isinstance(args.checkpoint, list) else [args.checkpoint]
    ensemble = getattr(args, "ensemble", DEFAULT_ENSEMBLE)
    if ensemble not in ENSEMBLE_MODES:
        raise SystemExit(f"--ensemble must be one of {ENSEMBLE_MODES}")

    model = build_model(cfg, device)
    score_kw = dict(ambiguity_policy=ambiguity_policy, ensemble=ensemble, device=device,
                    batch_size=args.batch_size, num_workers=args.num_workers)
    bins_df = None
    # seed support is computed for the coverage builder (the shipped one) only
    track_support = len(checkpoints) > 1 and region_builder == "coverage" and num_tiles <= 1
    min_support = getattr(args, "min_seed_support", None)
    min_support = MODE_MIN_SEED_SUPPORT[args.mode] if min_support is None else int(min_support)
    if min_support < 1:
        raise SystemExit("--min-seed-support must be >= 1")
    vote_k = min(min_support, len(checkpoints)) if track_support else 1
    if min_support > 1 and vote_k < min_support:
        # never silent: the mode thresholds assume the vote (MODE_THRESHOLDS comment)
        why = (f"only {len(checkpoints)} checkpoint(s)" if len(checkpoints) < min_support
               else "the vote needs --region-builder coverage and --num-tiles 1")
        print(f"[warn] seed vote of {min_support} requested but {why}: "
              f"{'voting at ' + str(vote_k) if vote_k > 1 else 'no vote applied'}. The '{args.mode}' "
              "thresholds were calibrated with the vote over 3 checkpoints; expect more calls.", file=sys.stderr)
    vote_dropped = 0
    write_seed_scores = bool(getattr(args, "write_seed_scores", False)) and len(checkpoints) > 1
    per_seed_regions: List[List[pd.DataFrame]] = [[] for _ in checkpoints]
    if windows_df.empty:  # every contig is short: nothing for the coarse pass
        windows_df = pd.DataFrame(columns=["genome_id", "contig_id", "start", "end", "fasta_path",
                                           "prob", "pred_label", "region_pred_label"])
        regions_df = pd.DataFrame(columns=REGION_COLUMNS)
    else:
        coarse_per_ckpt: List[np.ndarray] = []
        probs = score_windows(model, checkpoints, windows_df, window_len=cfg.data.window_len,
                              per_checkpoint_out=coarse_per_ckpt, **score_kw)
        windows_df = windows_df.copy()
        windows_df["prob"] = probs
        windows_df["pred_label"] = (windows_df["prob"] >= threshold).astype(int)
        if write_seed_scores:
            for k, p in enumerate(coarse_per_ckpt):
                windows_df[f"prob_ckpt{k + 1}"] = p

    if windows_df.empty:
        pass
    elif num_tiles > 1:
        bin_size = max(1, cfg.data.stride // num_tiles)
        bins_df, regions_df = build_region_table_multi_tile(
            windows_df,
            bin_size=bin_size,
            seed_threshold=seed_threshold,
            extend_threshold=extend_threshold,
            window_len=cfg.data.window_len,
            stride=cfg.data.stride,
        )
    elif region_builder == "coverage":
        # per-checkpoint scores are aligned to this (pre-builder) row order
        coarse_coords = windows_df[["genome_id", "contig_id", "start", "end", "fasta_path"]]
        windows_df, regions_df = build_region_table_coverage(
            windows_df,
            seed_threshold=seed_threshold,
            extend_threshold=extend_threshold,
            min_coverage_frac=coverage_frac,
            max_gap_bp=max_gap_bp,
            bridge_prob=bridge_prob,
        )
        if len(checkpoints) > 1:
            for k, p in enumerate(coarse_per_ckpt):
                _, sr = build_region_table_coverage(
                    coarse_coords.assign(prob=p),
                    seed_threshold=seed_threshold, extend_threshold=extend_threshold,
                    min_coverage_frac=coverage_frac, max_gap_bp=max_gap_bp, bridge_prob=bridge_prob,
                )
                per_seed_regions[k].append(sr)
            if vote_k > 1:
                regions_df, n = apply_seed_vote(regions_df, [r[-1] for r in per_seed_regions], vote_k)
                vote_dropped += n
                if n:
                    windows_df = label_windows(windows_df, regions_df)
    else:
        windows_df, regions_df = build_region_table(
            windows_df,
            seed_threshold=seed_threshold,
            extend_threshold=extend_threshold,
        )

    fine_windows_df = None
    refine_counts = None
    if refine_window_len and not regions_df.empty:
        contig_lengths = dict(iter_fasta_lengths(str(input_path)))
        neighborhoods = plan_neighborhoods(regions_df, contig_lengths, refine_margin)
        fine_windows_df = plan_fine_windows(
            neighborhoods, contig_lengths, genome_id, str(input_path),
            window_len=refine_window_len, stride=refine_stride,
        )
        if not fine_windows_df.empty:
            fine_windows_df = fine_windows_df.copy()
            fine_windows_df["prob"] = score_windows(
                model, checkpoints, fine_windows_df, window_len=refine_window_len, **score_kw
            )
            # one threshold for seed and extend at the fine scale; no merging
            # here -- merging belongs to the coarse track (merge coarse, split fine)
            fine_windows_df, fine_regions_df = build_region_table_coverage(
                fine_windows_df,
                seed_threshold=refine_threshold,
                extend_threshold=refine_threshold,
                min_coverage_frac=coverage_frac,
                max_gap_bp=0,
            )
        else:
            fine_regions_df = pd.DataFrame(columns=["contig_id", "start", "end"])
        regions_df = refine_regions(
            regions_df,
            fine_windows_df,
            fine_regions_df,
            edge_alpha=edge_alpha,
            split_min_gap_bp=split_min_gap_bp,
            split_max_prob=split_max_prob,
            kind_suffix=f"_refine{refine_window_len}_a{edge_alpha:g}"
            + (f"_ramp{alpha_ramp[0]}-{alpha_ramp[1]}" if alpha_ramp else ""),
            alpha_ramp_bp=alpha_ramp,
        )
        refine_counts = {k: int(v) for k, v in regions_df["refine_status"].value_counts().items()}
        # region_pred_label must describe the refined calls, not the coarse ones
        windows_df = label_windows(windows_df, regions_df)

    short_meta = {"mode": short_contigs, **short_report}
    if short_contigs == "scan" and not short_windows_df.empty:
        st = getattr(args, "short_contig_threshold", None)
        calibrated = st is None

        def short_thresholds(scale: int) -> Tuple[float, float]:
            t = float(st) if st is not None else SHORT_CONTIG_THRESHOLDS[args.mode][scale]
            return t, t

        if not calibrated:
            print(f"[warn] --short-contig-threshold {float(st):g} overrides the calibrated "
                  "short-contig thresholds.", file=sys.stderr)
        short_regions, short_windows, used = [], [], {}
        for scale, sw in short_windows_df.groupby("scale"):
            short_seed, short_extend = short_thresholds(int(scale))
            used[str(int(scale))] = {"seed": short_seed, "extend": short_extend}
            sw = sw.drop(columns="scale").reset_index(drop=True)
            short_per_ckpt: List[np.ndarray] = []
            sw["prob"] = score_windows(model, checkpoints, sw, window_len=int(scale),
                                       per_checkpoint_out=short_per_ckpt, **score_kw)
            sw["pred_label"] = (sw["prob"] >= short_seed).astype(int)
            if write_seed_scores:
                for k, p in enumerate(short_per_ckpt):
                    sw[f"prob_ckpt{k + 1}"] = p
            if track_support:
                for k, p in enumerate(short_per_ckpt):
                    _, ssr = build_region_table_coverage(
                        sw[["genome_id", "contig_id", "start", "end", "fasta_path"]].assign(prob=p),
                        seed_threshold=short_seed, extend_threshold=short_extend,
                        min_coverage_frac=coverage_frac, max_gap_bp=0,
                    )
                    per_seed_regions[k].append(ssr)
            # per-contig windows share one length, so the coverage builder's
            # per-contig window_len/stride inference holds; no merging on short contigs
            sw, sr = build_region_table_coverage(
                sw, seed_threshold=short_seed, extend_threshold=short_extend,
                min_coverage_frac=coverage_frac, max_gap_bp=0,
            )
            if vote_k > 1:
                sr, n = apply_seed_vote(sr, [r[-1] for r in per_seed_regions], vote_k)
                vote_dropped += n
                if n:
                    sw = label_windows(sw, sr)
            sr["boundary_kind"] = sr["boundary_kind"].astype(str) + f"_short{int(scale)}"
            short_regions.append(sr)
            short_windows.append(sw)
        short_regions_df = pd.concat(short_regions, ignore_index=True)
        if "refine_status" in regions_df.columns:
            short_regions_df = short_regions_df.assign(
                refine_status="short_contig", coarse_start=short_regions_df["start"],
                coarse_end=short_regions_df["end"], n_pieces=1)
        regions_df = pd.concat([regions_df, short_regions_df], ignore_index=True) \
            if not regions_df.empty else short_regions_df
        windows_df = pd.concat([windows_df, *short_windows], ignore_index=True) \
            if not windows_df.empty else pd.concat(short_windows, ignore_index=True)
        short_meta.update({"thresholds": used,
                           "scales": list(SHORT_CONTIG_SCALES), "stride": SHORT_CONTIG_STRIDE,
                           "n_regions": int(len(short_regions_df)), "calibrated": calibrated})

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    windows_out = outdir / "windows.tsv"
    regions_out = outdir / "regions.tsv"
    metadata_out = outdir / "run_metadata.json"

    support_meta = None
    if track_support:
        seed_regs = [pd.concat([x for x in r if len(x)], ignore_index=True) if any(len(x) for x in r)
                     else pd.DataFrame(columns=["contig_id", "start", "end"]) for r in per_seed_regions]
        regions_df = regions_df.copy()
        regions_df["seed_support"] = seed_support(regions_df, seed_regs) if len(regions_df) else []
        support_meta = {"n_checkpoints": len(checkpoints),
                        "min_seed_support": vote_k, "requested_min_seed_support": min_support,
                        "calls_dropped_by_vote": vote_dropped,
                        "calls_by_support": {str(k): int(v) for k, v in
                                             regions_df["seed_support"].value_counts().sort_index().items()}}

    windows_df.to_csv(windows_out, sep="\t", index=False)
    regions_df.to_csv(regions_out, sep="\t", index=False)
    write_metadata(
        cfg=cfg,
        checkpoint_path=checkpoints if len(checkpoints) > 1 else checkpoints[0],
        ensemble_mode=(ensemble if len(checkpoints) > 1 else "single"),
        fasta_path=str(input_path),
        mode=args.mode,
        threshold=threshold,
        seed_threshold=seed_threshold,
        extend_threshold=extend_threshold,
        device=str(device),
        windows_df=windows_df,
        regions_df=regions_df,
        out_path=metadata_out,
        num_tiles=num_tiles,
        ambiguity_policy=ambiguity_policy,
        grid_offset=max(0, int(getattr(args, "grid_offset", 0) or 0)),
        region_builder=region_builder,
        coverage_frac=coverage_frac,
        max_gap_bp=max_gap_bp,
        bridge_prob=bridge_prob,
        trained_window_len=trained_window_len,
        trained_stride=trained_stride,
        short_contigs=short_meta,
        seed_support=support_meta,
        extra={"input": input_summary,
               "checkpoint_sha256": {str(c): sha256(str(c)) for c in checkpoints},
               "outputs": ["windows.tsv", "regions.tsv", "regions.bed", "regions.gff3", "run_metadata.json"]
               + (["regions.fna"] if getattr(args, "write_region_fasta", False) else [])},
        refine=(None if not refine_window_len else {
            "window_len": refine_window_len, "stride": refine_stride,
            "threshold": refine_threshold, "margin_bp": refine_margin,
            "edge_alpha": edge_alpha, "edge_alpha_ramp_bp": list(alpha_ramp) if alpha_ramp else None,
            "split_min_gap_bp": split_min_gap_bp,
            "split_max_prob": split_max_prob, "n_fine_windows": int(0 if fine_windows_df is None else len(fine_windows_df)),
            "status_counts": refine_counts or {},
        }),
    )

    print(f"[done] wrote {windows_out}", file=sys.stderr)
    if fine_windows_df is not None:
        fine_out = outdir / "fine_windows.tsv"
        fine_windows_df.to_csv(fine_out, sep="\t", index=False)
        print(f"[done] wrote {fine_out}", file=sys.stderr)
    print(f"[done] wrote {regions_out}", file=sys.stderr)
    print(f"[done] wrote {metadata_out}", file=sys.stderr)

    if bins_df is not None:
        bins_out = outdir / "bins.tsv"
        bins_df.to_csv(bins_out, sep="\t", index=False)
        print(f"[done] wrote {bins_out}", file=sys.stderr)

    write_bed(regions_df, outdir / "regions.bed")
    write_gff3(regions_df, outdir / "regions.gff3", args.mode, package_version())
    print(f"[done] wrote {outdir / 'regions.bed'} and {outdir / 'regions.gff3'}", file=sys.stderr)
    if getattr(args, "write_region_fasta", False):
        write_region_fasta(regions_df, str(input_path), outdir / "regions.fna")
        print(f"[done] wrote {outdir / 'regions.fna'}", file=sys.stderr)
    print_run_summary(args.mode, input_summary, regions_df, short_meta, support_meta, ambiguity_policy,
                      len(checkpoints))


def print_run_summary(mode, inp, regions_df, short_meta, support_meta, ambiguity_policy, n_ckpt) -> None:
    """One human-readable block on stderr: what went in, what came out, what was skipped."""
    lines = [f"prophound {package_version()} | mode {mode} | {n_ckpt} checkpoint(s)",
             f"input: {inp['n_contigs']} contig(s), {inp['total_bp']:,} bp; "
             f"ambiguity codes {inp['ambiguous_bases']:,} ({inp['ambiguous_ppm']:g} ppm), policy {ambiguity_policy}"]
    if short_meta and short_meta.get("n_short"):
        lines.append(f"short contigs: {short_meta['n_short']} under the window length ({short_meta['bp_short']:,} bp); "
                     f"{short_meta['n_too_short']} never scanned ({short_meta['bp_too_short']:,} bp)")
    n = len(regions_df)
    bp = int((regions_df.end - regions_df.start).sum()) if n else 0
    lines.append(f"calls: {n} ({bp:,} bp, {100 * bp / max(1, inp['total_bp']):.1f}% of input)")
    if support_meta:
        lines.append(f"seed vote: min support {support_meta['min_seed_support']}, "
                     f"{support_meta['calls_dropped_by_vote']} call(s) dropped; kept by support "
                     f"{support_meta['calls_by_support']}")
    lines.append("calls are candidates for follow-up, not annotations")
    print("\n".join("[summary] " + l for l in lines), file=sys.stderr)


def main() -> None:
    run_predict(build_arg_parser().parse_args())


if __name__ == "__main__":
    main()
