from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd
import torch
from pyfaidx import Fasta
from torch.utils.data import DataLoader, Dataset
from tqdm import tqdm

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

MODE_THRESHOLDS = {
    "exploratory": 0.95,
    "default": 0.995,
    "conservative": 0.999,
}


class GenomeWindowDataset(Dataset):
    """Window dataset for direct inference on a FASTA file."""

    def __init__(self, windows_df: pd.DataFrame, window_len: int):
        self.windows_df = windows_df.reset_index(drop=True)
        self.window_len = window_len
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
            sequence = sequence + ("N" * (self.window_len - len(sequence)))
        elif len(sequence) > self.window_len:
            sequence = sequence[: self.window_len]
        return onehot_encode(sequence), torch.tensor(0.0, dtype=torch.float32)

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
        logits = model(batch_x).view(-1)
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


def build_multi_tile_windows(
    fasta_path: str,
    genome_id: str,
    data_config: DataConfig,
    num_tiles: int,
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
        return build_windows(fasta_path, genome_id, data_config)

    tile_stride = max(1, data_config.stride // num_tiles)
    tiles = []
    for tile_index in range(num_tiles):
        tile_df = build_windows(
            fasta_path, genome_id, data_config, offset=tile_index * tile_stride
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
      column == 0), i.e. exactly the same single regular grid and mean-based
      bin aggregation `build_region_table` would use. This carries over the
      frozen v1 detection behavior with no new false positives and no new
      dilution -- seeding is not where the finer tiling should help anyway.
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
        for row in tile0_df.itertuples(index=False):
            cb_start = int(row.start // stride)
            cb_end = min(coarse_n_bins, int(math.ceil(row.end / stride)))
            coarse_sum[cb_start:cb_end] += float(row.prob)
            coarse_count[cb_start:cb_end] += 1
        coarse_covered = coarse_count > 0
        coarse_scores = np.zeros(coarse_n_bins, dtype=np.float64)
        coarse_scores[coarse_covered] = coarse_sum[coarse_covered] / coarse_count[coarse_covered]

        coarse_seed_mask = coarse_covered & (coarse_scores >= seed_threshold)
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
    checkpoint_path: str,
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
) -> None:
    metadata = {
        "run_name": cfg.run_name,
        "checkpoint": checkpoint_path,
        "input_fasta": fasta_path,
        "mode": mode,
        "threshold": threshold,
        "seed_threshold": seed_threshold,
        "extend_threshold": extend_threshold,
        "device": device,
        "window_len": cfg.data.window_len,
        "stride": cfg.data.stride,
        "num_tiles": num_tiles,
        "n_windows": int(len(windows_df)),
        "n_regions": int(len(regions_df)),
    }
    if num_tiles > 1:
        metadata["region_builder"] = "multi_tile_hysteresis"
        metadata["boundary_refinement"] = "multi_tile_ramp_midpoint_crossing"
        metadata["bin_size"] = max(1, cfg.data.stride // num_tiles)
    else:
        metadata["region_builder"] = "strict_hysteresis_stride_support"
        metadata["boundary_refinement"] = "interpolated_from_stride_support_track"
    out_path.write_text(json.dumps(metadata, indent=2) + "\n")


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run frozen MGE-Explorer inference on a genome FASTA."
    )
    parser.add_argument("--config", required=True)
    parser.add_argument("--checkpoint", required=True)
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
    parser.add_argument("--device", default="auto")
    parser.add_argument("--batch-size", type=int, default=128)
    parser.add_argument("--num-workers", type=int, default=0)
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
    return parser


def run_predict(args: argparse.Namespace) -> None:
    """Run inference end to end from a parsed argument namespace.

    Shared by this module's own CLI and by `mge_detect`'s `predict` subcommand,
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

    num_tiles = max(1, getattr(args, "num_tiles", 1))

    input_path = Path(args.input_fasta)
    genome_id = input_path.stem
    windows_df = build_multi_tile_windows(
        str(input_path), genome_id=genome_id, data_config=cfg.data, num_tiles=num_tiles
    )
    if windows_df.empty:
        raise RuntimeError("No inference windows were generated from the input FASTA.")

    dataset = GenomeWindowDataset(windows_df, window_len=cfg.data.window_len)
    loader = DataLoader(
        dataset,
        batch_size=args.batch_size,
        shuffle=False,
        num_workers=args.num_workers,
        collate_fn=collate_batch,
        pin_memory=(device.type == "cuda"),
    )

    model = build_model(cfg, device)
    checkpoint = torch.load(args.checkpoint, map_location=device, weights_only=False)
    model.load_state_dict(checkpoint["model"])
    probs = predict_probs(model, loader, device)
    dataset.close()

    windows_df = windows_df.copy()
    windows_df["prob"] = probs
    windows_df["pred_label"] = (windows_df["prob"] >= threshold).astype(int)

    bins_df = None
    if num_tiles > 1:
        bin_size = max(1, cfg.data.stride // num_tiles)
        bins_df, regions_df = build_region_table_multi_tile(
            windows_df,
            bin_size=bin_size,
            seed_threshold=seed_threshold,
            extend_threshold=extend_threshold,
            window_len=cfg.data.window_len,
            stride=cfg.data.stride,
        )
    else:
        windows_df, regions_df = build_region_table(
            windows_df,
            seed_threshold=seed_threshold,
            extend_threshold=extend_threshold,
        )

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    windows_out = outdir / "windows.tsv"
    regions_out = outdir / "regions.tsv"
    metadata_out = outdir / "run_metadata.json"

    windows_df.to_csv(windows_out, sep="\t", index=False)
    regions_df.to_csv(regions_out, sep="\t", index=False)
    write_metadata(
        cfg=cfg,
        checkpoint_path=args.checkpoint,
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
    )

    print(f"[done] wrote {windows_out}")
    print(f"[done] wrote {regions_out}")
    print(f"[done] wrote {metadata_out}")

    if bins_df is not None:
        bins_out = outdir / "bins.tsv"
        bins_df.to_csv(bins_out, sep="\t", index=False)
        print(f"[done] wrote {bins_out}")


def main() -> None:
    run_predict(build_arg_parser().parse_args())


if __name__ == "__main__":
    main()
