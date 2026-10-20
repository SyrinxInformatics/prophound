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
from .options import (  # noqa: F401 -- MODE_THRESHOLDS, DEFAULT_EDGE_SNAP_FRAC, MAX_CALL_N_FRAC are imported from here by scripts
    AMBIGUITY_POLICIES,
    DEFAULT_AMBIGUITY_POLICY,
    DEFAULT_BRIDGE_PROB,
    DEFAULT_COVERAGE_FRAC,
    DEFAULT_EDGE_SNAP_FRAC,
    DEFAULT_MAX_GAP_BP,
    DEFAULT_REFINE_THRESHOLD,
    ENSEMBLE_MODES,
    MAX_CALL_N_FRAC,
    MODE_MIN_SEED_SUPPORT,
    MODE_THRESHOLDS,
    SHORT_CONTIG_MODES,
    SHORT_CONTIG_SCALES,
    SHORT_CONTIG_STRIDE,
    SHORT_CONTIG_THRESHOLDS,
    add_predict_arguments,
)
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
    ambiguity_policy: str = DEFAULT_AMBIGUITY_POLICY,
    grid_offset: int = 0,
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
    metadata["region_builder"] = "coverage_fraction"
    metadata["boundary_refinement"] = "none_bin_edges"
    metadata["coverage_frac"] = coverage_frac
    if max_gap_bp > 0:
        metadata["fragment_merge"] = {"max_gap_bp": max_gap_bp, "bridge_prob": bridge_prob}
    if refine:
        metadata["refine"] = refine
    out_path.write_text(json.dumps(metadata, indent=2) + "\n")


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


def call_n_fraction(regions_df: pd.DataFrame, fasta_path: str) -> np.ndarray:
    """Fraction of each call's bases that are N."""
    if regions_df.empty:
        return np.zeros(0)
    fa = Fasta(fasta_path, as_raw=True, sequence_always_upper=True)
    try:
        return np.array([fa[r.contig_id][int(r.start):int(r.end)].count("N") / max(1, int(r.end) - int(r.start))
                         for r in regions_df.itertuples(index=False)])
    finally:
        fa.close()


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
    window_len_arg = args.window_len
    stride_arg = args.stride
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

    coverage_frac = float(args.coverage_frac)
    max_gap_bp = int(args.max_gap_bp)
    bridge_prob = float(args.bridge_prob)
    refine_window_len = int(args.refine_window_len or 0)
    refine_stride = int(args.refine_stride)
    edge_alpha = float(args.edge_alpha)
    alpha_ramp = args.edge_alpha_ramp
    alpha_ramp = tuple(int(x) for x in alpha_ramp) if alpha_ramp else None
    if alpha_ramp and not 0 < alpha_ramp[0] < alpha_ramp[1]:
        raise SystemExit("--edge-alpha-ramp needs 0 < LO_BP < HI_BP")
    split_min_gap_bp = int(args.split_min_gap_bp)
    split_max_prob = float(args.split_max_prob)
    refine_threshold = args.refine_threshold
    refine_threshold = DEFAULT_REFINE_THRESHOLD if refine_threshold is None else float(refine_threshold)
    refine_margin = args.refine_margin
    refine_margin = cfg.data.window_len if refine_margin is None else int(refine_margin)
    if refine_window_len:
        if refine_window_len >= cfg.data.window_len:
            raise SystemExit("--refine-window-len must be smaller than --window-len")
        if refine_stride <= 0 or not 0.0 <= edge_alpha <= 1.0 or refine_margin < 0:
            raise SystemExit("need --refine-stride > 0, 0 <= --edge-alpha <= 1, --refine-margin >= 0")

    edge_snap_frac = float(args.edge_snap_frac)
    if not 0.0 <= edge_snap_frac <= 1.0:
        raise SystemExit("--edge-snap-frac must be between 0 and 1")

    input_path = Path(args.input_fasta)
    input_summary = validate_fasta(str(input_path))   # fails early with a clear message
    genome_id = input_path.stem
    windows_df = build_windows(
        str(input_path),
        genome_id=genome_id,
        data_config=cfg.data,
        offset=max(0, int(args.grid_offset or 0)),
    )
    short_contigs = args.short_contigs
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

    ambiguity_policy = args.ambiguity_policy
    checkpoints = args.checkpoint if isinstance(args.checkpoint, list) else [args.checkpoint]
    ensemble = args.ensemble
    if ensemble not in ENSEMBLE_MODES:
        raise SystemExit(f"--ensemble must be one of {ENSEMBLE_MODES}")

    model = build_model(cfg, device)
    score_kw = dict(ambiguity_policy=ambiguity_policy, ensemble=ensemble, device=device,
                    batch_size=args.batch_size, num_workers=args.num_workers)
    track_support = len(checkpoints) > 1
    min_support = args.min_seed_support
    min_support = MODE_MIN_SEED_SUPPORT[args.mode] if min_support is None else int(min_support)
    if min_support < 1:
        raise SystemExit("--min-seed-support must be >= 1")
    vote_k = min(min_support, len(checkpoints)) if track_support else 1
    if min_support > 1 and vote_k < min_support:
        # never silent: the mode thresholds assume the vote (MODE_THRESHOLDS comment)
        print(f"[warn] seed vote of {min_support} requested but only {len(checkpoints)} checkpoint(s): "
              f"{'voting at ' + str(vote_k) if vote_k > 1 else 'no vote applied'}. The '{args.mode}' "
              "thresholds were calibrated with the vote over 3 checkpoints; expect more calls.", file=sys.stderr)
    vote_dropped = 0
    write_seed_scores = bool(args.write_seed_scores) and len(checkpoints) > 1
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

    if not windows_df.empty:
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
        st = args.short_contig_threshold
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

    # before the N filter, so n_frac describes the reported coordinates
    snap_bp = int(edge_snap_frac * cfg.data.window_len)
    regions_df, n_snapped = snap_contig_edges(regions_df, dict(iter_fasta_lengths(str(input_path))), snap_bp)
    if n_snapped:
        windows_df = label_windows(windows_df, regions_df)
    edge_meta = {"frac": edge_snap_frac, "snap_bp": snap_bp, "calls_moved": n_snapped,
                 "calls_on_contig_edge": int((regions_df["contig_edge"] != "none").sum()) if len(regions_df) else 0}

    # calls made mostly of N (imputed scaffold gaps) are not candidates
    n_frac = call_n_fraction(regions_df, str(input_path))
    n_dropped = int((n_frac >= MAX_CALL_N_FRAC).sum())
    regions_df = regions_df.assign(n_frac=np.round(n_frac, 4))
    if n_dropped:
        regions_df = regions_df[regions_df.n_frac < MAX_CALL_N_FRAC].reset_index(drop=True)
        windows_df = label_windows(windows_df, regions_df)
        print(f"[warn] {n_dropped} call(s) dropped: at least {MAX_CALL_N_FRAC:.0%} N (imputed gap sequence).",
              file=sys.stderr)
    gap_filter_meta = {"max_call_n_frac": MAX_CALL_N_FRAC, "calls_dropped": n_dropped,
                       "calls_with_any_n": int((regions_df.n_frac > 0).sum())}

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
        ambiguity_policy=ambiguity_policy,
        grid_offset=max(0, int(args.grid_offset or 0)),
        coverage_frac=coverage_frac,
        max_gap_bp=max_gap_bp,
        bridge_prob=bridge_prob,
        trained_window_len=trained_window_len,
        trained_stride=trained_stride,
        short_contigs=short_meta,
        seed_support=support_meta,
        extra={"input": input_summary, "n_gap_filter": gap_filter_meta, "edge_snap": edge_meta,
               "checkpoint_sha256": {str(c): sha256(str(c)) for c in checkpoints},
               "outputs": ["windows.tsv", "regions.tsv", "regions.bed", "regions.gff3", "run_metadata.json"]
               + (["regions.fna"] if args.write_region_fasta else [])},
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

    write_bed(regions_df, outdir / "regions.bed")
    write_gff3(regions_df, outdir / "regions.gff3", args.mode, package_version())
    print(f"[done] wrote {outdir / 'regions.bed'} and {outdir / 'regions.gff3'}", file=sys.stderr)
    if args.write_region_fasta:
        write_region_fasta(regions_df, str(input_path), outdir / "regions.fna")
        print(f"[done] wrote {outdir / 'regions.fna'}", file=sys.stderr)
    print_run_summary(args.mode, input_summary, regions_df, short_meta, support_meta, ambiguity_policy,
                      len(checkpoints), gap_filter_meta)


def print_run_summary(mode, inp, regions_df, short_meta, support_meta, ambiguity_policy, n_ckpt,
                      gap_meta=None) -> None:
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
    if gap_meta and (gap_meta["calls_dropped"] or gap_meta["calls_with_any_n"]):
        lines.append(f"N gaps: {gap_meta['calls_dropped']} call(s) dropped as >= {gap_meta['max_call_n_frac']:.0%} N; "
                     f"{gap_meta['calls_with_any_n']} kept call(s) contain some N (see n_frac)")
    lines.append("calls are candidates for follow-up, not annotations")
    print("\n".join("[summary] " + l for l in lines), file=sys.stderr)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run prophage inference on a genome FASTA (same as `prophound predict`).")
    add_predict_arguments(parser)
    return parser


def main(argv: List[str] | None = None) -> None:
    # `python -m src.predict_genome` launchers spell the input `--input_fasta`
    argv = sys.argv[1:] if argv is None else list(argv)
    argv = ["--input" + a[len("--input_fasta"):] if a == "--input_fasta" or a.startswith("--input_fasta=") else a
            for a in argv]
    run_predict(build_arg_parser().parse_args(argv))


if __name__ == "__main__":
    main()
