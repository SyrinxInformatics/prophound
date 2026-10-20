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
"""Defaults, mode tables and the `predict` argument parser.

The one place these are defined. Kept free of torch so `prophound --help` stays fast.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

# Per mode: the seed threshold, and how many checkpoints must make an overlapping call
# of their own for an ensemble call to be kept (the seed vote). Exploratory keeps every
# call. Derivation and history: docs/cli-and-packaging.md "Seed vote per mode".
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
# The model's N channel is effectively untrained, so 'n-channel' can zero a window's
# score; 'impute' keeps the input in distribution (docs/scientific-validity.md "1a").
DEFAULT_AMBIGUITY_POLICY = "impute"

# Contigs shorter than --window-len are tiled with the largest scale that fits ('scan')
# or left unscanned ('skip'); either way they are reported.
SHORT_CONTIG_MODES = ("skip", "scan")
DEFAULT_SHORT_CONTIGS = "scan"
SHORT_CONTIG_SCALES = (5000, 3000)   # below ~3 kb a window cannot judge
SHORT_CONTIG_STRIDE = 500
# Seed and extend threshold per mode and scale: the lowest whose false calls per Mb on
# shredded genomes, with the mode's vote, do not exceed the intact 10 kb pipeline's
# (docs/cli-and-packaging.md).
SHORT_CONTIG_THRESHOLDS = {
    "exploratory": {5000: 0.95, 3000: 0.95},
    "default": {5000: 0.95, 3000: 0.95},
    "conservative": {5000: 0.99, 3000: 0.99},
}

ENSEMBLE_MODES = ("max", "mean")
DEFAULT_ENSEMBLE = "max"
# The coverage builder and fragment merging at the measured operating point
# (docs/model-quality-baseline.md "Operating point").
DEFAULT_COVERAGE_FRAC = 0.1
DEFAULT_MAX_GAP_BP = 50000      # 0 turns merging off
DEFAULT_BRIDGE_PROB = 0.1
# The fine refinement pass is off by default (0); the other settings are the opt-in
# "tight boundaries" setting (docs/downstream-postprocessing.md).
DEFAULT_REFINE_WINDOW_LEN = 0
DEFAULT_REFINE_STRIDE = 250
DEFAULT_REFINE_THRESHOLD = 0.95
DEFAULT_EDGE_ALPHA = 1.0        # 1 keeps the coarse edges
DEFAULT_SPLIT_MIN_GAP_BP = 4000
# Below DEFAULT_BRIDGE_PROB, so splitting cannot undo a merge.
DEFAULT_SPLIT_MAX_PROB = 0.05
# A call at least this fraction N is dropped: `impute` turns a scaffold gap into random
# sequence, which can be called (docs/v1-scope.md A.3).
MAX_CALL_N_FRAC = 0.5
# A call edge within this fraction of --window-len of a contig end moves onto the end:
# an element cut by a breakpoint is only partly scored there
# (docs/claims-and-evidence.md "G5b"). 0 turns snapping off.
DEFAULT_EDGE_SNAP_FRAC = 0.5

DEVICES = ("auto", "cpu", "cuda", "rocm", "mps", "xpu", "intel")
DEFAULT_BATCH_SIZE = 128


@dataclass(frozen=True)
class PredictOptions:
    """A `predict` run's settings, resolved against the checkpoint and validated once."""

    input_fasta: str
    outdir: str
    mode: str
    config: str
    checkpoints: tuple[str, ...]
    device: str
    threshold: float
    seed_threshold: float
    extend_threshold: float
    requested_min_seed_support: int
    ensemble: str
    write_region_fasta: bool
    write_seed_scores: bool        # effective: needs several checkpoints
    ambiguity_policy: str
    grid_offset: int
    batch_size: int
    num_workers: int
    coverage_frac: float
    max_gap_bp: int
    bridge_prob: float
    window_len: int
    stride: int
    trained_window_len: int
    trained_stride: int
    edge_snap_frac: float
    short_contigs: str
    short_contig_threshold: float | None
    refine_window_len: int          # 0 = refinement off
    refine_stride: int
    refine_threshold: float
    refine_margin: int
    edge_alpha: float
    edge_alpha_ramp: tuple[int, int] | None
    split_min_gap_bp: int
    split_max_prob: float

    @property
    def n_checkpoints(self) -> int:
        return len(self.checkpoints)

    @property
    def min_seed_support(self) -> int:
        """The vote actually applied: capped at the number of checkpoints, 1 with one."""
        if self.n_checkpoints < 2:
            return 1
        return min(self.requested_min_seed_support, self.n_checkpoints)

    @property
    def snap_bp(self) -> int:
        return int(self.edge_snap_frac * self.window_len)

    def short_threshold(self, scale: int) -> float:
        """Seed and extend threshold for a short-contig scale (one value serves both)."""
        if self.short_contig_threshold is not None:
            return float(self.short_contig_threshold)
        return SHORT_CONTIG_THRESHOLDS[self.mode][scale]

    @classmethod
    def from_args(cls, args: argparse.Namespace, trained_window_len: int, trained_stride: int) -> PredictOptions:
        """Resolve mode defaults and window overrides; exit with a message on an invalid combination."""
        seed_default = args.threshold if args.threshold is not None else MODE_THRESHOLDS[args.mode]
        seed = args.seed_threshold if args.seed_threshold is not None else seed_default
        extend = args.extend_threshold if args.extend_threshold is not None else min(seed, 0.99)
        window_len = int(args.window_len) if args.window_len else trained_window_len
        stride = int(args.stride) if args.stride else trained_stride
        if window_len <= 0 or stride <= 0:
            raise SystemExit("--window-len and --stride must be positive")
        ramp = (int(args.edge_alpha_ramp[0]), int(args.edge_alpha_ramp[1])) if args.edge_alpha_ramp else None
        if ramp and not 0 < ramp[0] < ramp[1]:
            raise SystemExit("--edge-alpha-ramp needs 0 < LO_BP < HI_BP")
        refine_threshold = DEFAULT_REFINE_THRESHOLD if args.refine_threshold is None else float(args.refine_threshold)
        refine_margin = window_len if args.refine_margin is None else int(args.refine_margin)
        refine_window_len = int(args.refine_window_len or 0)
        if refine_window_len:
            if refine_window_len >= window_len:
                raise SystemExit("--refine-window-len must be smaller than --window-len")
            if args.refine_stride <= 0 or not 0.0 <= args.edge_alpha <= 1.0 or refine_margin < 0:
                raise SystemExit("need --refine-stride > 0, 0 <= --edge-alpha <= 1, --refine-margin >= 0")
        if not 0.0 <= args.edge_snap_frac <= 1.0:
            raise SystemExit("--edge-snap-frac must be between 0 and 1")
        min_support = MODE_MIN_SEED_SUPPORT[args.mode] if args.min_seed_support is None else int(args.min_seed_support)
        if min_support < 1:
            raise SystemExit("--min-seed-support must be >= 1")
        checkpoints = tuple(args.checkpoint) if isinstance(args.checkpoint, list) else (args.checkpoint,)
        return cls(
            input_fasta=str(Path(args.input_fasta)), outdir=str(args.outdir), mode=args.mode, config=args.config,
            checkpoints=checkpoints, device=args.device,
            threshold=float(seed_default), seed_threshold=float(seed), extend_threshold=float(extend),
            requested_min_seed_support=min_support, ensemble=args.ensemble,
            write_region_fasta=bool(args.write_region_fasta),
            write_seed_scores=bool(args.write_seed_scores) and len(checkpoints) > 1,
            ambiguity_policy=args.ambiguity_policy, grid_offset=max(0, int(args.grid_offset or 0)),
            batch_size=args.batch_size, num_workers=args.num_workers,
            coverage_frac=float(args.coverage_frac), max_gap_bp=int(args.max_gap_bp),
            bridge_prob=float(args.bridge_prob), window_len=window_len, stride=stride,
            trained_window_len=trained_window_len, trained_stride=trained_stride,
            edge_snap_frac=float(args.edge_snap_frac), short_contigs=args.short_contigs,
            short_contig_threshold=args.short_contig_threshold, refine_window_len=refine_window_len,
            refine_stride=int(args.refine_stride), refine_threshold=refine_threshold, refine_margin=refine_margin,
            edge_alpha=float(args.edge_alpha), edge_alpha_ramp=ramp,
            split_min_gap_bp=int(args.split_min_gap_bp), split_max_prob=float(args.split_max_prob),
        )


def add_predict_arguments(parser: argparse.ArgumentParser) -> None:
    """Add every `predict` option to `parser`, in basic / advanced / experimental --help sections."""
    basic = parser.add_argument_group("basic options")
    advanced = parser.add_argument_group("advanced options (thresholds, ensemble, input handling, outputs)")
    experimental = parser.add_argument_group("experimental options (not validated for default use; see docs/)")

    basic.add_argument("--input", required=True, dest="input_fasta", help="Nucleotide FASTA to scan.")
    basic.add_argument("--outdir", required=True, help="Directory to write outputs into.")
    basic.add_argument(
        "--mode",
        choices=sorted(MODE_THRESHOLDS),
        default="default",
        help=(
            "Operating mode: seed threshold plus seed vote. exploratory 0.90, every "
            "call kept (~36 calls/genome); default 0.95, all 3 checkpoints agree "
            "(~14); conservative 0.995, all 3 agree (~8)."
        ),
    )
    basic.add_argument(
        "--config",
        required=True,
        help="Training config YAML for the checkpoint (defines model/window architecture).",
    )
    basic.add_argument(
        "--checkpoint",
        required=True,
        nargs="+",
        help=(
            "Trained checkpoint(s). Passing more than one runs a seed ensemble "
            "over them (see --ensemble) -- recommended, since it removes the "
            "dependence of recall on which replicate was packaged."
        ),
    )
    advanced.add_argument("--threshold", type=float, default=None, help="Override the mode's score threshold.")
    advanced.add_argument("--seed-threshold", type=float, default=None)
    advanced.add_argument("--extend-threshold", type=float, default=None)
    advanced.add_argument(
        "--min-seed-support",
        type=int,
        default=None,
        help=(
            "Keep an ensemble call only if at least this many checkpoints make an "
            "overlapping call on their own (default: the mode's; exploratory 1, "
            "default 3, conservative 3). Needs several --checkpoint files; 1 keeps "
            "every call."
        ),
    )
    basic.add_argument(
        "--write-region-fasta",
        action="store_true",
        help="Also write the called regions' sequences to regions.fna.",
    )
    advanced.add_argument(
        "--write-seed-scores",
        action="store_true",
        help=(
            "With several --checkpoint files, also write each checkpoint's score to "
            "windows.tsv (prob_ckpt1, ...) for offline re-analysis."
        ),
    )
    basic.add_argument("--device", default="auto", choices=DEVICES)
    advanced.add_argument(
        "--ambiguity-policy",
        choices=AMBIGUITY_POLICIES,
        default=DEFAULT_AMBIGUITY_POLICY,
        help=(
            "How IUPAC ambiguity codes (N, R, Y, ...) are encoded. 'impute' (default "
            "since 2026-09-24) replaces each with a deterministic concrete base and is "
            "a no-op on ambiguity-free sequence. 'n-channel' is the training-time "
            "encoding and reproduces numbers published before then, but two contiguous "
            "ambiguity codes can take a confident window from 1.000 to 0.000. 'uniform' "
            "spreads each ambiguous base as 0.25 across A/C/G/T."
        ),
    )
    experimental.add_argument(
        "--grid-offset",
        type=int,
        default=0,
        help=(
            "Shift the inference window grid by this many bp. Default 0. For measuring "
            "phase sensitivity: marginal calls can appear or vanish when the grid moves "
            "(see docs/scientific-validity.md 'Input-robustness standing check')."
        ),
    )
    advanced.add_argument("--batch-size", type=int, default=DEFAULT_BATCH_SIZE)
    advanced.add_argument("--num-workers", type=int, default=0)
    advanced.add_argument(
        "--coverage-frac",
        type=float,
        default=DEFAULT_COVERAGE_FRAC,
        help=(
            "Minimum fraction of a bin's covering windows that must score at least the "
            "extend threshold for the bin to count as inside a call (0 -> union of high "
            "windows, 1 -> intersection). Default %(default)s."
        ),
    )
    advanced.add_argument(
        "--max-gap-bp",
        type=int,
        default=DEFAULT_MAX_GAP_BP,
        help=(
            "Merge two calls separated by at "
            "most this many bp into one, repairing elements fragmented by an "
            "internal dip in the score track (worst for 60 kb+ elements, where "
            "it lifts median IoU from ~0.38 to ~0.85). Default %(default)s; 0 "
            "disables merging. Always pair with --bridge-prob -- an ungated "
            "merge bridges genuinely separate neighbouring elements. Above "
            "~50000 the gate binds before the distance budget does. See "
            "docs/gap-filling.md."
        ),
    )
    advanced.add_argument(
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
    experimental.add_argument(
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
    advanced.add_argument(
        "--edge-snap-frac",
        type=float,
        default=DEFAULT_EDGE_SNAP_FRAC,
        help=(
            "Move a call edge within this fraction of --window-len of a contig "
            "end onto that end (short contigs included): "
            "on fragmented input an element cut by a breakpoint is only partly "
            "scored near the end. Default %(default)s, i.e. 5 kb at 10 kb windows; "
            "0 disables. Edges on a contig end are flagged in contig_edge either way."
        ),
    )
    experimental.add_argument(
        "--stride",
        type=int,
        default=None,
        help="Inference window stride in bp (default: the checkpoint's trained stride, 1000).",
    )
    advanced.add_argument(
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
    advanced.add_argument(
        "--short-contig-threshold",
        type=float,
        default=None,
        help="Override the calibrated seed/extend threshold for --short-contigs scan.",
    )
    experimental.add_argument(
        "--refine-window-len",
        type=int,
        default=DEFAULT_REFINE_WINDOW_LEN,
        help=(
            "EXPERIMENTAL fine-scale pass (default 0 = off). Re-score each call "
            "+/- --refine-margin at this window length, then split calls whose fine "
            "track shows a background gap (conflated elements) and place edges "
            "between the fine core and the coarse edge (--edge-alpha). Never drops "
            "a call. See docs/adaptive-multiscale-detection.md."
        ),
    )
    experimental.add_argument("--refine-stride", type=int, default=DEFAULT_REFINE_STRIDE,
                              help="Fine-pass stride in bp. Default %(default)s.")
    experimental.add_argument(
        "--refine-threshold",
        type=float,
        default=None,
        help=(
            "Fine-pass score threshold (default 0.95, calibrated at 5 kb for the "
            "tight-boundaries setting; 0.6-0.8 keeps recall unchanged)."
        ),
    )
    experimental.add_argument("--refine-margin", type=int, default=None,
                              help="bp re-scored either side of each call (default: --window-len).")
    experimental.add_argument(
        "--edge-alpha",
        type=float,
        default=DEFAULT_EDGE_ALPHA,
        help=(
            "Edge placement in [0, 1]: edge = fine_core + alpha * (coarse_edge - "
            "fine_core). 1 keeps coarse edges (overextend ~half a window); 0 uses "
            "the fine core (truncates). Default %(default)s."
        ),
    )
    experimental.add_argument(
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
    experimental.add_argument(
        "--split-min-gap-bp",
        type=int,
        default=DEFAULT_SPLIT_MIN_GAP_BP,
        help="Minimum fine-track gap (bp) that can split a call. Default %(default)s.",
    )
    experimental.add_argument(
        "--split-max-prob",
        type=float,
        default=DEFAULT_SPLIT_MAX_PROB,
        help=(
            "Split only across a gap whose fine windows average below this -- the "
            "mirror of --bridge-prob, so splitting cannot undo gap-filling. "
            "Default %(default)s."
        ),
    )
    experimental.add_argument(
        "--ensemble",
        choices=ENSEMBLE_MODES,
        default=DEFAULT_ENSEMBLE,
        help=(
            "How to combine multiple --checkpoint score tracks. 'max' "
            "(default) preserves borderline detections; 'mean' dilutes them "
            "below the seed threshold and measurably costs recall."
        ),
    )
