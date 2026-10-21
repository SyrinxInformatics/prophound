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
import json
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
# The fine refinement pass, on by default since 0.1.0 (docs/att-boundaries.md "The combined
# pipeline"; the setting is from docs/downstream-postprocessing.md). 0 turns it off.
DEFAULT_REFINE_WINDOW_LEN = 5000
DEFAULT_REFINE_STRIDE = 250
DEFAULT_REFINE_THRESHOLD = 0.6
DEFAULT_EDGE_ALPHA = 0.6        # 1 keeps the coarse edges
DEFAULT_EDGE_ALPHA_RAMP = (15000, 60000)
DEFAULT_SPLIT_MIN_GAP_BP = 4000
# Below DEFAULT_BRIDGE_PROB, so splitting cannot undo a merge.
DEFAULT_SPLIT_MAX_PROB = 0.05
# A call at least this fraction N is dropped: `impute` turns a scaffold gap into random
# sequence, which can be called (docs/release-scope-0.1.0.md A.3).
MAX_CALL_N_FRAC = 0.5
# A call edge within this fraction of --window-len of a contig end moves onto the end:
# an element cut by a breakpoint is only partly scored there
# (docs/claims-and-evidence.md "G5b"). 0 turns snapping off.
DEFAULT_EDGE_SNAP_FRAC = 0.5

# The model shipped with the package: weights-only checkpoints, model.yaml and a manifest of
# SHA-256s (made by scripts/strip_checkpoint.py from the training checkpoints).
BUNDLED_MODEL_DIR = Path(__file__).resolve().parents[1] / "model"


def bundled_model() -> tuple[str, tuple[str, ...], dict]:
    """The bundled model's config path, checkpoint paths and manifest."""
    manifest = json.loads((BUNDLED_MODEL_DIR / "manifest.json").read_text())
    return (str(BUNDLED_MODEL_DIR / manifest["config"]),
            tuple(str(BUNDLED_MODEL_DIR / c["file"]) for c in manifest["checkpoints"]), manifest)


DEVICES = ("auto", "cpu", "cuda", "rocm", "mps", "xpu", "intel")
DEFAULT_BATCH_SIZE = 128


class PredictError(Exception):
    """A failure the user can fix: printed as one `error:` line, without a traceback."""


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
    att_search: bool = True         # att boundaries where found (docs/att-boundaries.md)
    composition_edges: str = "auto"  # auto = on in Default and Conservative
    verbosity: int = 0              # -1 quiet, 0 normal, 1 verbose
    force: bool = False             # replace outputs already in outdir
    bundled_model: bool = False     # checkpoints are the ones shipped with the package
    trust_checkpoint: bool = False  # allow full (pickled) training checkpoints

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
            raise PredictError("--window-len and --stride must be positive")
        ramp = (int(args.edge_alpha_ramp[0]), int(args.edge_alpha_ramp[1]))
        if ramp == (0, 0):   # turns the ramp off
            ramp = None
        elif not 0 < ramp[0] < ramp[1]:
            raise PredictError("--edge-alpha-ramp needs 0 < LO_BP < HI_BP (or 0 0 for no ramp)")
        refine_threshold = DEFAULT_REFINE_THRESHOLD if args.refine_threshold is None else float(args.refine_threshold)
        refine_margin = window_len if args.refine_margin is None else int(args.refine_margin)
        refine_window_len = int(args.refine_window_len or 0)
        if refine_window_len:
            if refine_window_len >= window_len:
                raise PredictError("--refine-window-len must be smaller than --window-len")
            if args.refine_stride <= 0 or not 0.0 <= args.edge_alpha <= 1.0 or refine_margin < 0:
                raise PredictError("need --refine-stride > 0, 0 <= --edge-alpha <= 1, --refine-margin >= 0")
        if not 0.0 <= args.edge_snap_frac <= 1.0:
            raise PredictError("--edge-snap-frac must be between 0 and 1")
        min_support = MODE_MIN_SEED_SUPPORT[args.mode] if args.min_seed_support is None else int(args.min_seed_support)
        if min_support < 1:
            raise PredictError("--min-seed-support must be >= 1")
        bundled_config, bundled_checkpoints, _ = bundled_model()
        bundled = not args.checkpoint
        if bundled:
            checkpoints = bundled_checkpoints
        else:
            checkpoints = tuple(args.checkpoint) if isinstance(args.checkpoint, list) else (args.checkpoint,)
        config = args.config or bundled_config
        inputs = args.input_fasta if isinstance(args.input_fasta, list) else [args.input_fasta]
        return cls(
            input_fasta=str(Path(inputs[0])) if inputs and inputs[0] else "", outdir=str(args.outdir),
            mode=args.mode, config=config, checkpoints=checkpoints, device=args.device,
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
            att_search=getattr(args, "att_search", "on") == "on",
            composition_edges=getattr(args, "composition_edges", "auto"),
            verbosity=-1 if getattr(args, "quiet", False) else int(bool(getattr(args, "verbose", False))),
            force=bool(getattr(args, "force", False)),
            bundled_model=bundled, trust_checkpoint=bool(getattr(args, "trust_checkpoint", False)),
        )


PREDICT_EPILOG = """\
outputs (in --outdir):
  regions.tsv           start here: the candidates, best first, 1-based: rank, contig, start, end,
                        edges (att site = exact to the base; estimated = can be kb off), a range for
                        each end, model support and notes
  regions.detailed.tsv  the same calls with every column, 0-based like BED (for scripts)
  regions.filtered.tsv  candidates set aside: the detailed columns plus filter (seed_vote, n_gap)
  regions.bed/.gff3     the calls for genome browsers; BED score = 1000 x mean_prob
  windows.tsv           every window's score (prob; prob_ckpt1.. with --write-seed-scores)
  fine_windows.tsv      the refinement pass's window scores
  run_metadata.json     settings, hashes, counts, timing and the command
  regions.fna           the calls' sequences (--write-region-fasta)

a batch (several --input, a directory, or --input-list) writes these per genome in
outdir/<genome_id>/, plus at the top: combined regions.tsv, regions.detailed.tsv, regions.filtered.tsv,
regions.bed/.gff3 (named genome_id:call_id), batch_summary.tsv and run_metadata.json.

calls are candidates for follow-up, not annotations. How to read them:
https://github.com/SyrinxInformatics/mge-detect/blob/main/docs/interpreting-results.md
"""


def add_predict_arguments(parser: argparse.ArgumentParser) -> None:
    """Add every `predict` option to `parser`, in basic / advanced / experimental --help sections."""
    basic = parser.add_argument_group("basic options")
    advanced = parser.add_argument_group("advanced options (thresholds, ensemble, input handling, outputs)")
    experimental = parser.add_argument_group("experimental options (not validated for default use)")

    basic.add_argument(
        "--input", dest="input_fasta", action="append", metavar="PATH",
        help=("Genome to scan: FASTA, GenBank/RefSeq, EMBL, or GFF3 with an embedded ##FASTA section (as Prokka "
              "and Bakta write), optionally gzip, bzip2 or xz compressed; the format is read from the content. "
              "Repeat it, or give a directory, to scan several genomes in one batch (one subdirectory each)."),
    )
    basic.add_argument(
        "--input-list", metavar="FILE",
        help="A file listing genomes to scan as a batch, one path per line (relative to the file; # comments).",
    )
    basic.add_argument(
        "--outdir", required=True,
        help="Directory for the outputs (created if missing). A directory that already holds outputs is refused.",
    )
    basic.add_argument(
        "--force", action="store_true",
        help="Replace the outputs already in --outdir. Every file prophound writes there is removed first.",
    )
    basic.add_argument(
        "--mode",
        choices=sorted(MODE_THRESHOLDS),
        default="default",
        help=(
            "How many calls to make (default: default). exploratory: seed threshold 0.90, every call kept, "
            "about 32 calls per genome. default: 0.95, a call is kept only if all 3 checkpoints make it, "
            "about 11 per genome. conservative: 0.995, all 3 checkpoints, about 5 per genome. The call "
            "counts are from a held-out test on 150 complete bacterial genomes."
        ),
    )
    basic.add_argument(
        "--device",
        default="auto",
        choices=DEVICES,
        help=(
            "Where to run the model. auto (default) picks CUDA, then Apple MPS, then Intel XPU, then CPU. "
            "rocm runs on an AMD GPU through PyTorch's CUDA interface; intel is the same as xpu."
        ),
    )
    basic.add_argument(
        "--write-region-fasta",
        action="store_true",
        help="Also write the called regions' sequences to regions.fna.",
    )
    chatter = basic.add_mutually_exclusive_group()
    chatter.add_argument("-q", "--quiet", action="store_true",
                         help="Print only warnings, errors and the summary.")
    chatter.add_argument("-v", "--verbose", action="store_true",
                         help="Also print each stage's timing, and the full traceback of an unexpected error.")

    advanced.add_argument(
        "--checkpoint",
        nargs="+",
        default=None,
        help=(
            "Model checkpoint file(s) to use instead of the bundled model (its three checkpoints are the "
            "default). The modes' seed vote needs three; one checkpoint makes more calls."
        ),
    )
    advanced.add_argument(
        "--config", default=None,
        help="Model config (YAML) for --checkpoint files (default: the bundled model's).",
    )
    advanced.add_argument(
        "--trust-checkpoint", action="store_true",
        help=("Allow full training checkpoints, which load by unpickling arbitrary Python objects; use only "
              "for files you trust. Weights-only files (scripts/strip_checkpoint.py) need no flag."),
    )
    advanced.add_argument(
        "--threshold", type=float, default=None,
        help="Override the mode's seed threshold (sets --seed-threshold, and --extend-threshold through it).",
    )
    advanced.add_argument(
        "--seed-threshold", type=float, default=None,
        help="Score a window needs to start a call (default: the mode's threshold).",
    )
    advanced.add_argument(
        "--extend-threshold", type=float, default=None,
        help="Score windows need to extend a call (default: the seed threshold, at most 0.99).",
    )
    advanced.add_argument(
        "--min-seed-support",
        type=int,
        default=None,
        help=(
            "Keep a call only if at least this many checkpoints make an overlapping call on their own "
            "(default: the mode's; 1 for exploratory, 3 otherwise). 1 keeps every call."
        ),
    )
    advanced.add_argument(
        "--write-seed-scores",
        action="store_true",
        help="With several checkpoints, also write each one's window scores to windows.tsv (prob_ckpt1, ...).",
    )
    advanced.add_argument(
        "--ambiguity-policy",
        choices=AMBIGUITY_POLICIES,
        default=DEFAULT_AMBIGUITY_POLICY,
        help=(
            "How IUPAC ambiguity codes (N, R, Y, ...) reach the model. impute (default) replaces each with a "
            "fixed, reproducible base and changes nothing on sequence without them. n-channel is the "
            "training-time encoding, but two ambiguous bases in a row can drop a window's score from 1 to 0. "
            "uniform gives each ambiguous base 0.25 on A, C, G and T."
        ),
    )
    advanced.add_argument(
        "--batch-size", type=int, default=DEFAULT_BATCH_SIZE,
        help="Windows scored at once (default %(default)s). Lower it if the device runs out of memory.",
    )
    advanced.add_argument(
        "--num-workers", type=int, default=0,
        help="Extra processes reading windows from the FASTA (default 0: read in the main process).",
    )
    advanced.add_argument(
        "--coverage-frac",
        type=float,
        default=DEFAULT_COVERAGE_FRAC,
        help=(
            "Fraction of the windows covering a position that must reach the extend threshold for it to be "
            "inside a call (default %(default)s; 0 is the union of high windows, 1 the intersection)."
        ),
    )
    advanced.add_argument(
        "--max-gap-bp",
        type=int,
        default=DEFAULT_MAX_GAP_BP,
        help=(
            "Merge calls separated by at most this many bp when the gap still looks element-like (see "
            "--bridge-prob). This repairs large elements whose score dips partway through. "
            "Default %(default)s; 0 disables."
        ),
    )
    advanced.add_argument(
        "--bridge-prob",
        type=float,
        default=DEFAULT_BRIDGE_PROB,
        help=(
            "Merge across a gap only if its mean window score stays at least this high everywhere "
            "(default %(default)s). 0 merges on distance alone, which can join separate neighbouring elements."
        ),
    )
    advanced.add_argument(
        "--edge-snap-frac",
        type=float,
        default=DEFAULT_EDGE_SNAP_FRAC,
        help=(
            "Move a call edge within this fraction of --window-len of a contig end onto the end (default "
            "%(default)s, 5 kb at 10 kb windows; 0 disables). On a fragmented assembly an element cut by a "
            "contig break is only partly scored near the break. The contig_edge column flags calls that "
            "touch a contig end either way."
        ),
    )
    advanced.add_argument(
        "--short-contigs",
        choices=SHORT_CONTIG_MODES,
        default=DEFAULT_SHORT_CONTIGS,
        help=(
            "Contigs shorter than --window-len. scan (default) tiles each with 5 kb or 3 kb windows, "
            "whichever fits, at thresholds set to give the same false-call rate as full-length windows; "
            "skip leaves them unscanned. Contigs under 3 kb are never scanned. A warning reports them either way."
        ),
    )
    advanced.add_argument(
        "--short-contig-threshold",
        type=float,
        default=None,
        help="Override the short-contig seed and extend threshold.",
    )

    experimental.add_argument(
        "--ensemble",
        choices=ENSEMBLE_MODES,
        default=DEFAULT_ENSEMBLE,
        help=(
            "How to combine the checkpoints' window scores: max (default) or mean. mean pulls borderline "
            "windows under the seed threshold and costs recall."
        ),
    )
    experimental.add_argument(
        "--grid-offset",
        type=int,
        default=0,
        help="Shift the window grid by this many bp (default 0), to test how calls depend on window placement.",
    )
    experimental.add_argument(
        "--window-len",
        type=int,
        default=None,
        help=(
            "Window length in bp (default: the trained length, 10000). The mode thresholds were set at the "
            "trained length, and scores on non-phage sequence rise as windows shrink."
        ),
    )
    experimental.add_argument(
        "--stride",
        type=int,
        default=None,
        help="Window step in bp (default: the trained step, 1000).",
    )
    advanced.add_argument(
        "--att-search",
        choices=("on", "off"),
        default="on",
        help=(
            "Look for an att site at each call: a direct repeat with one copy near each edge, left by "
            "integration. Where one is found, the call's edges are its two copies, exact to the base "
            "(boundary = att). Default on."
        ),
    )
    advanced.add_argument(
        "--composition-edges",
        choices=("auto", "on", "off"),
        default="auto",
        help=(
            "Move each estimated edge to where the sequence's composition changes from host-like to "
            "element-like. auto (default): on in default and conservative mode, off in exploratory, where "
            "it cost recall."
        ),
    )
    advanced.add_argument(
        "--refine-window-len",
        type=int,
        default=DEFAULT_REFINE_WINDOW_LEN,
        help=(
            "Fine-scale pass (default %(default)s; 0 turns it off). Re-scores each call and its surroundings "
            "with windows of this length, splits calls whose fine scores show a background gap, and places "
            "edges between the fine and coarse estimates (--edge-alpha). Never drops a call."
        ),
    )
    advanced.add_argument("--refine-stride", type=int, default=DEFAULT_REFINE_STRIDE,
                              help="Fine-pass window step in bp (default %(default)s).")
    advanced.add_argument(
        "--refine-threshold",
        type=float,
        default=None,
        help=f"Fine-pass score threshold (default {DEFAULT_REFINE_THRESHOLD}).",
    )
    advanced.add_argument("--refine-margin", type=int, default=None,
                              help="bp re-scored on either side of each call (default: --window-len).")
    advanced.add_argument(
        "--edge-alpha",
        type=float,
        default=DEFAULT_EDGE_ALPHA,
        help=(
            "Where refined edges go, from 0 (the fine core, which truncates) to 1 (the coarse edge, which "
            "overextends by about half a window). Default %(default)s."
        ),
    )
    advanced.add_argument(
        "--edge-alpha-ramp",
        type=int,
        nargs=2,
        metavar=("LO_BP", "HI_BP"),
        default=DEFAULT_EDGE_ALPHA_RAMP,
        help=(
            "Scale --edge-alpha with size: 0 for pieces whose fine core is at most LO_BP, --edge-alpha at "
            "HI_BP and above, log-linear between (default 15000 60000; 0 0 for no ramp). Small elements suit "
            "the fine core, large ones the coarse edge."
        ),
    )
    advanced.add_argument(
        "--split-min-gap-bp",
        type=int,
        default=DEFAULT_SPLIT_MIN_GAP_BP,
        help="Shortest fine-score gap (bp) that can split a call (default %(default)s).",
    )
    advanced.add_argument(
        "--split-max-prob",
        type=float,
        default=DEFAULT_SPLIT_MAX_PROB,
        help=(
            "Split only across a gap whose fine windows average below this (default %(default)s). It sits "
            "under --bridge-prob, so splitting cannot undo a merge."
        ),
    )
