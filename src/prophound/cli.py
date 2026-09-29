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


import argparse
import sys


# Must match prophound.models.prophage.src.predict_genome.MODE_THRESHOLDS.
# Kept as a local constant so `--help` doesn't have to import torch.
PROPHAGE_MODES = ("conservative", "default", "exploratory")
# Must match predict_genome.REGION_BUILDERS / DEFAULT_COVERAGE_FRAC, same reason.
REGION_BUILDERS = ("strict", "coverage")
ENSEMBLE_MODES = ("max", "mean")
DEFAULT_ENSEMBLE = "max"
# Must match predict_genome's DEFAULT_* block. Defaults changed 2026-09-21 to
# the measured operating point -- docs/model-quality-baseline.md.
DEFAULT_REGION_BUILDER = "coverage"
DEFAULT_COVERAGE_FRAC = 0.1
DEFAULT_MAX_GAP_BP = 50000
DEFAULT_BRIDGE_PROB = 0.1
# mirrors predict_genome's DEFAULT_REFINE_* / DEFAULT_EDGE_ALPHA / DEFAULT_SPLIT_*
DEFAULT_REFINE_WINDOW_LEN = 0
DEFAULT_REFINE_STRIDE = 250
DEFAULT_EDGE_ALPHA = 1.0
DEFAULT_SPLIT_MIN_GAP_BP = 4000
DEFAULT_SPLIT_MAX_PROB = 0.05


def _add_predict_subparser(subparsers: "argparse._SubParsersAction") -> None:
    parser = subparsers.add_parser(
        "predict",
        help="Detect prophage-like regions in a nucleotide FASTA.",
        description=(
            "Score sliding windows of a genome FASTA with the DNA-only prophage "
            "model and emit windows.tsv, regions.tsv, and run_metadata.json."
        ),
    )
    parser.add_argument(
        "--input", required=True, dest="input_fasta", help="Nucleotide FASTA to scan."
    )
    parser.add_argument("--outdir", required=True, help="Directory to write outputs into.")
    parser.add_argument(
        "--mode",
        choices=PROPHAGE_MODES,
        default="default",
        help=(
            "Operating mode: seed threshold plus seed vote. exploratory 0.90, every "
            "call kept (~36 calls/genome); default 0.95, all 3 checkpoints agree "
            "(~14); conservative 0.995, all 3 agree (~8)."
        ),
    )
    parser.add_argument(
        "--config",
        required=True,
        help="Training config YAML for the checkpoint (defines model/window architecture).",
    )
    parser.add_argument(
        "--checkpoint",
        required=True,
        nargs="+",
        help=(
            "Trained checkpoint(s). Passing more than one runs a seed ensemble "
            "over them (see --ensemble) -- recommended, since it removes the "
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
    parser.add_argument("--threshold", type=float, default=None, help="Override the mode's score threshold.")
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
            "every call."
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
            "windows.tsv (prob_ckpt1, ...) for offline re-analysis."
        ),
    )
    parser.add_argument(
        "--device",
        default="auto",
        choices=["auto", "cpu", "cuda", "rocm", "mps", "xpu", "intel"],
    )
    parser.add_argument(
        "--ambiguity-policy",
        choices=("n-channel", "uniform", "impute"),
        default="impute",
        help=(
            "How IUPAC ambiguity codes (N, R, Y, ...) are encoded. 'impute' (default "
            "since 2026-09-24) replaces each with a deterministic concrete base and is "
            "a no-op on ambiguity-free sequence. 'n-channel' is the training-time "
            "encoding and reproduces numbers published before then, but two contiguous "
            "ambiguity codes can take a confident window from 1.000 to 0.000. 'uniform' "
            "spreads each ambiguous base as 0.25 across A/C/G/T."
        ),
    )
    parser.add_argument(
        "--grid-offset",
        type=int,
        default=0,
        help=(
            "Shift the inference window grid by this many bp. Default 0. For measuring "
            "phase sensitivity: marginal calls can appear or vanish when the grid moves "
            "(see docs/scientific-validity.md 'Input-robustness standing check')."
        ),
    )
    parser.add_argument("--batch-size", type=int, default=128)
    parser.add_argument("--num-workers", type=int, default=0)
    parser.add_argument(
        "--num-tiles",
        type=int,
        default=1,
        help=(
            "EXPERIMENTAL, not recommended: score this many phase-shifted "
            "window tilings and merge them for finer region boundaries "
            "(default 1 = single grid, matches the frozen v1 behavior). Has "
            "shown real boundary-accuracy gains but also an unexplained "
            "region-recall regression and at least one unresolved seeding "
            "discrepancy vs. the single-grid baseline -- see "
            "docs/boundary-refinement.md before using this for anything but "
            "further investigation. Cost scales ~linearly with this value; "
            "also writes an extra bins.tsv diagnostic file when > 1."
        ),
    )
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
        choices=("skip", "scan"),
        default="scan",
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


# --help grouping for `predict` (docs/claims-and-evidence.md usability audit):
# everything not listed here is "advanced".
BASIC_OPTIONS = ("input_fasta", "outdir", "mode", "config", "checkpoint", "device", "write_region_fasta")
EXPERIMENTAL_OPTIONS = ("num_tiles", "grid_offset", "window_len", "stride", "refine_window_len", "refine_stride",
                        "refine_threshold", "refine_margin", "edge_alpha", "edge_alpha_ramp", "split_min_gap_bp",
                        "split_max_prob", "region_builder", "ensemble")


def _group_help(parser: argparse.ArgumentParser) -> None:
    """Re-file `predict`'s options into basic / advanced / experimental sections of --help.

    Only the help layout changes: every action keeps its definition and default.
    """
    groups = {"basic": parser.add_argument_group("basic options"),
              "advanced": parser.add_argument_group("advanced options (thresholds, ensemble, input handling, outputs)"),
              "experimental": parser.add_argument_group(
                  "experimental options (not validated for default use; see docs/)")}
    for grp in parser._action_groups:
        if grp in groups.values():
            continue
        for action in list(grp._group_actions):
            if action.dest == "help":
                continue
            kind = ("basic" if action.dest in BASIC_OPTIONS else
                    "experimental" if action.dest in EXPERIMENTAL_OPTIONS else "advanced")
            grp._group_actions.remove(action)
            groups[kind]._group_actions.append(action)


def build_parser() -> argparse.ArgumentParser:
    from prophound import __version__
    parser = argparse.ArgumentParser(
        prog="prophound",
        description=(
            "Recall-first, database-free prophage screen for prokaryotic "
            "genomes. Calls are candidates for follow-up, not annotations."
        ),
    )
    parser.add_argument("--version", action="version", version=f"prophound {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)
    _add_predict_subparser(subparsers)
    _group_help(subparsers.choices["predict"])
    return parser


def main() -> None:
    args = build_parser().parse_args()

    if args.command == "predict":
        from prophound.models.prophage.src.predict_genome import run_predict

        run_predict(args)
    else:  # pragma: no cover - argparse enforces valid choices
        print(f"Unknown command: {args.command}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
