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


# Must match mge_detect.models.prophage.src.predict_genome.MODE_THRESHOLDS.
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
        help="Confidence mode (controls the default score threshold).",
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
        "--device",
        default="auto",
        choices=["auto", "cpu", "cuda", "rocm", "mps", "xpu", "intel"],
    )
    parser.add_argument(
        "--ambiguity-policy",
        choices=("n-channel", "uniform", "impute"),
        default="n-channel",
        help=(
            "How IUPAC ambiguity codes (N, R, Y, ...) are encoded. 'n-channel' (default) "
            "is the training-time encoding and is what every published number used. "
            "'uniform' spreads each ambiguous base as 0.25 across A/C/G/T instead. "
            "EXPERIMENTAL: under 'n-channel' two contiguous ambiguity codes can take a "
            "confident window from 1.000 to 0.000; 'uniform' avoids that and is a no-op "
            "on ambiguity-free sequence, but it changes scores wherever such codes occur."
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


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mge-detect",
        description=(
            "Database-independent detection of mobile genetic elements in "
            "prokaryotic sequences."
        ),
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    _add_predict_subparser(subparsers)
    return parser


def main() -> None:
    args = build_parser().parse_args()

    if args.command == "predict":
        from mge_detect.models.prophage.src.predict_genome import run_predict

        run_predict(args)
    else:  # pragma: no cover - argparse enforces valid choices
        print(f"Unknown command: {args.command}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
