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
DEFAULT_COVERAGE_FRAC = 0.3


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
    parser.add_argument("--checkpoint", required=True, help="Path to a trained model checkpoint.")
    parser.add_argument("--threshold", type=float, default=None, help="Override the mode's score threshold.")
    parser.add_argument("--seed-threshold", type=float, default=None)
    parser.add_argument("--extend-threshold", type=float, default=None)
    parser.add_argument(
        "--device",
        default="auto",
        choices=["auto", "cpu", "cuda", "rocm", "mps", "xpu", "intel"],
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
        default="strict",
        help=(
            "How per-window scores become regions. 'strict' (default) is the "
            "frozen v1 builder: a bin is inside only if the mean of all windows "
            "covering it >= the extend threshold, which truncates regions and drops "
            "sparse-signal elements. 'coverage' marks a bin inside if at least "
            "--coverage-frac of its covering windows individually score >= the "
            "extend threshold -- much higher recall and boundary accuracy, more "
            "calls per genome. See docs/multi-model-architecture.md. Not combinable "
            "with --num-tiles > 1."
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
