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

from prophound.models.prophage.src.options import add_predict_arguments


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
    predict = subparsers.add_parser(
        "predict",
        help="Detect prophage-like regions in a nucleotide FASTA.",
        description=(
            "Score sliding windows of a genome FASTA with the DNA-only prophage "
            "model and emit windows.tsv, regions.tsv, and run_metadata.json."
        ),
    )
    add_predict_arguments(predict)
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
