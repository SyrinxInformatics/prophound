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
"""`prophound predict` and `python -m ...predict_genome` share one parser (options.py)."""

import subprocess
import sys

import pytest

from prophound.cli import build_parser
from prophound.models.prophage.src import options


def _subparsers(parser):
    action = next(a for a in parser._actions if hasattr(a, "choices") and isinstance(a.choices, dict))
    return action.choices


def test_help_renders_for_every_subcommand(capsys):
    """argparse formats help lazily, so a stray `%` in help text only fails here."""
    parser = build_parser()
    for name in [None, *_subparsers(parser)]:
        with pytest.raises(SystemExit) as exc:
            parser.parse_args(["--help"] if name is None else [name, "--help"])
        assert exc.value.code == 0
    assert "--input" in capsys.readouterr().out


def test_help_does_not_import_torch():
    code = ("import sys; from prophound.cli import build_parser; "
            "build_parser().format_help(); print('torch' in sys.modules)")
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True).stdout
    assert out.strip() == "False"


def test_predict_defaults_come_from_options():
    args = build_parser().parse_args(["predict", "--input", "x.fna", "--outdir", "o", "--config", "c.yaml",
                                      "--checkpoint", "a.pt"])
    assert args.mode == "default"
    assert args.coverage_frac == options.DEFAULT_COVERAGE_FRAC
    assert args.max_gap_bp == options.DEFAULT_MAX_GAP_BP
    assert args.edge_snap_frac == options.DEFAULT_EDGE_SNAP_FRAC
    assert args.ambiguity_policy == options.DEFAULT_AMBIGUITY_POLICY
    assert sorted(_subparsers(build_parser())["predict"]._option_string_actions) == sorted(
        _standalone_parser()._option_string_actions)


def _standalone_parser():
    from prophound.models.prophage.src.predict_genome import build_arg_parser
    return build_arg_parser()


def test_standalone_main_accepts_input_fasta(monkeypatch):
    from prophound.models.prophage.src import predict_genome
    seen = {}
    monkeypatch.setattr(predict_genome, "run_predict", lambda args: seen.update(vars(args)))
    for argv, expected in ((["--input_fasta", "g.fna"], "g.fna"), (["--input_fasta=h.fna"], "h.fna")):
        with pytest.raises(SystemExit) as exc:
            predict_genome.main([*argv, "--outdir", "o", "--config", "c.yaml", "--checkpoint", "a.pt"])
        assert exc.value.code == 0 and seen["input_fasta"] == expected
