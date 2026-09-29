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
"""`prophound predict` mirrors predict_genome's flags without importing torch;
the two parsers must agree on every shared default and choice list."""

from prophound.cli import build_parser
from prophound.models.prophage.src.predict_genome import build_arg_parser


def _actions(parser):
    return {a.dest: a for a in parser._actions if a.dest != "help"}


def test_predict_subcommand_matches_predict_genome_defaults():
    top = build_parser()
    sub = next(a for a in top._actions if hasattr(a, "choices") and isinstance(a.choices, dict))
    cli, inner = _actions(sub.choices["predict"]), _actions(build_arg_parser())
    shared = sorted(set(cli) & set(inner))
    assert len(shared) > 20
    for dest in shared:
        assert cli[dest].default == inner[dest].default, dest
        if inner[dest].choices is not None:
            assert sorted(cli[dest].choices) == sorted(inner[dest].choices), dest
