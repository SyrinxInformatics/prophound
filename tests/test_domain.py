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
"""The non-prokaryote warning (src/domain.py, docs/robustness-checks.md check 8) and its override."""

import random

from prophound.cli import build_parser
from prophound.models.prophage.src import domain
from prophound.models.prophage.src.options import PredictOptions

W = domain.TELOMERE_WINDOW


def _random(n, seed=0):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def _contig(seq):
    return (len(seq), seq[:W], seq[-W:])


def test_bacterial_like_input_passes():
    assert domain.non_prokaryote_reasons([_contig(_random(50_000)), _contig(_random(20_000, 1))]) == []


def test_telomeres_at_two_ends_warn():
    body = _random(30_000)
    seq = "CCCTAA" * 30 + body + "TTAGGG" * 30
    reasons = domain.non_prokaryote_reasons([_contig(seq)])
    assert len(reasons) == 1 and "telomere repeats at 2 of 2 contig ends" in reasons[0]


def test_one_telomeric_end_is_not_enough():
    seq = "TTAGGG" * 30 + _random(30_000)
    assert domain.non_prokaryote_reasons([_contig(seq)]) == []


def test_short_contigs_are_not_searched():
    seq = "CCCTAA" * 30 + _random(2_000) + "TTAGGG" * 30   # under TELOMERE_MIN_CONTIG
    assert domain.non_prokaryote_reasons([_contig(seq)]) == []


def test_lowercase_and_yeast_style_telomeres_count():
    seq = ("tgtgggtgtggtgtgtggg" * 3) + _random(30_000) + ("cccacacacccacaccca" * 3)
    assert domain.non_prokaryote_reasons([_contig(seq)])


def test_size_over_the_largest_prokaryote_warns():
    reasons = domain.non_prokaryote_reasons([(18_000_000, "A" * W, "C" * W)])
    assert len(reasons) == 1 and "18.0 Mb" in reasons[0]
    assert domain.non_prokaryote_reasons([(16_100_000, "A" * W, "C" * W)]) == []


def _opts(*extra):
    args = build_parser().parse_args(["predict", "--input", "x.fna", "--outdir", "o", *extra])
    return PredictOptions.from_args(args, 10_000, 1_000)


def test_check_is_on_by_default_and_the_flag_skips_it():
    assert _opts().domain_check is True
    assert _opts("--no-domain-check").domain_check is False
