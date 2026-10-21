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
"""Att boundaries and edge ranges (src/att.py, docs/att-boundaries.md)."""

import numpy as np
import pandas as pd
import pytest

from prophound.models.prophage.src import att

CONTIG = "c1"
LEN = 200_000


def random_seq(n: int, seed: int = 0) -> list:
    return list(np.random.default_rng(seed).choice(list("ACGT"), size=n))


def plant(seq: list, pos: int, core: str) -> None:
    seq[pos:pos + len(core)] = list(core)


def make_maximal(seq: list, a: int, b: int, length: int) -> None:
    """Make the bases flanking two copies differ, so the repeat is exactly `length` long."""
    for x, y in ((a - 1, b - 1), (a + length, b + length)):
        if seq[x] == seq[y]:
            seq[y] = {"A": "C", "C": "G", "G": "T", "T": "A"}[seq[y]]


CORE = "GATTCGCTAGCATGCAATCGGATCCAGTTGACCATGGCTA"   # 40 bp, high entropy


def genome_with_pair(core: str = CORE, a: int = 80_000, b: int = 110_000, seed: int = 0) -> att.Genome:
    seq = random_seq(LEN, seed)
    plant(seq, a, core)
    plant(seq, b, core)
    make_maximal(seq, a, b, len(core))
    return att.Genome({CONTIG: "".join(seq)})


def test_planted_pair_gives_base_exact_boundaries():
    g = genome_with_pair()
    hit = att.find_att(g, CONTIG, 77_000, 113_000)   # coarse call overshoots by 3 kb each side
    assert hit is not None
    assert (hit["att_start"], hit["att_end"], hit["att_len"], hit["att_core"]) == (80_000, 110_040, 40, CORE)
    assert hit["evalue"] <= att.MAX_EVALUE


def test_random_sequence_has_no_att():
    g = att.Genome({CONTIG: "".join(random_seq(LEN, seed=3))})
    assert att.find_att(g, CONTIG, 77_000, 113_000) is None


def test_a_third_copy_rejects_the_repeat():
    seq = list(genome_with_pair().seqs[CONTIG])
    plant(seq, 180_000, CORE)
    assert att.find_att(att.Genome({CONTIG: "".join(seq)}), CONTIG, 77_000, 113_000) is None


def test_low_complexity_core_is_rejected():
    g = genome_with_pair(core="AT" * 20)
    assert att.find_att(g, CONTIG, 77_000, 113_000) is None


def test_core_longer_than_the_cap_is_rejected():
    core = "".join(random_seq(att.MAX_CORE_LEN + 50, seed=7))
    g = genome_with_pair(core=core)
    assert att.find_att(g, CONTIG, 77_000, 113_000) is None


def test_pair_outside_the_accepted_zones_is_ignored():
    # attL 12 kb inside the coarse start: inside the search zone, outside the accepted one
    g = genome_with_pair(a=89_000, b=110_000)
    assert att.find_att(g, CONTIG, 77_000, 113_000) is None


def test_equal_pairs_tie_to_the_leftmost():
    seq = random_seq(LEN)
    other = "CCTAGGATCAGTCAGGTACGATTGCACTGACGTTAGCATC"   # a second 40 bp core
    for core, a, b in ((CORE, 80_000, 110_000), (other, 79_000, 111_000)):
        plant(seq, a, core)
        plant(seq, b, core)
        make_maximal(seq, a, b, 40)
    hit = att.find_att(att.Genome({CONTIG: "".join(seq)}), CONTIG, 77_000, 113_000)
    assert (hit["att_start"], hit["att_core"]) == (79_000, other)


def calls_frame(rows):
    cols = ["genome_id", "contig_id", "start", "end", "region_len_bp", "contig_edge"]
    return pd.DataFrame([dict(zip(cols, r)) for r in rows])


def test_apply_boundaries_replaces_edges_and_sets_ranges():
    g = genome_with_pair()
    calls = calls_frame([("g", CONTIG, 77_000, 113_000, 36_000, "none"),
                         ("g", CONTIG, 150_000, 170_000, 20_000, "none")])
    out = att.apply_boundaries(calls, g, {CONTIG: LEN})
    a, m = out.iloc[0], out.iloc[1]
    assert (a.boundary, a.start, a.end, a.region_len_bp, a.model_start, a.model_end) == (
        "att", 80_000, 110_040, 30_040, 77_000, 113_000)
    assert (a.start_lo, a.start_hi, a.end_lo, a.end_hi) == (80_000, 80_000, 110_040, 110_040)
    assert m.boundary == "model" and m.att_core == "" and pd.isna(m.att_len)
    out_bp, in_bp = round(att.REL_OUT * 20_000), round(att.REL_IN * 20_000)
    assert (m.start_lo, m.start_hi, m.end_lo, m.end_hi) == (
        150_000 - out_bp, 150_000 + in_bp, 170_000 - in_bp, 170_000 + out_bp)


def test_contig_edge_calls_are_not_searched_and_their_end_edge_is_exact():
    g = genome_with_pair()
    calls = calls_frame([("g", CONTIG, 77_000, LEN, LEN - 77_000, "end")])
    out = att.apply_boundaries(calls, g, {CONTIG: LEN}).iloc[0]
    assert out.boundary == "model" and out.start == 77_000
    assert out.end_lo == out.end_hi == LEN
    assert out.start_lo == 77_000 - round(att.REL_OUT * (LEN - 77_000))


def test_ranges_are_clipped_to_the_contig():
    calls = calls_frame([("g", CONTIG, 100, 20_100, 20_000, "none")])
    out = att.apply_boundaries(calls, None, {CONTIG: 20_500}).iloc[0]
    assert out.start_lo == 0 and out.end_hi == 20_500


@pytest.mark.parametrize("col", ["short", "boundary_kind"])
def test_short_contig_calls_are_not_searched(col):
    g = genome_with_pair()
    calls = calls_frame([("g", CONTIG, 77_000, 113_000, 36_000, "none")])
    calls[col] = True if col == "short" else "coverage_short5000"
    assert att.apply_boundaries(calls, g, {CONTIG: LEN}).iloc[0].boundary == "model"


def test_split_pieces_are_not_searched():
    g = genome_with_pair()
    calls = calls_frame([("g", CONTIG, 77_000, 113_000, 36_000, "none")]).assign(
        coarse_start=77_000, coarse_end=113_000, n_pieces=2)
    assert att.apply_boundaries(calls, g, {CONTIG: LEN}).iloc[0].boundary == "model"


def test_search_uses_the_coarse_coordinates():
    """A refined call is searched where the rule was tuned: its coarse edges."""
    g = genome_with_pair()
    calls = calls_frame([("g", CONTIG, 82_000, 108_000, 26_000, "none")]).assign(
        coarse_start=77_000, coarse_end=113_000, n_pieces=1)
    out = att.apply_boundaries(calls, g, {CONTIG: LEN}).iloc[0]
    assert (out.boundary, out.start, out.end, out.model_start) == ("att", 80_000, 110_040, 82_000)
