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
"""The researcher's table and terminal list (src/report.py)."""

import pandas as pd

from prophound.models.prophage.src import report


def detailed(rows):
    cols = ["call_id", "genome_id", "contig_id", "start", "end", "boundary", "start_lo", "start_hi", "end_lo",
            "end_hi", "mean_prob", "max_prob", "contig_edge", "n_frac", "seed_support", "att_len", "refine_status"]
    return pd.DataFrame([dict(zip(cols, r)) for r in rows])


ROWS = [
    # call_id, genome, contig, start, end, boundary, s_lo, s_hi, e_lo, e_hi, mean, max, edge, n_frac,
    # support, att, refine
    ("prophound_1", "g", "c1", 1000, 21000, "model", 0, 10880, 11120, 22100, 0.80, 0.9990, "none", 0.0,
     3, None, "refined"),
    ("prophound_2", "g", "c1", 50000, 80000, "att", 50000, 50000, 80000, 80000, 0.60, 1.0, "none", 0.0,
     3, 40, "refined"),
    ("prophound_3", "g", "c2", 0, 30000, "model", 0, 0, 15180, 31650, 0.70, 1.0, "start", 0.25,
     2, None, "refined"),
]


def test_rank_is_max_prob_then_mean_prob():
    d = detailed(ROWS)
    assert list(report.rank_calls(d)) == [3, 2, 1]   # 1.0/0.70 first, 1.0/0.60 next, 0.999 last


def test_rank_is_per_genome():
    d = pd.concat([detailed(ROWS), detailed(ROWS).assign(genome_id="h")], ignore_index=True)
    assert list(report.rank_calls(d)) == [3, 2, 1, 3, 2, 1]


def test_simple_table_is_one_based_and_sorted_by_rank():
    t = report.simple_table(detailed(ROWS), n_checkpoints=3)
    assert list(t.columns) == report.SIMPLE_COLUMNS and "contig_id" not in t.columns
    assert list(t.call_id) == ["prophound_3", "prophound_2", "prophound_1"]
    top, att, last = t.iloc[0], t.iloc[1], t.iloc[2]
    assert (last.start, last.end, last.length_kb) == (1001, 21000, 20.0)    # 0-based 1000 -> 1-based 1001
    assert (last.start_range, last.end_range) == ("1-10881", "11120-22100")
    assert (att.edges, att.start_range, att.end_range) == ("att site", "50001", "80000")
    assert att.note == "att site: 40 bp repeat at both edges"
    assert top.support == "2/3" and top.edges == "estimated"
    assert top.note == "runs off the contig start; the element may continue; 25% N"


def test_batch_table_names_the_genome_and_ranks_within_it():
    d = pd.concat([detailed(ROWS), detailed(ROWS[:1]).assign(genome_id="h")], ignore_index=True)
    t = report.simple_table(d, n_checkpoints=3, with_genome=True)
    assert t.columns[0] == "genome" and list(t.genome) == ["g", "g", "g", "h"] and list(t["rank"]) == [1, 2, 3, 1]


def test_empty_run():
    empty = detailed([]).iloc[0:0]
    assert list(report.simple_table(empty, 3).columns) == report.SIMPLE_COLUMNS
    assert report.terminal_lines(report.simple_table(empty, 3), "default")[0].startswith("No candidate")


def test_terminal_lines_list_the_best_calls_and_say_how_to_read_them():
    lines = report.terminal_lines(report.simple_table(detailed(ROWS), 3), "default", limit=2)
    assert lines[0] == "3 candidate prophage regions (default mode), best first:"
    assert "1,-30,000" not in "\n".join(lines) and "1-30,000" in lines[2]
    assert lines[4] == "  ... and 1 more in regions.tsv"
    assert report.GUIDE_URL in lines[-1]


def test_note_flags_edges_nothing_moved_off_the_coarse_edge():
    d = detailed(ROWS[:1]).assign(start_source="coarse", end_source="composition")
    assert report.simple_table(d, 3).iloc[0].note == (
        "start not refined: the score stays high past it, so the true edge is less certain")
    both = detailed(ROWS[:1]).assign(start_source="coarse", end_source="coarse")
    assert report.simple_table(both, 3).iloc[0].note.startswith(
        "start and end not refined: the score stays high past them")
    moved = detailed(ROWS[:1]).assign(start_source="fine", end_source="composition")
    assert report.simple_table(moved, 3).iloc[0].note == ""
