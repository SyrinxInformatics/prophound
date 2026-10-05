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
"""Edge snapping: a call edge near a contig end moves onto it, by a distance set
from --window-len, and edges on a contig end are flagged whether or not they moved."""

import pandas as pd

from prophound.models.prophage.src.predict_genome import DEFAULT_EDGE_SNAP_FRAC, snap_contig_edges


def calls(rows):
    df = pd.DataFrame(rows, columns=["contig_id", "start", "end"])
    return df.assign(region_len_bp=df.end - df.start)


def test_snap_distance_follows_window_length():
    L = {"a": 100_000, "s": 8_000}
    r = calls([("a", 4_000, 50_000), ("a", 60_000, 94_000), ("a", 9_000, 30_000), ("s", 3_000, 4_000)])
    out10, moved10 = snap_contig_edges(r, L, int(DEFAULT_EDGE_SNAP_FRAC * 10_000))   # 5 kb
    # a contig no longer than 2 x the snap distance: any call on it spans the contig
    assert list(zip(out10.start, out10.end)) == [(0, 50_000), (60_000, 94_000), (9_000, 30_000), (0, 8_000)]
    assert moved10 == 2
    assert list(out10.contig_edge) == ["start", "none", "none", "both"]
    assert (out10.region_len_bp == out10.end - out10.start).all()
    out20, moved20 = snap_contig_edges(r, L, int(DEFAULT_EDGE_SNAP_FRAC * 20_000))   # 10 kb
    assert list(zip(out20.start, out20.end)) == [(0, 50_000), (60_000, 100_000), (0, 30_000), (0, 8_000)]
    assert moved20 == 4
    assert list(out20.contig_edge) == ["start", "end", "start", "both"]


def test_off_still_flags_edges_already_on_a_contig_end():
    r = calls([("a", 0, 100_000), ("a", 2_000, 30_000)])
    out, moved = snap_contig_edges(r, {"a": 100_000}, 0)
    assert moved == 0
    assert list(zip(out.start, out.end)) == [(0, 100_000), (2_000, 30_000)]
    assert list(out.contig_edge) == ["both", "none"]


def test_empty_calls():
    out, moved = snap_contig_edges(calls([]), {}, 5_000)
    assert moved == 0 and "contig_edge" in out.columns and out.empty
