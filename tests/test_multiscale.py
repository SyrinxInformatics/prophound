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
"""Fine-scale refinement (multiscale.py): the four downstream behaviours.

Synthetic score tracks only -- a window scores high when enough of it overlaps
an element -- so these pin the logic, not the model.
"""

import pandas as pd
import pytest

from mge_detect.models.prophage.src.multiscale import (
    plan_fine_windows,
    plan_neighborhoods,
    refine_regions,
)
from mge_detect.models.prophage.src.predict_genome import build_region_table_coverage

G, C, L = "g", "c1", 200_000
FINE_LEN, FINE_STRIDE = 3000, 250


def coarse(start, end, kind="coverage_frac_0.1"):
    return pd.DataFrame([{"genome_id": G, "contig_id": C, "start": start, "end": end,
                          "n_windows": 10, "region_len_bp": end - start, "mean_prob": 0.99,
                          "max_prob": 1.0, "boundary_kind": kind}])


def fine_track(elements, lo, hi, gap_prob=0.01, min_frac=0.5):
    """Fine windows over [lo, hi): p=1 if >= min_frac of the window is element,
    else `gap_prob` (the background level between elements)."""
    rows = []
    for s in range((lo // FINE_STRIDE) * FINE_STRIDE, hi - FINE_LEN + 1, FINE_STRIDE):
        e = s + FINE_LEN
        ov = max((min(e, b) - max(s, a) for a, b in elements), default=0)
        rows.append({"genome_id": G, "contig_id": C, "start": s, "end": e, "fasta_path": "x",
                     "prob": 1.0 if ov >= min_frac * FINE_LEN else gap_prob})
    return pd.DataFrame(rows)


def run(coarse_df, fine_df, alpha=1.0, gap=2000, max_prob=0.1):
    fine_df, fine_regions = build_region_table_coverage(
        fine_df, seed_threshold=0.995, extend_threshold=0.995, min_coverage_frac=0.1, max_gap_bp=0)
    return refine_regions(coarse_df, fine_df, fine_regions, edge_alpha=alpha,
                          split_min_gap_bp=gap, split_max_prob=max_prob)


def overlaps(row, a, b):
    return min(row.end, b) > max(row.start, a)


# --- conflation: two small elements under one coarse call ------------------

def test_split_separates_two_elements_with_background_gap():
    els = [(50_000, 56_000), (64_000, 70_000)]           # 6 kb each, 8 kb apart
    out = run(coarse(45_000, 75_000), fine_track(els, 35_000, 85_000))
    assert list(out.refine_status) == ["split", "split"]
    assert overlaps(out.iloc[0], *els[0]) and not overlaps(out.iloc[0], *els[1])
    assert overlaps(out.iloc[1], *els[1]) and not overlaps(out.iloc[1], *els[0])


def test_split_refused_when_gap_is_elevated_so_gap_filling_survives():
    # a dip inside one large element: the gap still scores above split_max_prob,
    # which is exactly the case the coarse merge bridged -- must not be undone
    els = [(50_000, 56_000), (64_000, 70_000)]
    out = run(coarse(45_000, 75_000), fine_track(els, 35_000, 85_000, gap_prob=0.3))
    assert len(out) == 1 and out.iloc[0].refine_status == "refined"


def test_split_refused_when_gap_is_shorter_than_minimum():
    els = [(50_000, 56_000), (57_000, 63_000)]           # 1 kb apart
    out = run(coarse(45_000, 68_000), fine_track(els, 35_000, 78_000), gap=5000)
    assert len(out) == 1


# --- truncation / overextension: edge placement ----------------------------

def test_alpha_one_without_split_reproduces_coarse_call():
    out = run(coarse(45_000, 65_000), fine_track([(50_000, 60_000)], 35_000, 75_000), alpha=1.0)
    assert (out.iloc[0].start, out.iloc[0].end) == (45_000, 65_000)


def test_edges_interpolate_between_fine_core_and_coarse_edge():
    c = coarse(45_000, 65_000)
    f = fine_track([(50_000, 60_000)], 35_000, 75_000)
    core = run(c, f, alpha=0.0).iloc[0]
    half = run(c, f, alpha=0.5).iloc[0]
    assert 45_000 <= core.start and core.end <= 65_000           # never wider than coarse
    assert half.start == round(core.start - 0.5 * (core.start - 45_000))
    assert half.end == round(core.end + 0.5 * (65_000 - core.end))
    assert half.start <= core.start and half.end >= core.end     # never narrower than core


def test_alpha_tightens_small_element_iou():
    el = (50_000, 56_000)
    c = coarse(45_000, 61_000)                             # ~half a coarse window each side
    f = fine_track([el], 35_000, 71_000)

    def iou(r):
        inter = max(0, min(r.end, el[1]) - max(r.start, el[0]))
        return inter / (max(r.end, el[1]) - min(r.start, el[0]))

    assert iou(run(c, f, alpha=0.3).iloc[0]) > iou(run(c, f, alpha=1.0).iloc[0])


# --- recall: nothing is ever dropped ---------------------------------------

def test_call_without_fine_support_is_kept_and_flagged():
    out = run(coarse(45_000, 65_000), fine_track([], 35_000, 75_000))
    r = out.iloc[0]
    assert (r.start, r.end, r.refine_status) == (45_000, 65_000, "no_fine_support")


def test_edge_alpha_out_of_range_rejected():
    with pytest.raises(ValueError):
        run(coarse(45_000, 65_000), fine_track([(50_000, 60_000)], 35_000, 75_000), alpha=1.5)


# --- fine-window planning --------------------------------------------------

def test_fine_windows_are_grid_aligned_inside_contig_and_deduplicated():
    regs = pd.concat([coarse(1_000, 8_000), coarse(9_000, 15_000), coarse(190_000, 199_500)])
    nb = plan_neighborhoods(regs, {C: L}, margin_bp=10_000)
    assert nb[C] == [(0, 25_000), (180_000, L)]            # overlapping neighborhoods merged
    w = plan_fine_windows(nb, {C: L}, G, "x", FINE_LEN, FINE_STRIDE)
    assert (w.start % FINE_STRIDE == 0).all()
    assert (w.start >= 0).all() and (w.end <= L).all()
    assert not w.duplicated(["contig_id", "start"]).any()


def test_contig_shorter_than_fine_window_gets_no_windows():
    nb = {C: [(0, 2_000)]}
    assert plan_fine_windows(nb, {C: 2_000}, G, "x", FINE_LEN, FINE_STRIDE).empty
