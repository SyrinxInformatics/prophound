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
"""The shared call pipeline (calls.py): filters mark candidates instead of deleting them,
and the offline rules are the CLI's."""

import numpy as np
import pandas as pd
import pytest

from prophound.models.prophage.src.calls import (
    FILTER_N_GAP,
    FILTER_SEED_VOTE,
    CallRules,
    build_pass,
    build_region_table_coverage,
    concat_nonempty,
    finalize,
)

CONTIG_LEN = 200_000


def track(*high_spans):
    """10 kb windows at 1 kb stride over one contig; prob 0.999 for windows starting in a span."""
    starts = np.arange(0, CONTIG_LEN - 10_000 + 1, 1_000)
    prob = np.full(len(starts), 0.01)
    for lo, hi in high_spans:
        prob[(starts >= lo) & (starts < hi)] = 0.999
    return pd.DataFrame({"genome_id": "g", "contig_id": "c", "start": starts, "end": starts + 10_000, "prob": prob})


RULES = CallRules(seed_threshold=0.95, extend_threshold=0.95, min_seed_support=3, max_gap_bp=0, snap_bp=5_000)
A, B = (40_000, 60_000), (120_000, 140_000)


def two_call_pass():
    """Call A is made by all three checkpoints, call B by one only."""
    w = track(A, B)
    per = [track(A, B).prob.to_numpy(), track(A).prob.to_numpy(), track(A).prob.to_numpy()]
    return w, build_pass(w, per, 0.95, 0.95, RULES, max_gap_bp=0)


def test_vote_failures_are_kept_aside_not_deleted():
    w, p = two_call_pass()
    _, unvoted = build_region_table_coverage(w, 0.95, 0.95, RULES.coverage_frac, max_gap_bp=0)
    assert len(unvoted) == 2
    assert list(p.calls.start) == [unvoted.start[0]]
    assert list(p.filtered.start) == [unvoted.start[1]] and list(p.filtered["filter"]) == [FILTER_SEED_VOTE]
    assert len(p.per_seed) == 3


def test_finalize_records_why_and_measures_filtered_like_calls():
    _, p = two_call_pass()
    per_seed = [[s] for s in p.per_seed]
    # B is also mostly N: it fails both filters; A is clean
    n_frac = lambda df: np.where(df.start.to_numpy() >= 100_000, 0.8, 0.0)  # noqa: E731
    out = finalize(p.calls, p.filtered, per_seed, {"c": CONTIG_LEN}, n_frac, RULES, 0.5)
    assert len(out.calls) == 1 and out.calls.seed_support.tolist() == [3]
    assert out.filtered["filter"].tolist() == [f"{FILTER_SEED_VOTE};{FILTER_N_GAP}"]
    assert out.filtered.seed_support.tolist() == [1]
    assert {"n_frac", "contig_edge"} <= set(out.filtered.columns)
    assert out.meta["seed_support"]["calls_dropped_by_vote"] == 1
    assert out.meta["n_gap_filter"]["calls_dropped"] == 0   # B never passed, so the N filter dropped nothing


def test_n_filter_moves_a_passing_call_to_filtered():
    _, p = two_call_pass()
    out = finalize(p.calls, p.filtered, [[s] for s in p.per_seed], {"c": CONTIG_LEN},
                   lambda df: np.full(len(df), 0.6), RULES, 0.5)
    assert out.calls.empty
    assert sorted(out.filtered["filter"]) == sorted([FILTER_N_GAP, f"{FILTER_SEED_VOTE};{FILTER_N_GAP}"])
    assert out.meta["n_gap_filter"]["calls_dropped"] == 1


def test_split_pieces_report_the_support_their_parent_call_passed():
    # one coarse call over 40-80 kb, made by all three checkpoints; two of them only
    # cover its left half, so the right piece alone overlaps one checkpoint's call
    seeds = [pd.DataFrame({"contig_id": ["c"], "start": [s], "end": [e]})
             for s, e in [(40_000, 80_000), (40_000, 55_000), (40_000, 55_000)]]
    pieces = pd.DataFrame({"genome_id": "g", "contig_id": "c", "start": [41_000, 62_000], "end": [54_000, 79_000],
                           "coarse_start": 40_000, "coarse_end": 80_000, "refine_status": "split", "n_pieces": 2})
    out = finalize(pieces, pieces.iloc[0:0].assign(filter=pd.Series(dtype=str)), [[s] for s in seeds],
                   {"c": CONTIG_LEN}, lambda df: np.zeros(len(df)), RULES, 0.5)
    assert out.calls.seed_support.tolist() == [3, 3]
    assert out.calls.start.tolist() == [41_000, 62_000]   # reported coordinates are the pieces'


def test_snap_applies_to_filtered_candidates():
    w = track((0, 20_000))
    per = [w.prob.to_numpy(), np.full(len(w), 0.01), np.full(len(w), 0.01)]
    p = build_pass(w, per, 0.95, 0.95, RULES, max_gap_bp=0)
    out = finalize(p.calls, p.filtered, [[s] for s in p.per_seed], {"c": CONTIG_LEN},
                   lambda df: np.zeros(len(df)), RULES, 0.5)
    assert out.filtered.start.tolist() == [0] and out.filtered.contig_edge.tolist() == ["start"]


@pytest.mark.parametrize("mode", ["exploratory", "default", "conservative"])
def test_offline_rules_match_the_cli(mode):
    from prophound.cli import build_parser
    from prophound.models.prophage.src.options import PredictOptions
    from prophound.models.prophage.src.predict_genome import call_rules

    args = build_parser().parse_args(["predict", "--input", "x.fna", "--outdir", "o", "--config", "c.yaml",
                                      "--mode", mode, "--checkpoint", "a.pt", "b.pt", "c.pt"])
    opts = PredictOptions.from_args(args, trained_window_len=10_000, trained_stride=1_000)
    assert call_rules(opts) == CallRules.for_mode(mode, window_len=10_000, n_checkpoints=3)


def test_single_checkpoint_never_votes():
    assert CallRules.for_mode("default", 10_000, 1).min_seed_support == 1


def test_concat_nonempty_skips_empty_frames():
    a = pd.DataFrame({"x": pd.Series([], dtype=object)})
    b = pd.DataFrame({"x": [1, 2]})
    assert concat_nonempty([a, b])["x"].dtype == np.int64
    assert concat_nonempty([a, a]) is a
