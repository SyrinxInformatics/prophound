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
"""Short contigs (< --window-len): tiled at the largest scale that fits, never
silently dropped from the report."""

from prophound.models.prophage.src.predict_genome import (
    MODE_THRESHOLDS,
    SHORT_CONTIG_SCALES,
    SHORT_CONTIG_THRESHOLDS,
    build_short_contig_windows,
)


def write_fasta(path, lengths):
    with open(path, "w") as fh:
        for name, n in lengths.items():
            fh.write(f">{name}\n{'ACGT' * (n // 4)}\n")


def test_short_contigs_get_largest_fitting_scale_and_are_reported(tmp_path):
    fa = tmp_path / "g.fna"
    write_fasta(fa, {"long": 12_000, "c7k": 7_000, "c4k": 4_000, "c2k": 2_000})
    w, rep = build_short_contig_windows(str(fa), "g", window_len=10_000)
    assert set(w.contig_id) == {"c7k", "c4k"}                       # long is the coarse pass's
    assert (w[w.contig_id == "c7k"].scale == 5000).all()
    assert (w[w.contig_id == "c4k"].scale == 3000).all()
    assert list(w[w.contig_id == "c7k"].start) == [0, 500, 1000, 1500, 2000]
    assert ((w.end - w.start) == w.scale).all()
    c7 = w[w.contig_id == "c7k"]
    assert (c7.end <= 7_000).all()
    assert rep == {"n_short": 3, "bp_short": 13_000, "n_scannable": 2,
                   "n_too_short": 1, "bp_too_short": 2_000}


def test_no_short_contigs_gives_empty_frame(tmp_path):
    fa = tmp_path / "g.fna"
    write_fasta(fa, {"long": 20_000})
    w, rep = build_short_contig_windows(str(fa), "g", window_len=10_000)
    assert w.empty and rep["n_short"] == 0


def test_short_contig_thresholds_cover_every_mode_and_scale_and_keep_mode_order():
    assert set(SHORT_CONTIG_THRESHOLDS) == set(MODE_THRESHOLDS)
    for scale in SHORT_CONTIG_SCALES:
        by_mode = [SHORT_CONTIG_THRESHOLDS[m][scale] for m in ("exploratory", "default", "conservative")]
        assert by_mode == sorted(by_mode)


def test_seed_support_counts_checkpoints_with_an_overlapping_call():
    import pandas as pd
    from prophound.models.prophage.src.predict_genome import seed_support

    calls = pd.DataFrame({"contig_id": ["a", "a", "b"], "start": [0, 50_000, 0], "end": [20_000, 60_000, 5_000]})
    s1 = pd.DataFrame({"contig_id": ["a", "b"], "start": [10_000, 4_999], "end": [15_000, 9_000]})
    s2 = pd.DataFrame({"contig_id": ["a"], "start": [19_999], "end": [30_000]})
    s3 = pd.DataFrame({"contig_id": ["a"], "start": [20_000], "end": [49_999]})   # touches neither a-call
    assert seed_support(calls, [s1, s2, s3]) == [2, 0, 1]


def test_apply_seed_vote_keeps_calls_with_enough_support():
    import pandas as pd
    from prophound.models.prophage.src.predict_genome import apply_seed_vote

    calls = pd.DataFrame({"contig_id": ["a", "a", "b"], "start": [0, 50_000, 0], "end": [20_000, 60_000, 5_000]})
    s1 = pd.DataFrame({"contig_id": ["a", "b"], "start": [10_000, 4_999], "end": [15_000, 9_000]})
    s2 = pd.DataFrame({"contig_id": ["a"], "start": [19_999], "end": [30_000]})
    s3 = pd.DataFrame({"contig_id": ["a"], "start": [20_000], "end": [49_999]})
    kept, dropped = apply_seed_vote(calls, [s1, s2, s3], 2)      # supports are [2, 0, 1]
    assert list(kept.start) == [0] and dropped == 2
    same, none = apply_seed_vote(calls, [s1, s2, s3], 1)
    assert same.equals(calls) and none == 0


def test_mode_vote_table_covers_every_mode():
    from prophound.models.prophage.src.predict_genome import (
        MODE_MIN_SEED_SUPPORT, MODE_THRESHOLDS, SHORT_CONTIG_THRESHOLDS)

    assert set(MODE_MIN_SEED_SUPPORT) == set(MODE_THRESHOLDS) == set(SHORT_CONTIG_THRESHOLDS)
    assert MODE_MIN_SEED_SUPPORT["exploratory"] == 1   # discovery mode never drops a call
