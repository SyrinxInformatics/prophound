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
"""Output files: call table layout, the outdir policy, atomic writes, metadata paths."""

from pathlib import Path

import pandas as pd

from prophound.models.prophage.src import predict_genome as pg


def calls():
    return pd.DataFrame({
        "genome_id": "g", "contig_id": ["c2", "c1", "c1"], "start": [500, 9_000, 100], "end": [900, 9_500, 400],
        "mean_prob": [0.123456, 0.5, 0.987654], "max_prob": [1.0, 0.99999, 0.5],
        "boundary_kind": ["coverage_frac_0.1_short5000", "coverage_frac_0.1_gap50000_br0.1", "coverage_frac_0.1"],
    })


def test_call_table_orders_numbers_and_rounds():
    out = pg.call_table(calls(), {"c1": 0, "c2": 1}, "prophound_", 10_000)
    assert list(out.call_id) == ["prophound_1", "prophound_2", "prophound_3"]
    assert list(zip(out.contig_id, out.start)) == [("c1", 100), ("c1", 9_000), ("c2", 500)]
    assert list(out.window_bp) == [10_000, 10_000, 5_000]
    assert "boundary_kind" not in out.columns and out.columns[0] == "call_id"
    assert list(out.mean_prob) == [0.9877, 0.5, 0.1235]


def test_call_table_handles_no_calls():
    out = pg.call_table(calls().iloc[0:0], {}, "filtered_", 10_000)
    assert out.empty and out.columns[0] == "call_id" and "window_bp" in out.columns


def test_home_directory_is_shown_as_tilde():
    home = str(Path.home())
    assert pg._shown(home + "/x/best.pt") == "~/x/best.pt"
    assert pg._shown("/opt/x.pt") == "/opt/x.pt"
    assert pg._command_line(["predict", "--input", home + "/a b.fna"]) == "prophound predict --input ~/'a b.fna'"


def test_outdir_with_outputs_is_refused_and_next_free_suggested(tmp_path):
    out = tmp_path / "out"
    out.mkdir()
    (out / "notes.txt").write_text("not ours")
    assert pg._held_outputs(out) == []                       # unrelated files are fine
    (out / "regions.tsv").write_text("")
    (tmp_path / "out_1").mkdir()
    (tmp_path / "out_1" / "windows.tsv").write_text("")
    assert pg._held_outputs(out) == ["regions.tsv"]
    assert pg._next_free(out) == tmp_path / "out_2"


def test_staged_outputs_commit_and_discard(tmp_path):
    (tmp_path / "regions.fna").write_text("stale")
    stage = pg.StagedOutputs(tmp_path, force=True)
    stage.path("regions.tsv").write_text("new")
    assert not (tmp_path / "regions.tsv").exists()          # nothing visible before commit
    stage.commit()
    assert (tmp_path / "regions.tsv").read_text() == "new"
    assert not (tmp_path / "regions.fna").exists()           # --force removed the stale optional file
    assert not list(tmp_path.glob(".*prophound-tmp"))

    stage = pg.StagedOutputs(tmp_path, force=False)
    stage.path("windows.tsv").write_text("partial")
    stage.discard()
    assert not (tmp_path / "windows.tsv").exists() and not list(tmp_path.glob(".*prophound-tmp"))
