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
"""Batch input: which files a batch scans, what it refuses up front, and its combined outputs."""

import json

import pandas as pd
import pytest

from prophound.models.prophage.src import predict_genome as pg
from prophound.models.prophage.src.batch import (
    GenomeRow,
    batch_outputs_held,
    collect_inputs,
    summary_table,
    unique_names,
)
from prophound.models.prophage.src.options import PredictError


def touch(path, text=">c\nACGT\n"):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return str(path)


def test_single_file_is_not_a_batch(tmp_path):
    f = touch(tmp_path / "g.fna")
    assert collect_inputs([f], None) == ([f], False)


def test_directory_scan_takes_genome_files_only(tmp_path):
    d = tmp_path / "genomes"
    for name in ("b.fna", "a.gbff.gz", "c.fa.xz", "notes.txt", "b.fna.fai", ".hidden.fna"):
        touch(d / name)
    paths, batch = collect_inputs([str(d)], None)
    assert batch and [p.split("/")[-1] for p in paths] == ["a.gbff.gz", "b.fna", "c.fa.xz"]


def test_input_list_resolves_relative_paths_and_skips_comments(tmp_path):
    touch(tmp_path / "data" / "x.fna")
    lst = tmp_path / "data" / "list.txt"
    lst.write_text("# genomes\nx.fna  # first\n\n")
    paths, batch = collect_inputs(None, str(lst))
    assert batch and paths == [str(tmp_path / "data" / "x.fna")]


@pytest.mark.parametrize("make, message", [
    (lambda t: ([touch(t / "a" / "g.fna"), touch(t / "b" / "g.fa")], None), "share a genome ID"),
    (lambda t: ([touch(t / "g.fna"), str(t / "missing.fna")], None), "not found"),
    (lambda t: ([str((t / "empty").mkdir() or t / "empty")], None), "no genome files"),
    (lambda t: (None, None), "no input"),
])
def test_bad_batches_are_refused_before_any_work(tmp_path, make, message):
    with pytest.raises(PredictError, match=message):
        collect_inputs(*make(tmp_path))


def test_summary_counts_stay_integers():
    table = summary_table([GenomeRow("a", "a.fna", "ok", "fasta", 2, 370_000, 5, 1, 3.2),
                           GenomeRow("b", "b.fna", "failed", error="nothing to scan")])
    text = table.to_csv(sep="\t", index=False)
    assert "\t370000\t" in text and ".0\t" not in text.split("\n")[1]


def test_combined_names_are_unique():
    calls = pd.DataFrame({"call_id": ["prophound_1", "prophound_1"], "genome_id": ["a", "b"]})
    assert list(unique_names(calls).call_id) == ["a:prophound_1", "b:prophound_1"]


def test_batch_outdir_check_sees_genome_directories(tmp_path):
    touch(tmp_path / "g1" / "regions.tsv", "")
    assert batch_outputs_held(tmp_path, ["g1", "g2"], pg._held_outputs) == ["g1/"]
    touch(tmp_path / "batch_summary.tsv", "")
    assert batch_outputs_held(tmp_path, ["g1"], pg._held_outputs) == ["batch_summary.tsv", "g1/"]


def test_force_clears_genome_directories_the_new_batch_does_not_rewrite(tmp_path):
    (tmp_path / "run_metadata.json").write_text(json.dumps(
        {"batch": True, "genomes": [{"genome_id": "old"}, {"genome_id": "kept"}]}))
    touch(tmp_path / "old" / "regions.tsv", "")
    touch(tmp_path / "kept" / "regions.tsv", "")
    touch(tmp_path / "old" / "my_notes.txt", "mine")       # not ours: stays, and so does its directory
    stage = pg.StagedOutputs(tmp_path, force=True, keep_dirs=("kept",))
    stage.path("regions.tsv").write_text("new")
    stage.commit()
    assert not (tmp_path / "old" / "regions.tsv").exists() and (tmp_path / "old" / "my_notes.txt").exists()
    assert (tmp_path / "kept" / "regions.tsv").exists()
