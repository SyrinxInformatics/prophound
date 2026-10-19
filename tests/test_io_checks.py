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
"""Input validation and extra output formats (io_checks.py)."""

import gzip

import pandas as pd
import pytest

from prophound.models.prophage.src.io_checks import InputError, validate_fasta, write_bed, write_gff3


def _write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return str(p)


def test_accepts_nucleotide_and_counts_ambiguity(tmp_path):
    p = _write(tmp_path, "g.fna", ">a desc\nACGTNNRY\nACGT\n>b\nGGGG\n")
    s = validate_fasta(p)
    assert s["n_contigs"] == 2 and s["total_bp"] == 16
    assert s["ambiguous_bases"] == 4 and s["ambiguous_by_code"] == {"N": 2, "R": 1, "Y": 1}


@pytest.mark.parametrize("text, needle", [
    ("", "empty"),
    ("ACGT\n", "does not look like FASTA"),
    (">a\nACGT\n>a\nACGT\n", "duplicate contig IDs: a"),
    (">a\n>b\nACGT\n", "no sequence: a"),
    (">p\nMKVLAAGIVGLLLAQEPTTRFELIKQEPLQ\n", "looks like protein"),
])
def test_rejects_unusable_input(tmp_path, text, needle):
    p = _write(tmp_path, "x.fna", text)
    with pytest.raises(InputError, match=needle):
        validate_fasta(p)


def test_rejects_missing_and_gzip(tmp_path):
    with pytest.raises(InputError, match="not found"):
        validate_fasta(str(tmp_path / "nope.fna"))
    gz = tmp_path / "g.fna.gz"
    with gzip.open(gz, "wt") as fh:
        fh.write(">a\nACGT\n")
    with pytest.raises(InputError, match="gzip"):
        validate_fasta(str(gz))


def test_bed_is_zero_based_and_gff3_one_based(tmp_path):
    r = pd.DataFrame({"contig_id": ["c1"], "start": [100], "end": [250], "max_prob": [0.9876],
                      "mean_prob": [0.5], "seed_support": [3], "boundary_kind": ["cov"]})
    write_bed(r, tmp_path / "r.bed")
    write_gff3(r, tmp_path / "r.gff3", "default", "0.1.0")
    assert (tmp_path / "r.bed").read_text() == "c1\t100\t250\tprophound_1\t988\t.\n"
    gff = (tmp_path / "r.gff3").read_text().splitlines()
    assert gff[0] == "##gff-version 3"
    f = gff[2].split("\t")
    assert f[2] == "prophage" and (f[3], f[4]) == ("101", "250")
    assert "seed_support=3" in f[8] and "note=candidate" in f[8]


def test_call_n_fraction_counts_gap_bases(tmp_path):
    from prophound.models.prophage.src.predict_genome import MAX_CALL_N_FRAC, call_n_fraction

    fa = tmp_path / "g.fna"
    fa.write_text(">c1\n" + "ACGT" * 250 + "N" * 1000 + "ACGT" * 250 + "\n")      # 1 kb, 1 kb N, 1 kb
    r = pd.DataFrame({"contig_id": ["c1", "c1", "c1"], "start": [0, 500, 1000], "end": [1000, 1500, 2000]})
    f = call_n_fraction(r, str(fa))
    assert list(f) == [0.0, 0.5, 1.0]
    assert list(f >= MAX_CALL_N_FRAC) == [False, True, True]
