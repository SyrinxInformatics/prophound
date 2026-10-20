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
"""The input layer: every supported format gives the same normalised FASTA, and
unusable inputs fail with a reason."""

import bz2
import gzip
import lzma
import os

import pytest
from pyfaidx import Fasta

from prophound.models.prophage.src.inputs import genome_id_for, prepare_input
from prophound.models.prophage.src.io_checks import InputError, validate_fasta

SEQS = {"chr1": "ACGTTGCA" * 40 + "NNRY", "plasmid1": "GGGCCCAT" * 25}


def fasta_text(width=60, newline="\n", lower=False):
    out = []
    for name, seq in SEQS.items():
        seq = seq.lower() if lower else seq
        out.append(f">{name} some description")
        out += [seq[i:i + width] for i in range(0, len(seq), width)]
    return newline.join(out) + newline


def records(prepared):
    fa = Fasta(prepared.fasta, as_raw=True)
    try:
        return {k: fa[k][:] for k in fa.keys()}  # noqa: SIM118 -- pyfaidx Fasta is not a dict
    finally:
        fa.close()


@pytest.fixture
def work(tmp_path):
    d = tmp_path / "work"
    d.mkdir()
    return d


def write(tmp_path, name, data):
    p = tmp_path / name
    p.write_bytes(data if isinstance(data, bytes) else data.encode())
    return str(p)


@pytest.mark.parametrize("variant", ["plain", "crlf", "lower", "uneven"])
def test_fasta_variants_normalise_to_the_same_sequence(tmp_path, work, variant):
    text = {"plain": fasta_text(), "crlf": fasta_text(newline="\r\n"), "lower": fasta_text(lower=True),
            # uneven line widths crashed pyfaidx before (`cat a.fna b.fna` from two sources)
            "uneven": ">chr1 x\n" + SEQS["chr1"][:70] + "\n" + SEQS["chr1"][70:] + "\n>plasmid1\n"
                      + SEQS["plasmid1"] + "\n"}[variant]
    p = prepare_input(write(tmp_path, "g.fna", text), work)
    assert records(p) == SEQS and p.format == "fasta" and p.compression is None
    with open(p.fasta) as fh:
        assert all(len(line.rstrip("\n")) <= 80 for line in fh)
    validate_fasta(p.fasta)                                     # passes on the normalised copy


@pytest.mark.parametrize("name, compress", [("g.fna.gz", gzip.compress), ("g.fna.bz2", bz2.compress),
                                            ("g.fna.xz", lzma.compress)])
def test_compressed_fasta(tmp_path, work, name, compress):
    p = prepare_input(write(tmp_path, name, compress(fasta_text().encode())), work)
    assert records(p) == SEQS and p.compression in ("gzip", "bzip2", "xz") and p.genome_id == "g"


def flatfile(tmp_path, fmt, compress=None):
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    recs = [SeqRecord(Seq(s), id=f"NC_00000{i}.1", name=f"LOCUS{i}", description=name,
                      annotations={"molecule_type": "DNA"}) for i, (name, s) in enumerate(SEQS.items(), 1)]
    p = tmp_path / f"g.{fmt}"
    SeqIO.write(recs, p, fmt)
    if compress:
        p.write_bytes(compress(p.read_bytes()))
    return str(p)


@pytest.mark.parametrize("fmt", ["genbank", "embl"])
def test_flat_files_use_the_versioned_accession(tmp_path, work, fmt):
    p = prepare_input(flatfile(tmp_path, fmt, compress=gzip.compress), work)
    assert p.format == fmt and p.compression == "gzip"
    got = records(p)
    assert sorted(got) == ["NC_000001.1", "NC_000002.1"]
    assert list(got.values()) == list(SEQS.values())
    if fmt == "genbank":
        assert p.contig_names == {"NC_000001.1": "LOCUS1", "NC_000002.1": "LOCUS2"}


def test_genbank_without_sequence_is_refused(tmp_path, work):
    text = ("LOCUS       NZ_X01000000             0 bp    DNA     linear   CON 01-JAN-2020\n"
            "DEFINITION  master record.\nACCESSION   NZ_X01000000\nVERSION     NZ_X01000000.1\n"
            "FEATURES             Location/Qualifiers\nCONTIG      join(NZ_X01000001.1:1..100)\n//\n")
    with pytest.raises(InputError, match="no sequence"):
        prepare_input(write(tmp_path, "master.gbff", text), work)


def test_gff3_with_embedded_fasta(tmp_path, work):
    text = "##gff-version 3\nchr1\tProkka\tCDS\t1\t30\t.\t+\t0\tID=x\n##FASTA\n" + fasta_text()
    p = prepare_input(write(tmp_path, "g.gff", text), work)
    assert p.format == "gff3" and records(p) == SEQS


@pytest.mark.parametrize("name, data, message", [
    ("g.gff", "##gff-version 3\nchr1\tx\tCDS\t1\t3\t.\t+\t0\tID=x\n", "without a ##FASTA section"),
    ("r.fastq", "@read1\nACGT\n+\nIIII\n", "FASTQ"),
    ("a.sam", "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n", "SAM"),
    ("a.bam", gzip.compress(b"BAM\x01rest"), "BAM"),
    ("x.txt", "hello world\n", "unrecognised format"),
    ("empty.fna", "", "empty"),
])
def test_unusable_inputs_are_refused(tmp_path, work, name, data, message):
    with pytest.raises(InputError, match=message):
        prepare_input(write(tmp_path, name, data), work)


def test_protein_is_refused_by_the_validator(tmp_path, work):
    p = prepare_input(write(tmp_path, "p.faa", ">p1\n" + "MKVLHEQRSTWFPLIVDE" * 10 + "\n"), work)
    with pytest.raises(InputError, match="protein"):
        validate_fasta(p.fasta, label="p.faa")


def test_rna_u_is_read_as_t(tmp_path, work):
    # before, impute replaced every U with a random base, silently
    p = prepare_input(write(tmp_path, "r.fna", ">r\nACGUACGUuu\n"), work)
    assert records(p) == {"r": "ACGTACGTTT"} and p.u_as_t == 4


def test_read_only_input_directory_and_no_index_left_behind(tmp_path, work):
    src = tmp_path / "ro"
    src.mkdir()
    path = write(src, "g.fna", fasta_text())
    os.chmod(src, 0o555)
    try:
        p = prepare_input(path, work)
        Fasta(p.fasta).close()                                   # builds the .fai in the work dir
        assert sorted(os.listdir(src)) == ["g.fna"]
    finally:
        os.chmod(src, 0o755)


@pytest.mark.parametrize("name, expected", [("GCF_1.2_ASM_genomic.fna", "GCF_1.2_ASM_genomic"),
                                            ("x.fna.gz", "x"), ("x.gbff.bz2", "x"), ("x.gff3", "x")])
def test_genome_id(name, expected):
    assert genome_id_for(name) == expected
