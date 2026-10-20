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
"""Read any supported genome file into one plain, indexed FASTA.

The format is detected from the content, not the file name: FASTA, GenBank/RefSeq,
EMBL, or GFF3 with an embedded `##FASTA` section, each optionally gzip, bzip2 or xz
compressed. Reads and alignments are refused with a reason. The copy is written into
a work directory the caller owns and removes, so the user's directory never gets a
`.fai` index and may be read-only. Sequence is upper-cased and wrapped at 80 columns,
and RNA `U` is read as `T`.
"""

from __future__ import annotations

import bz2
import gzip
import hashlib
import lzma
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, Iterable, Iterator, TextIO, Tuple

from .io_checks import InputError

LINE_WIDTH = 80
COMPRESSION = (
    (b"\x1f\x8b", "gzip", gzip.open),
    (b"BZh", "bzip2", bz2.open),
    (b"\xfd7zXZ\x00", "xz", lzma.open),
)
COMPRESSION_SUFFIXES = (".gz", ".bgz", ".bz2", ".xz")
SAM_HEADER_TAGS = ("@HD\t", "@SQ\t", "@RG\t", "@PG\t", "@CO\t")


@dataclass
class PreparedInput:
    """The input as the pipeline reads it."""

    fasta: str                     # the normalised FASTA in the work directory
    source: str                    # the user's file
    genome_id: str
    format: str                    # fasta, genbank, embl or gff3
    compression: str | None
    sha256: str                    # of the user's file, as given
    u_as_t: int = 0                # RNA U bases read as T
    contig_names: Dict[str, str] = field(default_factory=dict)   # contig ID -> record name, where they differ


def genome_id_for(path: str) -> str:
    """The file name without a compression suffix and its last extension (`x.fna.gz` -> `x`)."""
    name = Path(path).name
    for suffix in COMPRESSION_SUFFIXES:
        if name.lower().endswith(suffix):
            name = name[: -len(suffix)]
            break
    return Path(name).stem


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _text_opener(path: Path) -> Tuple[Callable[[], TextIO], str | None]:
    """A function that opens `path` as text, decompressing if needed, and the compression's name."""
    with open(path, "rb") as fh:
        magic = fh.read(6)
    for prefix, name, opener in COMPRESSION:
        if magic.startswith(prefix):
            with opener(path, "rb") as fh:
                try:
                    head = fh.read(4)
                except (OSError, EOFError, lzma.LZMAError) as e:
                    raise InputError(f"{path} looks {name}-compressed but cannot be read: {e}") from None
            if head == b"BAM\x01":
                raise InputError(f"{path} is a BAM alignment file; give the assembled genome instead")
            return (lambda o=opener: o(path, "rt", encoding="latin-1", newline=None)), name
    return (lambda: open(path, encoding="latin-1", newline=None)), None


def _first_line(fh: TextIO) -> str:
    for line in fh:
        if line.strip():
            return line
    return ""


def _detect(first: str, path: Path) -> str:
    if first.startswith(">"):
        return "fasta"
    if first.startswith("LOCUS"):
        return "genbank"
    if first.startswith("ID   "):
        return "embl"
    if first.startswith("##gff-version"):
        return "gff3"
    if first.startswith(SAM_HEADER_TAGS):
        raise InputError(f"{path} is a SAM alignment file; give the assembled genome instead")
    if first.startswith("@"):
        raise InputError(f"{path} looks like FASTQ sequencing reads; prophound scans assembled contigs, "
                         "so assemble the reads first")
    raise InputError(f"{path}: unrecognised format (expected FASTA, GenBank, EMBL or GFF3 with sequence); "
                     f"first line: {first.strip()[:60]!r}")


def _fasta_records(lines: Iterable[str], path: Path) -> Iterator[Tuple[str, Iterator[str]]]:
    """(header, sequence lines) per record from FASTA text."""
    it = iter(lines)
    header = None
    for line in it:
        if line.startswith(">"):
            header = line.rstrip()
            break
        if line.strip():
            raise InputError(f"{path} does not look like FASTA (first line does not start with '>')")
    while header is not None:
        nxt: list = []

        def body(nxt=nxt):
            for line in it:
                if line.startswith(">"):
                    nxt.append(line.rstrip())
                    return
                yield line
        yield header, body()
        header = nxt[0] if nxt else None


class _Writer:
    """Writes records as upper-case FASTA wrapped at LINE_WIDTH, reading U as T."""

    def __init__(self, fh: TextIO):
        self.fh, self.u_as_t = fh, 0

    def record(self, header: str, chunks: Iterable[str]) -> None:
        self.fh.write(header + "\n")
        buf = ""
        for chunk in chunks:
            seq = "".join(chunk.split()).upper()
            if "U" in seq:
                self.u_as_t += seq.count("U")
                seq = seq.replace("U", "T")
            buf += seq
            while len(buf) >= LINE_WIDTH:
                self.fh.write(buf[:LINE_WIDTH] + "\n")
                buf = buf[LINE_WIDTH:]
        if buf:
            self.fh.write(buf + "\n")


def _write_fasta(fh: TextIO, out: _Writer, path: Path) -> None:
    for header, body in _fasta_records(fh, path):
        out.record(header, body)


def _write_gff3(fh: TextIO, out: _Writer, path: Path) -> None:
    for line in fh:
        if line.startswith("##FASTA"):
            _write_fasta(fh, out, path)
            return
    raise InputError(f"{path} is GFF3 without a ##FASTA section, so it holds no sequence; "
                     "give the genome FASTA (or a GFF3 with the sequence embedded, as Prokka and Bakta write)")


def _write_flatfile(fh: TextIO, out: _Writer, path: Path, fmt: str, names: Dict[str, str]) -> None:
    from Bio import BiopythonParserWarning, SeqIO
    from Bio.Seq import UndefinedSequenceError

    try:
        with warnings.catch_warnings():
            # header-formatting complaints; the sequence itself is validated afterwards
            warnings.simplefilter("ignore", BiopythonParserWarning)
            for rec in SeqIO.parse(fh, fmt):
                try:
                    seq = str(rec.seq)
                except UndefinedSequenceError:
                    seq = ""
                if not seq:
                    raise InputError(f"{path}: record {rec.id} has no sequence (a master or contig-only record?); "
                                     "download the version with sequence, or the FASTA")
                if rec.name and rec.name != rec.id:
                    names[rec.id] = rec.name
                out.record(f">{rec.id} {rec.description}".rstrip(), [seq])
    except InputError:
        raise
    except Exception as e:  # Biopython's parse errors
        raise InputError(f"cannot read {path} as {fmt.upper()}: {e}") from None


def prepare_input(source: str, workdir: Path) -> PreparedInput:
    """Convert `source` to a normalised FASTA in `workdir`; raise InputError on an unusable file."""
    path = Path(source)
    if not path.exists():
        raise InputError(f"input not found: {source}")
    if path.is_dir():
        raise InputError(f"input is a directory: {source}")
    if path.stat().st_size == 0:
        raise InputError(f"input is empty: {source}")
    try:
        open_text, compression = _text_opener(path)
        with open_text() as fh:
            first = _first_line(fh)
    except PermissionError:
        raise InputError(f"cannot read {source}: permission denied") from None
    except (OSError, EOFError, lzma.LZMAError) as e:
        raise InputError(f"cannot read {source}: {e}") from None
    if not first:
        raise InputError(f"input has no content: {source}")
    fmt = _detect(first, path)

    out_path = workdir / "input.fna"
    names: Dict[str, str] = {}
    with open_text() as fh, open(out_path, "w") as out_fh:
        writer = _Writer(out_fh)
        if fmt == "fasta":
            _write_fasta(fh, writer, path)
        elif fmt == "gff3":
            _write_gff3(fh, writer, path)
        else:
            _write_flatfile(fh, writer, path, fmt, names)
    return PreparedInput(fasta=str(out_path), source=source, genome_id=genome_id_for(source), format=fmt,
                         compression=compression, sha256=_sha256(path), u_as_t=writer.u_as_t, contig_names=names)
