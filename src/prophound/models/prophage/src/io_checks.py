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
"""Input validation, provenance and extra output formats for `predict`.

docs/claims-and-evidence.md (CLI usability audit, 2026-09-28): fail early with a
clear message on inputs the pipeline cannot use, record what produced a run, and
write regions in formats other tools read.
"""

from __future__ import annotations

import hashlib
from collections import Counter
from pathlib import Path

import pandas as pd

NUCLEOTIDE = set("ACGTUN")
IUPAC_AMBIGUOUS = set("RYKMSWBDHVN")
# fraction of letters that must be nucleotide or IUPAC codes; protein FASTA is
# mostly other letters (E, F, I, L, P, Q, ...)
MIN_NUCLEOTIDE_FRACTION = 0.9


class InputError(SystemExit):
    """Raised for an unusable input, with a message meant for the user."""

    def __init__(self, msg: str):
        super().__init__(f"error: {msg}")


def package_version() -> str:
    try:
        from importlib.metadata import version
        return version("prophound")
    except Exception:   # not installed (e.g. run from the HPCC source tree)
        return "unknown (not installed)"


def sha256(path: str, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def validate_fasta(path: str) -> dict:
    """Check a nucleotide FASTA in one streaming pass and summarise it.

    Rejects: a missing or empty file, gzip (pyfaidx needs plain text or bgzip,
    and would fail with an obscure error), text that is not FASTA, protein
    FASTA, records with no sequence, and duplicate contig IDs (pyfaidx keys
    on the first word of the header). Returns contig count, total bp and
    per-code counts of ambiguity bases for the run summary.
    """
    p = Path(path)
    if not p.exists():
        raise InputError(f"input FASTA not found: {path}")
    if p.stat().st_size == 0:
        raise InputError(f"input FASTA is empty: {path}")
    with open(p, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        raise InputError(f"{path} is gzip-compressed; decompress it first (e.g. `gunzip -k {p.name}`)")

    ids, dup, empty = set(), [], []
    letters, ambiguous = Counter(), Counter()
    n_contigs, total_bp, cur, cur_len, seen_header = 0, 0, None, 0, False
    with open(p, "r", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line:
                continue
            if line.startswith(">"):
                if cur is not None and cur_len == 0:
                    empty.append(cur)
                cur = line[1:].split()[0] if line[1:].split() else ""
                if not cur:
                    raise InputError(f"{path}: a header line has no sequence ID")
                if cur in ids:
                    dup.append(cur)
                ids.add(cur)
                n_contigs += 1
                cur_len, seen_header = 0, True
                continue
            if not seen_header:
                raise InputError(f"{path} does not look like FASTA (first line does not start with '>')")
            s = line.upper()
            cur_len += len(s)
            total_bp += len(s)
            letters.update(s)
    if cur is not None and cur_len == 0:
        empty.append(cur)
    if n_contigs == 0:
        raise InputError(f"{path} contains no FASTA records")
    if dup:
        raise InputError(f"{path} has duplicate contig IDs: {', '.join(sorted(set(dup))[:5])}"
                         + (" ..." if len(set(dup)) > 5 else ""))
    if empty:
        raise InputError(f"{path} has records with no sequence: {', '.join(empty[:5])}"
                         + (" ..." if len(empty) > 5 else ""))
    alpha = sum(v for k, v in letters.items() if k.isalpha())
    nuc = sum(v for k, v in letters.items() if k in NUCLEOTIDE | IUPAC_AMBIGUOUS)
    if alpha and nuc / alpha < MIN_NUCLEOTIDE_FRACTION:
        raise InputError(f"{path} looks like protein, not nucleotide, sequence "
                         f"({100 * nuc / alpha:.0f}% nucleotide letters); prophound needs nucleotide FASTA")
    for k in IUPAC_AMBIGUOUS:
        if letters.get(k):
            ambiguous[k] = letters[k]
    n_amb = sum(ambiguous.values())
    return {"n_contigs": n_contigs, "total_bp": total_bp, "ambiguous_bases": n_amb,
            "ambiguous_ppm": round(1e6 * n_amb / max(1, total_bp), 2),
            "ambiguous_by_code": dict(sorted(ambiguous.items()))}


def write_bed(regions: pd.DataFrame, path: Path) -> None:
    """BED6: 0-based half-open, as regions.tsv. Score is max_prob scaled to 0-1000."""
    with open(path, "w") as fh:
        for i, r in enumerate(regions.itertuples(index=False), 1):
            score = int(round(1000 * float(getattr(r, "max_prob", 0) or 0)))
            fh.write(f"{r.contig_id}\t{int(r.start)}\t{int(r.end)}\tprophound_{i}\t{score}\t.\n")


def write_gff3(regions: pd.DataFrame, path: Path, mode: str, version: str) -> None:
    """GFF3, 1-based inclusive. Type is the Sequence Ontology `prophage` term;
    calls are candidates, which the `note` attribute says."""
    with open(path, "w") as fh:
        fh.write("##gff-version 3\n")
        fh.write(f"# prophound {version}, mode {mode}. Calls are candidate prophage regions, not annotations.\n")
        for i, r in enumerate(regions.itertuples(index=False), 1):
            attrs = [f"ID=prophound_{i}", "note=candidate"]
            for k in ("max_prob", "mean_prob", "seed_support", "n_frac", "boundary_kind"):
                v = getattr(r, k, None)
                if v is not None and not (isinstance(v, float) and pd.isna(v)):
                    attrs.append(f"{k}={v:.4f}" if isinstance(v, float) else f"{k}={v}")
            score = float(getattr(r, "max_prob", 0) or 0)
            fh.write(f"{r.contig_id}\tprophound\tprophage\t{int(r.start) + 1}\t{int(r.end)}\t{score:.4f}\t.\t.\t"
                     + ";".join(attrs) + "\n")


def write_region_fasta(regions: pd.DataFrame, fasta_path: str, path: Path) -> None:
    from pyfaidx import Fasta
    fa = Fasta(fasta_path, as_raw=True)
    try:
        with open(path, "w") as fh:
            for i, r in enumerate(regions.itertuples(index=False), 1):
                seq = fa[r.contig_id][int(r.start):int(r.end)]
                fh.write(f">prophound_{i} {r.contig_id}:{int(r.start) + 1}-{int(r.end)}\n")
                for j in range(0, len(seq), 80):
                    fh.write(seq[j:j + 80] + "\n")
    finally:
        fa.close()
