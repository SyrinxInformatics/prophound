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
"""Sequence-only signals that an input genome is not prokaryotic (no genes, no database).

prophound is trained on bacterial genomes; on eukaryotic sequence it calls large blocks
(analysis/robustness/). These signals feed a command-line warning (docs/cli-and-packaging.md):

- long-ORF coverage: the share of bases inside a stop-to-stop open reading frame of at least
  MIN_ORF_CODONS codons in any of six frames. Prokaryotic genomes are densely coding; most
  eukaryotic genomes are not. `orf_excess` subtracts the same share on a mononucleotide-shuffled
  copy, so a high-GC genome's chance ORFs do not count as coding.
- telomeres: the share of contig ends (contigs of at least TELOMERE_MIN_CONTIG bp) carrying a known
  eukaryotic telomere repeat in tandem within TELOMERE_WINDOW bp of the end.
- total length, GC, CpG observed/expected and the soft-masked (lowercase) share.

The shipped warning uses two of them, chosen on a labeled panel (analysis/robustness/domain_check.tsv): telomere repeats at
MIN_TELOMERE_ENDS or more contig ends, or more than MAX_PROKARYOTE_BP in total. On 166 bacteria and
5 archaea neither fired; on 12 eukaryotes one or the other fired on 8.
"""

from __future__ import annotations

import re

import numpy as np

MIN_ORF_CODONS = 100
TELOMERE_WINDOW = 500
TELOMERE_MIN_CONTIG = 10_000
SHUFFLE_SEED = 0
# the warning rule: either signal flags the input
MIN_TELOMERE_ENDS = 2
MAX_PROKARYOTE_BP = 17_000_000  # the largest known prokaryotic genome is about 16 Mb

_CODE = np.full(256, 4, dtype=np.uint8)
for _i, _b in enumerate(b"ACGT"):
    _CODE[_b] = _i
    _CODE[ord(chr(_b).lower())] = _i
_COMPLEMENT = np.array([3, 2, 1, 0, 4], dtype=np.uint8)
# forward-strand stop codons TAA, TAG, TGA as 16*b0 + 4*b1 + b2 with A0 C1 G2 T3
_STOPS = np.array([48, 50, 56])

_TELOMERES = re.compile(
    "|".join([
        r"(?:TTAGGG){4,}", r"(?:CCCTAA){4,}",          # vertebrates, many fungi and protists
        r"(?:TTTAGGG){4,}", r"(?:CCCTAAA){4,}",        # plants, green algae
        r"(?:TTAGG){5,}", r"(?:CCTAA){5,}",            # insects and other arthropods
        r"(?:TG{1,3}){10,}", r"(?:C{1,3}A){10,}",      # budding yeasts (irregular)
        r"(?:TTACA?G{2,8}){4,}", r"(?:C{2,8}T?GTAA){4,}",  # fission yeast
        r"(?:TT[TC]AGGG){4,}", r"(?:CCCT[AG]AA){4,}",  # Plasmodium
    ])
)


def telomere_end(end: str) -> bool:
    """Whether a contig end (its first or last TELOMERE_WINDOW bp, any case) carries a telomere repeat."""
    return bool(_TELOMERES.search(end.upper()))


def _orf_covered(codes: np.ndarray) -> np.ndarray:
    """Boolean mask over `codes` (one strand, 0-3 bases, 4 = other) of bases in a long stop-to-stop ORF."""
    n = len(codes)
    covered = np.zeros(n + 1, dtype=np.int32)
    for frame in range(3):
        m = (n - frame) // 3
        if m <= MIN_ORF_CODONS:
            continue
        c = codes[frame:frame + 3 * m].reshape(m, 3).astype(np.int32)
        codon = 16 * c[:, 0] + 4 * c[:, 1] + c[:, 2]
        breaks = np.isin(codon, _STOPS) | (c.max(axis=1) > 3)
        idx = np.concatenate(([-1], np.flatnonzero(breaks), [m]))
        lengths = np.diff(idx) - 1
        keep = lengths >= MIN_ORF_CODONS
        starts = frame + 3 * (idx[:-1][keep] + 1)
        ends = frame + 3 * idx[1:][keep]
        np.add.at(covered, starts, 1)
        np.add.at(covered, ends, -1)
    return np.cumsum(covered[:-1]) > 0


def orf_coverage(codes: np.ndarray) -> float:
    """Share of non-ambiguous bases inside a long ORF on either strand."""
    valid = codes < 4
    if not valid.any():
        return 0.0
    fwd = _orf_covered(codes)
    rev = _orf_covered(_COMPLEMENT[codes[::-1]])[::-1]
    return float(((fwd | rev) & valid).sum() / valid.sum())


def _shuffled(codes: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    out = codes.copy()
    valid = np.flatnonzero(codes < 4)
    out[valid] = rng.permutation(codes[valid])
    return out


def genome_signals(seqs: dict[str, str]) -> dict:
    """The domain signals for one input genome, given its contigs (id -> sequence)."""
    rng = np.random.default_rng(SHUFFLE_SEED)
    total = valid = gc = cg = c_count = g_count = lower = 0
    orf_bp = shuf_bp = 0.0
    ends = tel_ends = 0
    for seq in seqs.values():
        raw = seq.encode("ascii", "replace")
        arr = np.frombuffer(raw, dtype=np.uint8)
        codes = _CODE[arr]
        ok = codes < 4
        n_ok = int(ok.sum())
        total += len(arr)
        valid += n_ok
        lower += int(((arr >= 97) & (arr <= 122)).sum())
        c_n, g_n = int((codes == 1).sum()), int((codes == 2).sum())
        c_count += c_n
        g_count += g_n
        gc += c_n + g_n
        cg += int(((codes[:-1] == 1) & (codes[1:] == 2)).sum())
        if n_ok:
            orf_bp += orf_coverage(codes) * n_ok
            shuf_bp += orf_coverage(_shuffled(codes, rng)) * n_ok
        if len(arr) >= TELOMERE_MIN_CONTIG:
            for end in (seq[:TELOMERE_WINDOW], seq[-TELOMERE_WINDOW:]):
                ends += 1
                tel_ends += telomere_end(end)
    orf = orf_bp / valid if valid else 0.0
    shuf = shuf_bp / valid if valid else 0.0
    cpg_oe = (cg * valid / (c_count * g_count)) if c_count and g_count else 0.0
    return {
        "total_bp": total,
        "n_contigs": len(seqs),
        "gc": gc / valid if valid else 0.0,
        "cpg_oe": cpg_oe,
        "lower_frac": lower / total if total else 0.0,
        "orf_cov": orf,
        "orf_cov_shuffled": shuf,
        "orf_excess": orf - shuf,
        "contig_ends": ends,
        "telomere_ends": tel_ends,
        "telomere_frac": tel_ends / ends if ends else 0.0,
    }


def non_prokaryote_reasons(contigs) -> list[str]:
    """Why an input looks non-prokaryotic, or [] (the warning rule).

    `contigs` yields (length, first TELOMERE_WINDOW bp, last TELOMERE_WINDOW bp) per contig, so a
    caller reads only the contig ends.
    """
    total = ends = tel_ends = 0
    for length, head, tail in contigs:
        total += length
        if length >= TELOMERE_MIN_CONTIG:
            ends += 2
            tel_ends += telomere_end(head) + telomere_end(tail)
    reasons = []
    if tel_ends >= MIN_TELOMERE_ENDS:
        reasons.append(f"eukaryotic telomere repeats at {tel_ends} of {ends} contig ends")
    if total > MAX_PROKARYOTE_BP:
        reasons.append(f"{total / 1e6:.1f} Mb in total, more than any known prokaryotic genome (about 16 Mb)")
    return reasons
