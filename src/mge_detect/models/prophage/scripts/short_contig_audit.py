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


"""Corpus-wide input-quality audit: short contigs and ambiguity codes (docs/cli-and-packaging.md).

`predict_genome.build_windows` skips any contig shorter than `window_len`
(10 kb) and nothing warns, so an assembly's short contigs are silently never
examined. Verified end to end 2026-09-22 on a mixed assembly: an 8 kb contig
carrying real prophage sequence never reached `windows.tsv` while `regions.tsv`
looked like a complete result.

This job quantifies the blind spot over the whole corpus and, crucially, over
the geNomad truth: how many known elements sit wholly or partly inside contigs
the tool would never look at. Pure CPU, no inference, no GPU.

Reports per genome and in aggregate:
  * contigs and bp below `--window_len`, as a fraction of the assembly
  * positive elements (from `--beds_dir`) inside a skipped contig -- these are
    unreachable recall, not model error
  * the worst-affected genomes, which is where fragmented/MAG assemblies show up

It also counts **IUPAC ambiguity codes** on the same pass, since it is the same
scan and ambiguity is the other half of the input-quality axis (`docs/v1-scope.md`
gating audit A.2). Ambiguity codes collapse a window's score to zero
(`docs/scientific-validity.md` "1a"), so how common they are across the corpus,
and at what per-genome rate, decides whether the encoding fix is a headline
correctness issue or an edge case. Job `816446` could not answer this: only 3 of
its 42 eval genomes carry any codes at all.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import pandas as pd
from pyfaidx import Fasta


def contig_lengths(path: str, count_ambiguity: bool = False):
    """Return ({contig: length}, ambiguity_bases, contigs_with_ambiguity).

    Ambiguity is counted on the same pass -- reading every FASTA twice for two
    cheap statistics would double a corpus-wide scan for no reason.
    """
    fa = Fasta(path, as_raw=True, sequence_always_upper=True)
    try:
        lens, n_amb, amb_contigs = {}, 0, 0
        for k in fa.keys():
            if count_ambiguity:
                seq = str(fa[k][:])
                lens[k] = len(seq)
                c = sum(1 for b in seq if b not in "ACGT")
                n_amb += c
                amb_contigs += 1 if c else 0
            else:
                lens[k] = len(fa[k])
        return lens, n_amb, amb_contigs
    finally:
        fa.close()  # one descriptor per genome otherwise leaks (docs/environment.md)


def read_bed_positives(bed_path: Path):
    out = []
    for line in open(bed_path):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 4 or parts[3] != "positive":
            continue
        out.append((parts[0], int(parts[1]), int(parts[2])))
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sequences_dir", required=True)
    ap.add_argument("--beds_dir", help="optional: beds_pos, to count unreachable truth elements")
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--window_len", type=int, default=10000)
    ap.add_argument("--fasta_suffix", default=".fna")
    ap.add_argument("--limit", type=int, default=0, help="cap genomes scanned (0 = all)")
    ap.add_argument("--no_ambiguity", action="store_true",
                    help="skip the ambiguity count (faster: lengths then come from the .fai)")
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    fastas = sorted(Path(args.sequences_dir).glob(f"*{args.fasta_suffix}"))
    if args.limit:
        fastas = fastas[: args.limit]
    print(f"[scan] {len(fastas)} genomes, window_len={args.window_len}", flush=True)

    rows = []
    for i, fp in enumerate(fastas):
        genome = fp.name[: -len(args.fasta_suffix)]
        try:
            lens, n_amb, amb_contigs = contig_lengths(str(fp), count_ambiguity=not args.no_ambiguity)
        except Exception as e:
            print(f"[skip] {genome}: {e}", flush=True)
            continue
        short = {k: v for k, v in lens.items() if v < args.window_len}
        total_bp = sum(lens.values())
        row = dict(genome_id=genome, n_contigs=len(lens), n_short_contigs=len(short),
                   total_bp=total_bp, short_bp=sum(short.values()),
                   frac_bp_skipped=sum(short.values()) / max(1, total_bp),
                   n_ambiguity=n_amb, ambiguity_ppm=n_amb / max(1, total_bp) * 1e6,
                   contigs_with_ambiguity=amb_contigs)
        if args.beds_dir:
            bed = Path(args.beds_dir) / f"{genome}.bed"
            # an element lives on exactly one contig, so it is either reachable
            # or entirely unreachable -- there is no partial case
            whole = total_el = 0
            if bed.exists():
                for contig, _s, _e in read_bed_positives(bed):
                    total_el += 1
                    if contig in short:
                        whole += 1
            row.update(n_elements=total_el, elements_in_skipped_contigs=whole)
        rows.append(row)
        if (i + 1) % 200 == 0:
            print(f"[scan] {i + 1}/{len(fastas)}", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(out / "per_genome.csv", index=False)
    agg = {
        "n_genomes": int(len(df)),
        "window_len": args.window_len,
        "genomes_with_short_contigs": int((df.n_short_contigs > 0).sum()),
        "total_contigs": int(df.n_contigs.sum()),
        "total_short_contigs": int(df.n_short_contigs.sum()),
        "total_bp": int(df.total_bp.sum()),
        "total_bp_skipped": int(df.short_bp.sum()),
        "frac_bp_skipped_overall": float(df.short_bp.sum() / max(1, df.total_bp.sum())),
        "genomes_losing_over_1pct": int((df.frac_bp_skipped > 0.01).sum()),
        "genomes_losing_over_10pct": int((df.frac_bp_skipped > 0.10).sum()),
    }
    if "elements_in_skipped_contigs" in df:
        agg["total_elements"] = int(df.n_elements.sum())
        agg["elements_in_skipped_contigs"] = int(df.elements_in_skipped_contigs.sum())
        agg["frac_elements_unreachable"] = float(
            df.elements_in_skipped_contigs.sum() / max(1, df.n_elements.sum()))
    if "n_ambiguity" in df and (df.n_ambiguity > 0).any():
        aff = df[df.n_ambiguity > 0]
        agg["ambiguity"] = {
            "genomes_with_any": int(len(aff)),
            "frac_genomes_with_any": float((df.n_ambiguity > 0).mean()),
            "total_ambiguity_bases": int(df.n_ambiguity.sum()),
            "overall_ppm": float(df.n_ambiguity.sum() / max(1, df.total_bp.sum()) * 1e6),
            "median_ppm_among_affected": float(aff.ambiguity_ppm.median()),
            "max_ppm": float(df.ambiguity_ppm.max()),
            "genomes_over_100ppm": int((df.ambiguity_ppm > 100).sum()),
            "genomes_over_1000ppm": int((df.ambiguity_ppm > 1000).sum()),
            "note": ("A window is ~10 kb, so even 100 ppm puts ~1 ambiguous base in a typical "
                     "window, and two contiguous codes take a confident window from 1.000 to "
                     "0.000. genomes_over_100ppm is the count at material risk under the "
                     "n-channel default. EDL933, the known failure, is 1,182 ppm."),
        }
    agg["interpretation"] = (
        "frac_elements_unreachable is recall the tool cannot achieve at any threshold "
        "or operating point, because the sequence is never examined. If it is non-trivial "
        "the short-contig path needs fixing before any recall figure is quoted for "
        "fragmented assemblies; if it is ~0 on this corpus, the risk is specific to MAG "
        "and draft input rather than to the benchmark, and a warning suffices."
    )
    (out / "results.json").write_text(json.dumps(agg, indent=2) + "\n")
    print(json.dumps(agg, indent=2), flush=True)
    worst = df.sort_values("frac_bp_skipped", ascending=False).head(10)
    print("\nworst-affected genomes:\n" + worst.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
