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
"""G5b inputs: the fresh Philympics genomes as FASTA plus a truth table in LAMBDA's layout.

docs/claims-and-evidence.md "G5b preregistration" fixes the genome set and the
element filter; this only converts. Philympics (linsalrob/ProphagePredictionComparisons,
commit pinned in the prereg) ships GenBank files through Git LFS and curated
coordinates in genbank/prophageCoords.tsv (1-based, inclusive).

Writes, so that g5_score.py runs on G5b unchanged:
  fasta/<genome>.fna   every record of the GenBank file (chromosome, plasmids,
                       draft contigs), headers = bare accession without version
  g5b_truth.tsv        assembly, sequence, start, end (0-based half-open); elements
                       with >= 5 genes or >= 5 kb; `<`/`>` fuzzy bounds taken at face value
  g5b_genomes.tsv      assembly, sequence, organism, genome_length, corpus_status, genus
  genome_list.txt

Example:
  python scripts/g5b_prepare.py --genbank_dir analysis/g5b/genbank \\
      --coords analysis/g5b/prophageCoords.tsv --out_dir analysis/g5b
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import pandas as pd
from Bio import SeqIO

MIN_GENES, MIN_BP = 5, 5000


def bare(acc: str) -> str:
    return str(acc).split(".")[0]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genbank_dir", required=True)
    ap.add_argument("--coords", required=True)
    ap.add_argument("--out_dir", required=True)
    a = ap.parse_args()
    out = Path(a.out_dir)
    (out / "fasta").mkdir(parents=True, exist_ok=True)
    coords = pd.read_csv(a.coords, sep="\t")
    genomes, truth = [], []
    for gb in sorted(Path(a.genbank_dir).glob("*.gb.gz")):
        name = gb.name[: -len(".gb.gz")]
        seen, total, org = set(), 0, ""
        with gzip.open(gb, "rt") as fh, open(out / "fasta" / f"{name}.fna", "w") as fo:
            for rec in SeqIO.parse(fh, "genbank"):
                acc = bare(rec.id)
                if acc in seen:
                    raise SystemExit(f"{gb.name}: accession {acc} repeats once the version is dropped")
                seen.add(acc)
                org = org or rec.annotations.get("organism", "")
                fo.write(f">{acc}\n{str(rec.seq).upper()}\n")
                total += len(rec.seq)
        c = coords[gb.name == coords.FILE]
        missing = set(c.ACCESSION.map(bare)) - seen
        if missing:
            raise SystemExit(f"{gb.name}: truth on accessions absent from the GenBank file: {sorted(missing)}")
        for r in c.itertuples():
            # PROPHAGE START / END / GENES IN PROPHAGE; a fuzzy bound (`>1753287`, GenBank style:
            # the element extends past it) is taken at the stated position
            s, e, n = (int(str(x).lstrip("<>")) for x in (r._5, r._6, r._7))
            s, e = min(s, e), max(s, e)
            if n >= MIN_GENES or e - s + 1 >= MIN_BP:
                truth.append((name, bare(r.ACCESSION), s - 1, e))
        genomes.append((name, ";".join(sorted(seen)), org, total, "not_in_corpus", org.split(" ")[0]))
    G = pd.DataFrame(genomes, columns=["assembly", "sequence", "organism", "genome_length", "corpus_status", "genus"])
    T = pd.DataFrame(truth, columns=["assembly", "sequence", "start", "end"]).drop_duplicates()
    G.to_csv(out / "g5b_genomes.tsv", sep="\t", index=False)
    T.to_csv(out / "g5b_truth.tsv", sep="\t", index=False)
    (out / "genome_list.txt").write_text("".join(f"{g}\n" for g in G.assembly))
    print(G[["assembly", "organism", "genome_length"]].assign(records=G.sequence.str.count(";") + 1).to_string(index=False))
    print(f"[done] {len(G)} genomes, {len(T)} truth elements (>= {MIN_GENES} genes or >= {MIN_BP} bp) "
          f"of {len(coords[coords.FILE.isin([g + '.gb.gz' for g in G.assembly])])} rows")


if __name__ == "__main__":
    main()
