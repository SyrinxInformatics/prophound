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
"""C5 of docs/research-integrity.md: how much of each evaluation element is in training.

For every element of every --set, the fraction of its canonical 21-mers that
occur in any training-split positive element (beds_pos intervals of train
genomes). This is the measure small_element_spikein.py's hold-out check uses
(`disjoint`), applied to the evaluation sets instead of the implants. Nucleotide
level only: it cannot see protein-level homology (F5).

A --set is LABEL=PATH. The file needs genome, contig, start and end columns:
`accession` or `genome_id` for the genome, and 1-based inclusive coordinates
(the eval manifests, literature_truth.tsv and insitu_truth.tsv all are). Each
set's genomes are read from --genomes_dir, or --spiked_dir for `spk_` genomes.

Example:
  python scripts/eval_containment.py --beds_dir beds_pos --sequences_dir sequences \\
      --train_split splits/train.txt --genomes_dir ../raw/genomes \\
      --set sizestrat42=scripts/eval_sizestrat_elements.csv --out analysis/integrity/containment.tsv
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from small_element_spikein import canon_kmers  # noqa: E402


def read_set(label: str, path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t" if path.endswith((".tsv", ".txt")) else ",")
    df = df.rename(columns={"genome_id": "accession"})
    need = {"accession", "contig", "start", "end"}
    if not need <= set(df.columns):
        raise SystemExit(f"{label}: {path} lacks columns {sorted(need - set(df.columns))}")
    keep = [c for c in ("name", "stratum", "class", "genomad_label", "split") if c in df.columns]
    return df[["accession", "contig", "start", "end", *keep]].assign(set=label)


def main() -> None:
    from pyfaidx import Fasta
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--beds_dir", required=True)
    ap.add_argument("--sequences_dir", required=True)
    ap.add_argument("--train_split", required=True)
    ap.add_argument("--genomes_dir", required=True)
    ap.add_argument("--spiked_dir", default=None)
    ap.add_argument("--set", action="append", required=True, metavar="LABEL=PATH")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    ev = pd.concat([read_set(*s.split("=", 1)) for s in a.set], ignore_index=True)
    codes, missing = [], []
    for g, sub in ev.groupby("accession", sort=False):
        d = Path(a.spiked_dir) if (a.spiked_dir and g.startswith("spk_")) else Path(a.genomes_dir)
        fp = d / f"{g}.fna"
        if not fp.exists():
            missing.append(g)
            codes += [(i, np.empty(0, np.uint64)) for i in sub.index]
            continue
        fa = Fasta(str(fp), as_raw=True, sequence_always_upper=True)
        try:
            for i, r in sub.iterrows():
                codes.append((i, canon_kmers(fa[r.contig][int(r.start) - 1:int(r.end)])))
        finally:
            fa.close()
    if missing:
        raise SystemExit(f"{len(missing)} genomes have no FASTA, e.g. {missing[:3]}")
    codes = dict(codes)
    owner = np.concatenate([np.full(codes[i].size, i) for i in ev.index])
    allc = np.concatenate([codes[i] for i in ev.index])
    uniq, inv = np.unique(allc, return_inverse=True)    # sort once; each genome is a binary search
    seen_u = np.zeros(uniq.size, bool)

    train = set(open(a.train_split).read().split())
    n_iv = 0
    for bed in sorted(Path(a.beds_dir).glob("*")):
        g = re.sub(r"\.(bed|tsv|txt)$", "", bed.name)
        if g not in train or not (Path(a.sequences_dir) / f"{g}.fna").exists():
            continue
        ivs = [(p[0], int(p[1]), int(p[2])) for p in (l.split() for l in open(bed) if l.strip() and not l.startswith("#"))
               if len(p) >= 4 and p[3] == "positive"]
        if not ivs:
            continue
        fa = Fasta(str(Path(a.sequences_dir) / f"{g}.fna"), as_raw=True, sequence_always_upper=True)
        try:
            km = np.unique(np.concatenate([canon_kmers(fa[c][max(0, s):e]) for c, s, e in ivs if c in fa.keys()]
                                          or [np.empty(0, np.uint64)]))
        finally:
            fa.close()
        pos = np.searchsorted(uniq, km)
        ok = pos < uniq.size
        ok[ok] = uniq[pos[ok]] == km[ok]
        seen_u[pos[ok]] = True
        n_iv += len(ivs)
    seen = seen_u[inv]
    hit = np.bincount(owner[seen], minlength=int(ev.index.max()) + 1)
    tot = np.bincount(owner, minlength=int(ev.index.max()) + 1)
    ev["n_kmers"] = tot[ev.index]
    ev["cont_train_positive"] = np.round(hit[ev.index] / np.maximum(1, tot[ev.index]), 4)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    ev.to_csv(a.out, sep="\t", index=False)
    print(f"[containment] {n_iv} training positive intervals; {len(ev)} evaluation elements", flush=True)
    for s, sub in ev.groupby("set"):
        q = sub.cont_train_positive
        print(f"[{s}] n={len(sub)} median {q.median():.3f}; >=0.1: {(q >= 0.1).sum()}, >=0.5: {(q >= 0.5).sum()}, "
              f">=0.9: {(q >= 0.9).sum()}", flush=True)


if __name__ == "__main__":
    main()
