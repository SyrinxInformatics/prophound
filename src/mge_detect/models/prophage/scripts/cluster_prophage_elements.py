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


"""Cluster prophage elements by sequence similarity, for clade-held-out splits.

docs/scientific-validity.md item 2a -- the novelty experiment. The project's
splits are by GENOME (make_splits.py), which is honest for "unseen genome" but
does not test the product claim: a test genome's prophage may be a close
homolog of a training prophage, so a model that memorised phage families would
still score well. Holding out whole CLADES is what tests generalization to
divergent elements.

Stage 1 (this script): group the geNomad-called elements into clades, and report
how many genomes a clade-held-out split would actually leave usable. That number
has to be inspected BEFORE any retraining is submitted -- a split that strands
most genomes is not worth 3 GPU-days.

Clustering is MinHash (Mash-style) over canonical k-mers, implemented here in
numpy on purpose: no mmseqs/mash/skani dependency, so it cannot fail on a module
that turns out not to be installed. Sketches are bottom-s hashes of canonical
k-mers; Jaccard between sketches estimates shared k-mer content, and clades are
single-linkage components above --min_jaccard.

Usage (see cluster_prophage_elements.slurm):
    python cluster_prophage_elements.py --beds_dir beds_pos --sequences_dir sequences \
        --out_dir analysis/clades --k 21 --sketch 2000 --min_jaccard 0.2
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

# Module scope on purpose: the HPCC login node cannot load miniconda, so a
# missing library is only discoverable inside a job -- fail in seconds, not
# after the sketch pass (docs/environment.md).
from pyfaidx import Fasta, IndexNotFoundError

BASES = np.full(256, 255, dtype=np.uint8)
for i, b in enumerate(b"ACGT"):
    BASES[b] = i
    BASES[b + 32] = i

MASK64 = np.uint64(0xFFFFFFFFFFFFFFFF)


def _mix64(x: np.ndarray) -> np.ndarray:
    """splitmix64 finalizer -- a cheap, well-mixed 64-bit hash."""
    x = x.astype(np.uint64)
    x = (x ^ (x >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
    x = (x ^ (x >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
    return x ^ (x >> np.uint64(31))


def sketch(seq: str, k: int, s: int) -> np.ndarray:
    """Bottom-s MinHash sketch of the sequence's canonical k-mers."""
    arr = BASES[np.frombuffer(seq.encode("ascii", "ignore"), dtype=np.uint8)]
    arr = arr[arr < 4]
    if arr.size < k:
        return np.empty(0, dtype=np.uint64)
    n = arr.size - k + 1
    fwd = np.zeros(n, dtype=np.uint64)
    rev = np.zeros(n, dtype=np.uint64)
    for i in range(k):
        fwd = (fwd << np.uint64(2)) | arr[i : i + n].astype(np.uint64)
        rev |= (np.uint64(3) - arr[i : i + n].astype(np.uint64)) << np.uint64(2 * i)
    canon = np.minimum(fwd, rev)
    h = _mix64(canon)
    h = np.unique(h)
    return h[:s] if h.size > s else h


def jaccard(a: np.ndarray, b: np.ndarray, s: int) -> float:
    """Bottom-s sketch Jaccard: merge, take the s smallest of the union."""
    if a.size == 0 or b.size == 0:
        return 0.0
    union = np.unique(np.concatenate([a, b]))[:s]
    shared = np.intersect1d(union, np.intersect1d(a, b, assume_unique=True), assume_unique=True).size
    return shared / union.size


def open_fasta(path: str, retries: int = 5, wait_s: float = 2.0) -> Fasta:
    """pyfaidx open that tolerates another process rewriting the .fai.

    pyfaidx rewrites a .fai in place when it is missing or older than its FASTA;
    two jobs opening the same sequences/ on one node can collide on that write
    (composition_baseline job 815767 died this way, 2026-09-21). Retry briefly;
    a persistent failure still raises.
    """
    for attempt in range(retries):
        try:
            return Fasta(path, as_raw=True, sequence_always_upper=True)
        except IndexNotFoundError as e:  # a missing/unreadable FASTA raises differently
            if attempt == retries - 1:
                raise
            print(f"[fasta] {path}: {e}; retry {attempt + 1}/{retries - 1}", flush=True)
            time.sleep(wait_s * (attempt + 1))


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
    ap.add_argument("--beds_dir", required=True)
    ap.add_argument("--sequences_dir", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--k", type=int, default=21)
    ap.add_argument("--sketch", type=int, default=2000)
    ap.add_argument("--min_jaccard", type=float, default=0.2)
    ap.add_argument("--fasta_suffix", default=".fna")
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    beds = sorted(Path(args.beds_dir).glob("*"))
    print(f"[data] {len(beds)} bed files in {args.beds_dir}", flush=True)

    rows, sketches = [], []
    missing_fasta = 0
    for bi, bed in enumerate(beds):
        genome = bed.name
        for ext in (".bed", ".tsv", ".txt"):
            if genome.endswith(ext):
                genome = genome[: -len(ext)]
        fa_path = Path(args.sequences_dir) / f"{genome}{args.fasta_suffix}"
        if not fa_path.exists():
            missing_fasta += 1
            continue
        positives = read_bed_positives(bed)
        if not positives:
            continue
        fa = open_fasta(str(fa_path))
        for contig, start, end in positives:
            try:
                seq = fa[contig][start:end]
            except Exception:
                continue
            sk = sketch(seq, args.k, args.sketch)
            if sk.size == 0:
                continue
            rows.append({"genome_id": genome, "contig_id": contig, "start": start,
                         "end": end, "len_bp": end - start, "n_hashes": int(sk.size)})
            sketches.append(sk)
        if (bi + 1) % 200 == 0:
            print(f"[sketch] {bi + 1}/{len(beds)} genomes, {len(rows)} elements", flush=True)

    print(f"[sketch] {len(rows)} elements sketched ({missing_fasta} genomes had no FASTA)", flush=True)
    if not rows:
        raise SystemExit("No elements found -- check --beds_dir/--sequences_dir and --fasta_suffix")
    df = pd.DataFrame(rows)

    # single-linkage over pairs above the threshold, via union-find
    n = len(sketches)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    print(f"[cluster] {n * (n - 1) // 2} pairs to score", flush=True)
    for i in range(n):
        for j in range(i + 1, n):
            if find(i) != find(j) and jaccard(sketches[i], sketches[j], args.sketch) >= args.min_jaccard:
                union(i, j)
        if (i + 1) % 250 == 0:
            print(f"[cluster] {i + 1}/{n}", flush=True)

    df["clade"] = [find(i) for i in range(n)]
    df["clade"] = df.groupby("clade").ngroup()
    df.to_csv(out / "elements_clades.csv", index=False)

    sizes = df.clade.value_counts()
    # A genome is usable for a clade-held-out split only if all of its elements
    # fall on the same side of the split; genomes spanning both sides have to be
    # dropped. This is the number that decides whether 2a is affordable.
    per_genome_clades = df.groupby("genome_id").clade.nunique()
    summary = {
        "n_elements": int(len(df)),
        "n_genomes_with_elements": int(df.genome_id.nunique()),
        "n_clades": int(df.clade.nunique()),
        "largest_clade": int(sizes.iloc[0]),
        "singleton_clades": int((sizes == 1).sum()),
        "median_clade_size": float(sizes.median()),
        "genomes_with_multiple_clades": int((per_genome_clades > 1).sum()),
        "max_clades_in_one_genome": int(per_genome_clades.max()),
        "params": {"k": args.k, "sketch": args.sketch, "min_jaccard": args.min_jaccard},
        "next_step": (
            "Inspect n_clades and genomes_with_multiple_clades before submitting "
            "any retraining. A clade-held-out split strands every genome that "
            "hosts elements on both sides; if that is most of the corpus, the "
            "design needs rethinking (e.g. hold out clades that co-occur least) "
            "rather than 3 GPU-days."
        ),
    }
    (out / "clade_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
