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
"""Shred genomes into short contigs to simulate a fragmented assembly.

docs/adaptive-multiscale-detection.md Part 4: `--short-contigs scan` scores
contigs under 10 kb at 5 kb / 3 kb windows, whose thresholds cannot be borrowed
from the 10 kb operating point (negatives drift up and per-window sensitivity
falls as windows shrink). The eval corpus has almost no short contigs, so this
makes them: every contig is cut into consecutive fragments whose lengths are
drawn uniformly from [min_len, max_len), seeded, covering the whole sequence.
The final remainder of each contig is kept whatever its length, so the
"too short to scan" path is exercised too.

Fragment ids are `<contig>__<offset>`; short_contig_sweep.py maps calls back
to genome coordinates with that offset.

Example:
  python scripts/shred_genomes.py --genome_list scripts/eval_genome_list_sizestrat.txt \\
      --genome_dir ../raw/genomes --out_dir analysis/short_contig_calib/fasta
"""

from __future__ import annotations

import argparse
import random
from pathlib import Path

from pyfaidx import Fasta

SEP = "__"


def shred(length: int, rng: random.Random, min_len: int, max_len: int):
    pos = 0
    while pos < length:
        n = rng.randrange(min_len, max_len)
        yield pos, min(length, pos + n)
        pos += n


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genome_list", required=True)
    ap.add_argument("--genome_dir", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--min_len", type=int, default=3000)
    ap.add_argument("--max_len", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    genomes = [l.strip() for l in open(args.genome_list) if l.strip()]
    n_frag = 0
    for g in genomes:
        rng = random.Random(f"{args.seed}:{g}")   # per-genome seed: stable if the list changes
        fa = Fasta(str(Path(args.genome_dir) / f"{g}.fna"), as_raw=True, sequence_always_upper=True)
        with open(out / f"{g}.fna", "w") as fh:
            for contig in fa.keys():
                seq = fa[contig][:]
                for s, e in shred(len(seq), rng, args.min_len, args.max_len):
                    fh.write(f">{contig}{SEP}{s}\n{seq[s:e]}\n")
                    n_frag += 1
        fa.close()
    print(f"[done] {len(genomes)} genomes -> {n_frag} fragments in {out}", flush=True)


if __name__ == "__main__":
    main()
