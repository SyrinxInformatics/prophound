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

import random
from pathlib import Path

base = Path.cwd()
contig_bed = base / "windows/bed/all_win.contig.sorted.bed"
orig_bed = base / "windows/bed/all_win.uni.bed"
labels_out = base / "windows/labels/all_512k.labels.tsv"
splits_dir = base / "splits"
tmp_dir = base / "tmp"

splits_dir.mkdir(parents=True, exist_ok=True)
labels_out.parent.mkdir(parents=True, exist_ok=True)
tmp_dir.mkdir(parents=True, exist_ok=True)

with open(contig_bed) as fin, open(labels_out, "w") as fout:
    for line in fin:
        if not line.strip():
            continue
        f = line.rstrip().split("\t")
        name = f[3]
        y = 1 if name.startswith("pos|") else 0
        fout.write(f"{name}\t{y}\n")

contigs_in_contig_bed = []
with open(contig_bed) as fin:
    for line in fin:
        if not line.strip():
            continue
        contigs_in_contig_bed.append(line.split("\t", 1)[0])
contigs_in_contig_bed = set(contigs_in_contig_bed)

contig_to_genome = {}
with open(orig_bed) as fin:
    for line in fin:
        if not line.strip():
            continue
        chrom = line.split("\t", 1)[0]
        parts = chrom.split("|", 1)
        if len(parts) != 2:
            continue
        genome, contig = parts
        if contig in contigs_in_contig_bed and contig not in contig_to_genome:
            contig_to_genome[contig] = genome

genomes = sorted(
    set(contig_to_genome[c] for c in contigs_in_contig_bed if c in contig_to_genome)
)
rng = random.Random(42)
rng.shuffle(genomes)
n = len(genomes)
val = n // 10
train = n - 2 * val
train_gen = set(genomes[:train])
val_gen = set(genomes[train : train + val])
test_gen = set(genomes[train + val :])

(base / "splits/train.genomes").write_text("\n".join(sorted(train_gen)) + "\n")
(base / "splits/val.genomes").write_text("\n".join(sorted(val_gen)) + "\n")
(base / "splits/test.genomes").write_text("\n".join(sorted(test_gen)) + "\n")


def write_split_bed(gen_set, out_path):
    with open(contig_bed) as fin, open(out_path, "w") as fout:
        for line in fin:
            if not line.strip():
                continue
            contig = line.split("\t", 1)[0]
            g = contig_to_genome.get(contig)
            if g in gen_set:
                fout.write(line)


write_split_bed(train_gen, base / "splits/train.contig.bed")
write_split_bed(val_gen, base / "splits/val.contig.bed")
write_split_bed(test_gen, base / "splits/test.contig.bed")
print("done")

