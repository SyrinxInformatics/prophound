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

ALL_BED=windows/bed/all_win.uni.bed
FA=raw/all_genomes.fa

awk -F'\t' 'BEGIN{OFS="\t"}{split($1,a,"|"); print a[2],$2,$3,$4}' "$ALL_BED" > windows/bed/all_win.contig.bed

cut -f1 windows/bed/all_win.contig.bed | sort -u > tmp/bed.contigs
grep '^>' "$FA" | cut -d' ' -f1 | sed 's/^>//' | sort -u > tmp/fa.contigs
comm -23 tmp/bed.contigs tmp/fa.contigs | wc -l    # expect 0

samtools faidx "$FA"
cut -f1,2 "$FA".fai > labels/all.contig.sizes
bedtools sort -g labels/all.contig.sizes -i windows/bed/all_win.contig.bed > windows/bed/all_win.contig.sorted.bed

bedtools getfasta -fi "$FA" -bed windows/bed/all_win.contig.sorted.bed -name -fo windows/fasta/all_512k.fa

grep -c '^>' windows/fasta/all_512k.fa   # expect 7940
