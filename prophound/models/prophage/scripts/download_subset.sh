#!/usr/bin/env bash
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

set -euo pipefail

echo "[A] Fetching RefSeq assembly summary…"
mkdir -p raw
if [ ! -f raw/assembly_summary_refseq.txt ]; then
    wget -q -O raw/assembly_summary_refseq.txt \
        https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt
fi

echo "[B] Selecting one representative complete bacterial genome per genus…"
awk -F'\t' '
BEGIN {OFS="\t"}
$1 !~ /^#/ && $12=="Complete Genome" && $25=="bacteria" {
    if (!seen[$8]++) {
        gsub(/ftp:/,"https:",$20); gsub(/\/$/,"",$20);
        print $20 "/" $1 "_" $16 "_genomic.fna.gz"
    }
}' raw/assembly_summary_refseq.txt > subset_urls.txt

echo "Will download $(wc -l < subset_urls.txt) genomes."

echo "[B2] Downloading genomes…"
mkdir -p raw/genomes
cat subset_urls.txt | parallel --bar -j8 \
   wget -qc --retry-connrefused --waitretry=5 --tries=5 --timeout=60 \
        -nd -P raw/genomes {}

echo "[B3] Decompressing genomes…"
find raw/genomes -name '*.gz' -print0 | xargs -0 gunzip

echo "[C] Writing genome_list.txt…"
find raw/genomes -name '*.fna' | sort -u > genome_list.txt
echo "Genome list contains $(wc -l < genome_list.txt) entries"

echo "Representative genomes downloaded and ready."
