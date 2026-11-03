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
cd /work/$USER/mge_explorer
IN="${1:-raw/bacterial_ftp_paths.txt}"
mkdir -p raw/genomes logs
grep -E '\.fna\.gz$' "$IN" > raw/urls.txt
cat raw/urls.txt | parallel --bar -j16 'wget -qc --retry-connrefused --waitretry=5 --tries=5 --timeout=60 -nd -P raw/genomes {} || echo "FAIL {}"'
find raw/genomes -name '*.gz' -print0 | xargs -0 -r gunzip -f
find raw/genomes -name '*.fna' | sort -u > genome_list.txt
rm -rf genome_chunks && mkdir -p genome_chunks
split -d -l 1000 genome_list.txt genome_chunks/GL_
echo "genomes:" $(wc -l < genome_list.txt)
echo "chunks :" $(ls genome_chunks/GL_* | wc -l)
