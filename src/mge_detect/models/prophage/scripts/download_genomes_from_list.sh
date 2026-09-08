#!/usr/bin/env bash
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
