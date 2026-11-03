#!/usr/bin/env bash
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
