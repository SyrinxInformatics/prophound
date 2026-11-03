#!/usr/bin/env bash
set -euo pipefail

ALL="windows/bed/all_win.contig.sorted.bed"
SIZES="labels/all.contig.sizes"
OUTDIR="windows/bed"
mkdir -p "$OUTDIR" tmp/clean

awk '$4 ~ /^pos\|/' "$ALL" > tmp/clean/pos.bed
awk '$4 !~ /^pos\|/' "$ALL" > tmp/clean/neg.bed

bedtools sort -g "$SIZES" -i tmp/clean/pos.bed > tmp/clean/pos.sorted.bed
bedtools sort -g "$SIZES" -i tmp/clean/neg.bed > tmp/clean/neg.sorted.bed

bedtools intersect -v -a tmp/clean/neg.sorted.bed -b tmp/clean/pos.sorted.bed > tmp/clean/neg.noPos.bed

bedtools cluster -i tmp/clean/neg.noPos.bed \
| awk 'BEGIN{OFS="\t"}{if(!seen[$NF]++){print $1,$2,$3,$4}}' \
> tmp/clean/neg.cleaned.bed

cat tmp/clean/pos.sorted.bed tmp/clean/neg.cleaned.bed \
| bedtools sort -g "$SIZES" -i - \
> "$OUTDIR/all_win.contig.cleaned.bed"

mv tmp/clean/pos.sorted.bed  "$OUTDIR/pos_win.cleaned.bed"
mv tmp/clean/neg.cleaned.bed "$OUTDIR/neg_win.cleaned.bed"

POS_N=$(wc -l < "$OUTDIR/pos_win.cleaned.bed")
NEG_N=$(wc -l < "$OUTDIR/neg_win.cleaned.bed")
ALL_N=$(wc -l < "$OUTDIR/all_win.contig.cleaned.bed")
echo "kept_positives=$POS_N kept_negatives=$NEG_N total_cleaned=$ALL_N"

bedtools intersect -u -a "$OUTDIR/neg_win.cleaned.bed" -b "$OUTDIR/pos_win.cleaned.bed" | wc -l | awk '{print "neg_vs_pos_overlaps="$1}'

bedtools intersect -u -a "$OUTDIR/neg_win.cleaned.bed" -b "$OUTDIR/neg_win.cleaned.bed" -sorted \
| awk 'NR>1{c=1} END{print "neg_vs_neg_overlaps=" (c?c-1:0)}'

