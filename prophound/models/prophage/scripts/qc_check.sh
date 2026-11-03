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
if command -v mawk >/dev/null; then AWK=mawk; else AWK=awk; fi
export LC_ALL=C

ALLC="windows/bed/all_win.contig.sorted.bed"
FA="windows/fasta/all_512k.fa"
LAB="windows/labels/all_512k.labels.tsv"
TR="splits/train.contig.bed"
VA="splits/val.contig.bed"
TE="splits/test.contig.bed"
SIZES="labels/all.contig.sizes"
POS_GENOMAD="positives/geNomad.norm.bed"
POS_CONTIG="positives/geNomad.norm.contig.sorted.bed"
TMP="tmp/qc"; mkdir -p "$TMP"

for f in "$ALLC" "$FA" "$LAB" "$TR" "$VA" "$TE"; do [ -s "$f" ] || { echo "MISSING $f"; exit 1; }; done
command -v bedtools >/dev/null || { echo "MISSING bedtools"; exit 1; }

if [ ! -s "$SIZES" ]; then
  if [ -s "${FA}.fai" ]; then cut -f1,2 "${FA}.fai" > "$SIZES"
  elif command -v samtools >/dev/null; then samtools faidx "$FA"; cut -f1,2 "${FA}.fai" > "$SIZES"
  else echo "MISSING $SIZES and samtools"; exit 1; fi
fi

wc -l "$ALLC" | awk '{print "bed_windows="$1}'
grep -c '^>' "$FA" | awk '{print "fasta_windows="$1}'
wc -l "$LAB" | awk '{print "label_rows="$1}'

cut -f4 "$ALLC" | sort > "$TMP/all.ids"
grep -E '^>' "$FA" | sed 's/^>//' | cut -d' ' -f1 | sed 's/:.*$//' | sort > "$TMP/fa.ids"
cut -f1 "$LAB" | sort > "$TMP/lab.ids"
D1=$(comm -3 "$TMP/all.ids" "$TMP/fa.ids" | wc -l); D2=$(comm -3 "$TMP/all.ids" "$TMP/lab.ids" | wc -l)
echo "idset_diff_bed_vs_fa=$D1"
echo "idset_diff_bed_vs_lab=$D2"
[ "$D1" -eq 0 ] && [ "$D2" -eq 0 ] || { echo "ID MISMATCH"; exit 1; }

awk -F'\t' '($2==0||$2==1){ok++} END{print "label_valid_rows="ok}' "$LAB"
POS=$(awk -F'\t' '$2==1{c++} END{print c+0}' "$LAB"); NEG=$(awk -F'\t' '$2==0{c++} END{print c+0}' "$LAB")
echo "labels_pos=$POS labels_neg=$NEG"

BAD=$($AWK 'BEGIN{bad=0;n=0} /^>/{if(n>0&&len!=512001)bad++; len=0; n++; next} {gsub(/\r/,""); if(!/^>/)len+=length($0)} END{if(n>0&&len!=512001)bad++; print bad}' "$FA")
echo "fasta_len_ne_512001=$BAD"

$AWK 'BEGIN{A=C=G=T=N=0} /^>/{next} {A+=gsub(/[Aa]/,""); C+=gsub(/[Cc]/,""); G+=gsub(/[Gg]/,""); T+=gsub(/[Tt]/,""); N+=gsub(/[^ACGTacgt]/,"")} END{tot=A+C+G+T+N; printf("gc_frac=%.4f n_frac=%.4f\n", (tot?(G+C)/tot:0), (tot?N/tot:0))}' "$FA"

cut -f4 "$TR" | sort > "$TMP/tr.ids"
cut -f4 "$VA" | sort > "$TMP/va.ids"
cut -f4 "$TE" | sort > "$TMP/te.ids"
TRN=$(wc -l < "$TMP/tr.ids"); VAN=$(wc -l < "$TMP/va.ids"); TEN=$(wc -l < "$TMP/te.ids")
echo "split_sizes train=$TRN val=$VAN test=$TEN"
comm -12 "$TMP/tr.ids" "$TMP/va.ids" | wc -l | awk '{print "split_overlap_train_val="$1}'
comm -12 "$TMP/tr.ids" "$TMP/te.ids" | wc -l | awk '{print "split_overlap_train_test="$1}'
comm -12 "$TMP/va.ids" "$TMP/te.ids" | wc -l | awk '{print "split_overlap_val_test="$1}'
cat "$TMP/tr.ids" "$TMP/va.ids" "$TMP/te.ids" | sort > "$TMP/all_split.ids"
comm -3 "$TMP/all.ids" "$TMP/all_split.ids" | wc -l | awk '{print "split_union_diff_all="$1}'

if [ ! -s "$POS_CONTIG" ]; then
  awk 'BEGIN{OFS="\t"}{split($1,a,"|"); print a[2],$2,$3,$4}' "$POS_GENOMAD" \
  | bedtools sort -g "$SIZES" -i - > "$POS_CONTIG"
fi
awk -F'\t' '$4 !~ /^pos\|/ {print}' "$ALLC" > "$TMP/neg_windows.contig.bed"
bedtools intersect -u -a "$TMP/neg_windows.contig.bed" -b "$POS_CONTIG" | wc -l | awk '{print "neg_overlap_with_positive="$1}'

FAIL=0
[ "$(awk '{print $1}' <<<"$(wc -l < "$ALLC")")" -eq "$(grep -c '^>' "$FA")" ] || FAIL=1
[ "$(wc -l < "$LAB")" -eq "$(wc -l < "$ALLC")" ] || FAIL=1
[ "$D1" -eq 0 ] && [ "$D2" -eq 0 ] || FAIL=1
[ "$(awk '{print $1}' <<<"$(awk -F'\t' '($2==0||$2==1){ok++} END{print ok+0}' "$LAB")")" -eq "$(wc -l < "$LAB")" ] || FAIL=1
OV1=$(comm -12 "$TMP/tr.ids" "$TMP/va.ids" | wc -l); OV2=$(comm -12 "$TMP/tr.ids" "$TMP/te.ids" | wc -l); OV3=$(comm -12 "$TMP/va.ids" "$TMP/te.ids" | wc -l)
[ "$OV1" -eq 0 ] && [ "$OV2" -eq 0 ] && [ "$OV3" -eq 0 ] || FAIL=1
[ "$(comm -3 "$TMP/all.ids" "$TMP/all_split.ids" | wc -l)" -eq 0 ] || FAIL=1
[ "$(bedtools intersect -u -a "$TMP/neg_windows.contig.bed" -b "$POS_CONTIG" | wc -l)" -eq 0 ] || FAIL=1

if [ "$FAIL" -eq 0 ]; then echo "QC_PASS"; else echo "QC_FAIL"; exit 1; fi

