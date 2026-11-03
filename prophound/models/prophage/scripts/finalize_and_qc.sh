#!/usr/bin/env bash
set -euo pipefail

ALLC="windows/bed/all_win.contig.cleaned.bed"
FA="windows/fasta/all_512k.cleaned.fa"
LAB="windows/labels/all_512k.cleaned.labels.tsv"
TR="splits/train.contig.cleaned.bed"
VA="splits/val.contig.cleaned.bed"
TE="splits/test.contig.cleaned.bed"
POS_BED="windows/bed/pos_win.cleaned.bed"
SOFT_NEG_BED="windows/bed/neg_win.cleaned.bed"
POS_GENOMAD="positives/geNomad.norm.bed"
SIZES="labels/all.contig.sizes"
HN_SCRIPT="01_scripts/make_hard_negatives.py"

REQ=( "$ALLC" "$FA" "$LAB" "$TR" "$VA" "$TE" "$POS_BED" "$SOFT_NEG_BED" "$POS_GENOMAD" "$HN_SCRIPT" )
for f in "${REQ[@]}"; do [ -s "$f" ] || { echo "MISSING $f"; exit 1; }; done
command -v bedtools >/dev/null || { echo "MISSING bedtools"; exit 1; }

mkdir -p tmp/qc windows/bed windows/fasta windows/labels windows/manifest release
TMP="tmp/qc"
DATE_TAG=$(date +"%Y%m%d_%H%M%S")
OUTDIR="release/session_${DATE_TAG}"
mkdir -p "$OUTDIR"

if [ ! -s "$SIZES" ]; then
  if [ ! -s "${FA}.fai" ]; then
    command -v samtools >/dev/null || { echo "MISSING samtools for faidx"; exit 1; }
    samtools faidx "$FA"
  fi
  cut -f1,2 "${FA}.fai" > "$SIZES"
fi

python3 "$HN_SCRIPT"
bedtools sort -g "$SIZES" -i windows/bed/neg_win.hard.strict.bed > "$TMP/neg_hard.sorted.bed"
mv "$TMP/neg_hard.sorted.bed" windows/bed/neg_win.hard.strict.bed

bedtools intersect -u -a windows/bed/neg_win.hard.strict.bed -b "$POS_BED" | wc -l | awk '{print "hard_vs_pos_overlaps="$1}' > "$TMP/hard_ov.txt"
bedtools intersect -a windows/bed/neg_win.hard.strict.bed -b windows/bed/neg_win.hard.strict.bed -wo \
| awk '($1!=$5)||($2!=$6)||($3!=$7)||($4!=$8)' | wc -l | awk '{print "hard_vs_hard_overlaps="$1}' >> "$TMP/hard_ov.txt"

cat "$SOFT_NEG_BED" windows/bed/neg_win.hard.strict.bed \
| bedtools sort -g "$SIZES" -i - > windows/bed/neg_win.with_hard.bed

bedtools intersect -u -a windows/bed/neg_win.with_hard.bed -b "$POS_BED" | wc -l | awk '{print "allneg_vs_pos_overlaps="$1}' > "$TMP/allneg_ov.txt"
bedtools intersect -a windows/bed/neg_win.with_hard.bed -b windows/bed/neg_win.with_hard.bed -wo \
| awk '($1!=$5)||($2!=$6)||($3!=$7)||($4!=$8)' | wc -l | awk '{print "allneg_vs_allneg_overlaps="$1}' >> "$TMP/allneg_ov.txt"

awk -F'\t' 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}' "$POS_BED" > "$TMP/pos.clean.bed"
awk -F'\t' 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}' windows/bed/neg_win.with_hard.bed > "$TMP/neg.clean.bed"
cat "$TMP/pos.clean.bed" "$TMP/neg.clean.bed" \
| bedtools sort -g "$SIZES" -i - > windows/bed/all_win.contig.final.bed

bedtools intersect -a windows/bed/all_win.contig.final.bed -b windows/bed/all_win.contig.final.bed -wo \
| awk '($1!=$5)||($2!=$6)||($3!=$7)||($4!=$8)' | wc -l > "$TMP/final_self_ov.n"

bedtools getfasta -fi raw/all_genomes.fa -bed windows/bed/all_win.contig.final.bed -name -fo windows/fasta/all_512k.final.fa
awk -F'\t' 'BEGIN{OFS="\t"}{print $4, ($4 ~ /^pos\|/)?1:0}' windows/bed/all_win.contig.final.bed > windows/labels/all_512k.final.labels.tsv

awk '{print $4,"\ttrain"}' "$TR" > "$TMP/split.map"
awk '{print $4,"\tval"}'   "$VA" >> "$TMP/split.map"
awk '{print $4,"\ttest"}'  "$TE" >> "$TMP/split.map"
awk 'NR==FNR{m[$1]=$2;next}{sp=m[$4]; if(sp=="train")print > "splits/train.contig.final.bed"; else if(sp=="val") print > "splits/val.contig.final.bed"; else if(sp=="test")print > "splits/test.contig.final.bed"}' \
  "$TMP/split.map" windows/bed/all_win.contig.final.bed

echo "bed_windows=$(wc -l < "$ALLC")"
echo "fasta_windows=$(grep -c '^>' "$FA")"
echo "label_rows=$(wc -l < "$LAB")"

cut -f4 "$ALLC" | sort > "$TMP/all.ids"
grep -E '^>' "$FA" | sed 's/^>//; s/:.*$//' | sort > "$TMP/fa.ids"
cut -f1 "$LAB" | sort > "$TMP/lab.ids"
D1=$(comm -3 "$TMP/all.ids" "$TMP/fa.ids" | wc -l)
D2=$(comm -3 "$TMP/all.ids" "$TMP/lab.ids" | wc -l)
echo "idset_diff_bed_vs_fa=$D1"
echo "idset_diff_bed_vs_lab=$D2"
[ "$D1" -eq 0 ] && [ "$D2" -eq 0 ] || { echo "ID MISMATCH"; exit 1; }

awk -F'\t' '($2==0||$2==1){ok++} END{print "label_valid_rows="ok}' "$LAB"
POS=$(awk -F'\t' '$2==1{c++} END{print c+0}' "$LAB")
NEG=$(awk -F'\t' '$2==0{c++} END{print c+0}' "$LAB")
echo "labels_pos=$POS labels_neg=$NEG"

awk 'BEGIN{bad=0;n=0} /^>/{if(n>0&&len!=512001)bad++; len=0; n++; next} {gsub(/\r/,""); if(!/^>/)len+=length($0)} END{if(n>0&&len!=512001)bad++; print "fasta_len_ne_512001="bad}' "$FA"
awk 'BEGIN{A=C=G=T=N=0} /^>/{next} {A+=gsub(/[Aa]/,""); C+=gsub(/[Cc]/,""); G+=gsub(/[Gg]/,""); T+=gsub(/[Tt]/,""); N+=gsub(/[^ACGTacgt]/,"")} END{tot=A+C+G+T+N; printf("gc_frac=%.4f n_frac=%.4f\n", (tot?(G+C)/tot:0), (tot?N/tot:0))}' "$FA"

bedtools sort -g "$SIZES" -i "$ALLC" > "$TMP/all.sorted.bed"
bedtools intersect -a "$TMP/all.sorted.bed" -b "$TMP/all.sorted.bed" -wo \
| awk '($1!=$5)||($2!=$6)||($3!=$7)||($4!=$8)' > "$TMP/self_overlaps.tsv"
SO=$(wc -l < "$TMP/self_overlaps.tsv")
echo "self_overlaps_all=$SO"

awk '$4 ~ /^pos\|/' "$TMP/all.sorted.bed" > "$TMP/pos.bed"
awk '$4 !~ /^pos\|/' "$TMP/all.sorted.bed" > "$TMP/neg.bed"
bedtools intersect -u -a "$TMP/neg.bed" -b "$TMP/pos.bed" | wc -l | awk '{print "neg_vs_pos_overlaps="$1}'
bedtools intersect -a "$TMP/neg.bed" -b "$TMP/neg.bed" -wo \
| awk '($1!=$5)||($2!=$6)||($3!=$7)||($4!=$8)' | wc -l | awk '{print "neg_vs_neg_overlaps="$1}'
bedtools intersect -a "$TMP/pos.bed" -b "$TMP/pos.bed" -wo \
| awk '($1!=$5)||($2!=$6)||($3!=$7)||($4!=$8)' | wc -l | awk '{print "pos_vs_pos_overlaps="$1}'

for SPL in train val test; do
  BED="splits/${SPL}.contig.cleaned.bed"
  [ -s "$BED" ] || { echo "MISSING $BED"; exit 1; }
  cut -f4 "$BED" > "$TMP/${SPL}.ids"
  N=$(wc -l < "$TMP/${SPL}.ids")
  P=$(awk 'NR==FNR{lab[$1]=$2;next} {if(lab[$1]==1) c++} END{print c+0}' "$LAB" "$TMP/${SPL}.ids")
  echo "${SPL}_n=$N ${SPL}_pos=$P ${SPL}_neg=$(($N-$P))"
done

comm -12 "$TMP/train.ids" "$TMP/val.ids"  | wc -l | awk '{print "split_overlap_train_val="$1}'
comm -12 "$TMP/train.ids" "$TMP/test.ids" | wc -l | awk '{print "split_overlap_train_test="$1}'
comm -12 "$TMP/val.ids"   "$TMP/test.ids" | wc -l | awk '{print "split_overlap_val_test="$1}'

HN=$(wc -l < windows/bed/neg_win.hard.strict.bed)
SN=$(awk 'BEGIN{n=0} {if($4 !~ /^negH\|/) n++} END{print n}' windows/bed/neg_win.with_hard.bed)
TP=$(wc -l < "$POS_BED")
AN=$(wc -l < windows/bed/neg_win.with_hard.bed)

{
  echo "# Final QC summary ${DATE_TAG}"
  echo "bed_windows=$(wc -l < "$ALLC")"
  echo "fasta_windows=$(grep -c '^>' "$FA")"
  echo "label_rows=$(wc -l < "$LAB")"
  echo "idset_diff_bed_vs_fa=$D1"
  echo "idset_diff_bed_vs_lab=$D2"
  echo "labels_pos=$POS labels_neg=$NEG"
  awk 'BEGIN{bad=0;n=0} /^>/{if(n>0&&len!=512001)bad++; len=0; n++; next} {gsub(/\r/,""); if(!/^>/)len+=length($0)} END{if(n>0&&len!=512001)bad++; print "fasta_len_ne_512001="bad}' "$FA"
  awk 'BEGIN{A=C=G=T=N=0} /^>/{next} {A+=gsub(/[Aa]/,""); C+=gsub(/[Cc]/,""); G+=gsub(/[Gg]/,""); T+=gsub(/[Tt]/,""); N+=gsub(/[^ACGTacgt]/,"")} END{tot=A+C+G+T+N; printf("gc_frac=%.4f n_frac=%.4f\n", (tot?(G+C)/tot:0), (tot?N/tot:0))}' "$FA"
  echo "self_overlaps_all=$SO"
  echo "hard_negatives=$HN soft_negatives=$SN positives=$TP all_negatives=$AN"
  cat "$TMP/hard_ov.txt"
  cat "$TMP/allneg_ov.txt"
  echo "final_self_overlaps=$(cat "$TMP/final_self_ov.n")"
  for SPL in train val test; do
    BED="splits/${SPL}.contig.cleaned.bed"
    N=$(wc -l < "$TMP/${SPL}.ids"); P=$(awk 'NR==FNR{lab[$1]=$2;next} {if(lab[$1]==1) c++} END{print c+0}' "$LAB" "$TMP/${SPL}.ids")
    echo "${SPL}_n=$N ${SPL}_pos=$P ${SPL}_neg=$(($N-$P))"
  done
  echo "split_overlap_train_val=$(comm -12 "$TMP/train.ids" "$TMP/val.ids" | wc -l)"
  echo "split_overlap_train_test=$(comm -12 "$TMP/train.ids" "$TMP/test.ids" | wc -l)"
  echo "split_overlap_val_test=$(comm -12 "$TMP/val.ids" "$TMP/test.ids" | wc -l)"
} | tee "$OUTDIR/QC_SUMMARY.txt" >/dev/null

cp -a \
  windows/bed/pos_win.cleaned.bed \
  windows/bed/neg_win.cleaned.bed \
  windows/bed/neg_win.hard.strict.bed \
  windows/bed/neg_win.with_hard.bed \
  windows/bed/all_win.contig.final.bed \
  windows/fasta/all_512k.final.fa \
  windows/labels/all_512k.final.labels.tsv \
  splits/train.contig.final.bed \
  splits/val.contig.final.bed \
  splits/test.contig.final.bed \
  "$OUTDIR"/

[ -s windows/manifest/all_512k_manifest.csv ] && cp windows/manifest/all_512k_manifest.csv "$OUTDIR"/ || true
cp -a positives/geNomad.norm.bed "$OUTDIR"/
[ -s positives/geNomad.norm.contig.sorted.bed ] && cp positives/geNomad.norm.contig.sorted.bed "$OUTDIR"/ || true

echo "Artifacts packaged in: $OUTDIR"

