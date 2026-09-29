WIN=256000

mkdir -p windows/bed windows/meta

awk -F'\t' -v OFS='\t' '{mid=int(($2+$3)/2); print $1, mid, mid+1, "pos|"NR}' positives/geNomad.norm.bed \
| bedtools slop -g labels/all.sizes -b $WIN \
| sort -k1,1 -k2,2n > windows/bed/pos_win.bed
NPOS=$(wc -l < windows/bed/pos_win.bed)
echo "positive windows: $NPOS"

cat labels/negatives/*.neg.bed > labels/negatives/all.neg.bed
bedtools shuffle -i windows/bed/pos_win.bed -g labels/all.sizes \
                 -incl labels/negatives/all.neg.bed -noOverlapping -chrom \
| awk -F'\t' -v OFS='\t' '{split($4,a,"|"); print $1,$2,$3,"neg|"a[2]}' \
> windows/bed/neg_win.bed
NNEG=$(wc -l < windows/bed/neg_win.bed)
echo "negative windows: $NNEG"

cat windows/bed/pos_win.bed windows/bed/neg_win.bed \
| sort -k1,1 -k2,2n > windows/bed/all_win.bed

