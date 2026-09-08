cd /work/pke835/mge_explorer

find genomad_out -type f -name '*_provirus.tsv' -print0 \
| while IFS= read -r -d '' t; do
    base=$(basename "$(dirname "$(dirname "$t")")")
    awk -F'\t' -v OFS='\t' -v G="$base" '
      FNR==1{
        seq=1; b=2; e=3;
        for(i=1;i<=NF;i++){
          if($i ~ /^(seq_name|sequence|contig)$/) seq=i
          if($i ~ /^(start|begin)$/)               b=i
          if($i ~ /^(end|stop)$/)                  e=i
        }
        next
      }
      ($b~/^[0-9]+$/ && $e~/^[0-9]+$/){
        print G "|" $seq, $b-1, $e, "geNomad"
      }' "$t"
  done \
| sort -k1,1 -k2,2n > positives/geNomad.bed

echo "Lines in positives/geNomad.bed: $(wc -l < positives/geNomad.bed)"
head positives/geNomad.bed

