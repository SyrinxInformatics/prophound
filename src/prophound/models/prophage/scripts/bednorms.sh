awk -F'\t' -v OFS='\t' '
{
  # $1 = genome|seq_name     where seq_name may be "contig" or "contig|provirus_..."
  split($1,a,"|");
  genome=a[1];
  contig=a[2];
  prov=(a[3] ~ /^provirus_/) ? a[3] : "NA";
  name = (prov=="NA" ? "geNomad" : "geNomad|" prov);
  print genome "|" contig, $2, $3, name
}' positives/geNomad.bed \
| sort -k1,1 -k2,2n > positives/geNomad.norm.bed

echo "Normalized: $(wc -l < positives/geNomad.norm.bed)"
head positives/geNomad.norm.bed

