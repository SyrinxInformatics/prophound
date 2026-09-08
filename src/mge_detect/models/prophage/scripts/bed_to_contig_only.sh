ALL_BED=windows/bed/all_win.uni.bed
FA=raw/all_genomes.fa

awk -F'\t' 'BEGIN{OFS="\t"}{split($1,a,"|"); print a[2],$2,$3,$4}' "$ALL_BED" > windows/bed/all_win.contig.bed

cut -f1 windows/bed/all_win.contig.bed | sort -u > tmp/bed.contigs
grep '^>' "$FA" | cut -d' ' -f1 | sed 's/^>//' | sort -u > tmp/fa.contigs
comm -23 tmp/bed.contigs tmp/fa.contigs | wc -l    # expect 0

samtools faidx "$FA"
cut -f1,2 "$FA".fai > labels/all.contig.sizes
bedtools sort -g labels/all.contig.sizes -i windows/bed/all_win.contig.bed > windows/bed/all_win.contig.sorted.bed

bedtools getfasta -fi "$FA" -bed windows/bed/all_win.contig.sorted.bed -name -fo windows/fasta/all_512k.fa

grep -c '^>' windows/fasta/all_512k.fa   # expect 7940
