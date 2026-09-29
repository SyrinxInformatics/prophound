cd /work/$USER/mge_explorer/plasmids/chunks
mkdir -p reruns
for CH in files_chunk_0001.list files_chunk_0002.list files_chunk_0003.list files_chunk_0004.list; do
  CHFILE="$CH"
  DONE="reruns/${CH}.done.txt"
  REM="reruns/${CH}.remain.txt"
  : > "$DONE"; : > "$REM"
  while read -r g; do
    b=$(basename "$g" .fna)
    t="/work/$USER/mge_explorer/plasmids/genome_chunks/genomad_out/$b/${b}_find_proviruses/${b}_provirus.tsv"
    if [ -s "$t" ]; then echo "$g" >> "$DONE"; else echo "$g" >> "$REM"; fi
  done < "$CHFILE"
  echo "$CH  done: $(wc -l < "$DONE")   remain: $(wc -l < "$REM")"
done
