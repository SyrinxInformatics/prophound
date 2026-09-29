#!/bin/bash
# Inventory the HPCC for the scientific-validity program (docs/scientific-validity.md).
# Read-only. Run on the login node, no SLURM needed:
#   ssh pke835@arc.utsa.edu 'bash -s' < validity_preflight.sh   # (or push it and run it)
# Purpose: find out what exists and what is installed BEFORE committing CPU/GPU
# hours, so the analysis jobs don't fail deep.
set -uo pipefail
B=/work/$USER/mge_explorer
D=$B/data_refactored

echo "=== roots ==="
for p in "$B" "$D"; do [ -d "$p" ] && echo "OK   $p" || echo "MISS $p"; done

echo; echo "=== inputs the validity jobs need ==="
for p in "$D/data/manifest_windows.filtered.csv" "$D/splits/train.txt" "$D/splits/val.txt" \
         "$D/splits/test.txt" "$D/sequences" "$D/beds_pos" "$D/beds_all" "$B/raw/genomes" "$B/genomad_out"; do
  # NOTE: no `du -sh` here. It walks the whole tree, and on genomad_out (one
  # subdirectory per genome) that takes long enough to look like a hang.
  if [ -e "$p" ]; then echo "OK   $p"; else echo "MISS $p"; fi
done

echo; echo "=== manifest schema + size ==="
# `wc -l` on a 1.3G csv is slow; the header is what matters, and rows are
# estimated from file size / mean line length of the first 1000 lines.
[ -f "$D/data/manifest_windows.filtered.csv" ] && {
  head -1 "$D/data/manifest_windows.filtered.csv"
  MF="$D/data/manifest_windows.filtered.csv"
  BYTES=$(stat -c %s "$MF" 2>/dev/null || stat -f %z "$MF")
  AVG=$(head -1000 "$MF" | wc -c)
  echo "approx rows: $(( BYTES / (AVG / 1000) ))  (estimated from size, not counted)"
}

echo; echo "=== positives: where element coordinates live ==="
echo "beds_pos files: $(ls "$D/beds_pos" 2>/dev/null | wc -l)"
echo "sequences/ sample (confirms --fasta_suffix):"; ls "$D/sequences" 2>/dev/null | head -3
ls "$D/beds_pos" 2>/dev/null | head -2
[ -d "$D/beds_pos" ] && head -3 "$D/beds_pos/$(ls "$D/beds_pos" | head -1)" 2>/dev/null

echo; echo "=== geNomad outputs (for 1b annotation; gene tables carry marker calls) ==="
# -print -quit stops at the first hit instead of walking the whole tree.
G1=$(find "$B/genomad_out" -name "*_genes.tsv" -print -quit 2>/dev/null)
if [ -n "$G1" ]; then echo "example: $G1"; echo "columns:"; head -1 "$G1"; else echo "no *_genes.tsv found (1b needs another annotation source)"; fi

echo; echo "=== conda env + python libs ==="
module load miniconda 2>/dev/null
conda activate /home/$USER/.conda/envs/mge_explorer_gpu 2>/dev/null && echo "env activated: $(which python)"
python - <<'PY' 2>&1
for m in ("numpy", "pandas", "pyfaidx", "sklearn", "scipy", "yaml", "torch"):
    try:
        mod = __import__(m)
        print(f"  OK   {m} {getattr(mod, '__version__', '')}")
    except Exception as e:
        print(f"  MISS {m}  ({type(e).__name__})")
PY

echo; echo "=== optional external tools (we avoid depending on these, but note them) ==="
for t in mmseqs mash skani fastANI cd-hit-est blastn prodigal hmmsearch genomad; do
  command -v "$t" >/dev/null 2>&1 && echo "  OK   $t" || echo "  --   $t (absent on PATH)"
done

echo; echo "=== partitions ==="
sinfo -s 2>/dev/null | head -20
echo; echo "=== preflight done ==="
