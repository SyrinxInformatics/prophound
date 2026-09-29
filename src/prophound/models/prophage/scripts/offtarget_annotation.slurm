#!/bin/bash
# docs/error-analysis.md: what are the off-target calls? Summarises geNomad's
# per-gene annotations inside every call in analysis/offtarget/calls42.csv
# (42-genome sample, ensemble, Default + Exploratory, with seed support) and in
# length-matched random controls. CPU only; reads existing geNomad outputs.
# Submit from /work/$USER/mge_explorer/data_refactored:
#   sbatch scripts/offtarget_annotation.slurm
#SBATCH -J offtarget_annotation
#SBATCH -p compute1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -t 02:00:00
#SBATCH -o analysis/_slurm_logs/%x_%j.out
#SBATCH -e analysis/_slurm_logs/%x_%j.err

set -uo pipefail
echo "=== job start ==="; date; hostname

module load miniconda
conda activate /home/$USER/.conda/envs/mge_explorer_gpu
cd /work/$USER/mge_explorer/data_refactored

CALLS=analysis/offtarget/calls42.csv
ELEMENTS=scripts/eval_sizestrat_elements.csv
GENOMAD=../genomad_out      # siblings of data_refactored/, not inside it
GENOMES=../raw/genomes
OUT=analysis/offtarget/annotation.json
mkdir -p analysis/_slurm_logs analysis/offtarget

# ---- preflight (docs/environment.md) ---------------------------------------
fail=0
for f in "$CALLS" "$ELEMENTS" scripts/offtarget_annotation.py; do
  [ -s "$f" ] || { echo "PREFLIGHT FAIL: missing or empty $f" >&2; fail=1; }
done
for d in "$GENOMAD" "$GENOMES"; do
  [ -d "$d" ] || { echo "PREFLIGHT FAIL: no directory $d" >&2; fail=1; }
done
if [ $fail -eq 0 ]; then
  miss=0
  for G in $(tail -n +2 "$CALLS" | cut -d, -f3 | sort -u); do
    [ -s "$GENOMAD/$G/${G}_annotate/${G}_genes.tsv" ] || { miss=$((miss+1)); [ $miss -le 3 ] && echo "PREFLIGHT: no genes.tsv for $G" >&2; }
    [ -s "$GENOMES/$G.fna" ] || { miss=$((miss+1)); [ $miss -le 3 ] && echo "PREFLIGHT: no FASTA for $G" >&2; }
  done
  [ $miss -eq 0 ] || { echo "PREFLIGHT FAIL: $miss missing genome inputs" >&2; fail=1; }
fi
[ $fail -eq 0 ] || { echo "PREFLIGHT FAILED -- nothing run" >&2; exit 1; }
echo "[preflight] all inputs present"
# ---------------------------------------------------------------------------

python scripts/offtarget_annotation.py --calls "$CALLS" --elements "$ELEMENTS" \
  --genomad_dir "$GENOMAD" --genomes_dir "$GENOMES" --out "$OUT"
rc=$?
echo "=== job end (exit $rc) ==="; date
exit $rc
