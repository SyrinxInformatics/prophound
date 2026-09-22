#!/bin/bash
# docs/scientific-validity.md "Input-robustness standing check", GRID-1: how often
# do calls change when the inference window grid shifts? Locally, a 777 bp shift
# turned 3 regions into 5 on one slice; this measures it over the held-out sample.
# Submit from /work/$USER/mge_explorer/data_refactored:
#   sbatch scripts/grid_phase_stability.slurm
#SBATCH -J grid_phase
#SBATCH -p gpu1v100
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -t 12:00:00
#SBATCH -o analysis/_slurm_logs/%x_%j.out
#SBATCH -e analysis/_slurm_logs/%x_%j.err

set -uo pipefail
echo "=== job start ==="; date; hostname; nvidia-smi --query-gpu=name --format=csv,noheader

module load miniconda
conda activate /home/$USER/.conda/envs/mge_explorer_gpu
cd /work/$USER/mge_explorer/data_refactored

LIST=scripts/eval_genome_list_sizestrat.txt
ELEMENTS=scripts/eval_sizestrat_elements.csv
OUT=analysis/grid_phase
CKPT_ALL="runs/prophage_v2_pool_r1/checkpoints/best.pt runs/prophage_v2_pool_r2/checkpoints/best.pt runs/prophage_v2_pool_r3/checkpoints/best.pt"
mkdir -p analysis/_slurm_logs "$OUT"

# ---- preflight: fail loudly before burning a queue slot -------------------
# 816440/816441 (2026-09-22) both "COMPLETED" in 3 seconds having done nothing:
# the genome list had never been pushed, so `while read < $LIST` looped zero
# times, and without this block the job still exited 0. Check every input.
GENOMES=../raw/genomes        # NOT raw/genomes -- it is a sibling of
                              # data_refactored, not inside it (docs/environment.md)
fail=0
for f in "$LIST" "$ELEMENTS"; do
  [ -s "$f" ] || { echo "PREFLIGHT FAIL: missing or empty $f" >&2; fail=1; }
done
[ -d "$GENOMES" ] || { echo "PREFLIGHT FAIL: no genome dir $GENOMES" >&2; fail=1; }
for C in $CKPT_ALL; do
  [ -s "$C" ] || { echo "PREFLIGHT FAIL: missing checkpoint $C" >&2; fail=1; }
done
if [ -s "$LIST" ]; then
  n=$(wc -l < "$LIST"); echo "[preflight] $n genomes listed"
  miss=0
  while read -r G; do
    [ -z "$G" ] && continue
    [ -s "$GENOMES/$G.fna" ] || { miss=$((miss+1)); [ $miss -le 3 ] && echo "PREFLIGHT: no FASTA $GENOMES/$G.fna" >&2; }
  done < "$LIST"
  [ $miss -eq 0 ] || { echo "PREFLIGHT FAIL: $miss listed genomes have no FASTA" >&2; fail=1; }
fi
[ $fail -eq 0 ] || { echo "PREFLIGHT FAILED -- nothing run" >&2; exit 1; }
echo "[preflight] all inputs present"
# --------------------------------------------------------------------------

# Five phases across the 1000 bp stride. Uses the three-checkpoint ensemble
# because that is the adopted default -- the question is whether the *shipped*
# configuration is phase-stable, not whether a single seed is.
CKPTS="$CKPT_ALL"
CFG=runs/prophage_v2_pool_r1/config_snapshot.yaml

for OFF in 0 200 400 600 800; do
  while read -r G; do
    [ -z "$G" ] && continue
    D="$OUT/off$OFF/$G"
    [ -s "$D/regions.tsv" ] && continue
    python -m src.predict_genome \
      --input_fasta "$GENOMES/$G.fna" --outdir "$D" \
      --config "$CFG" --checkpoint $CKPTS --ensemble max \
      --mode default --device cuda \
      --grid-offset "$OFF" \
      >>"$OUT/off$OFF.log" 2>&1 || echo "FAILED off$OFF $G" >>"$OUT/failures.txt"
  done < "$LIST"
  echo "[done] offset $OFF: $(ls "$OUT/off$OFF" 2>/dev/null | wc -l) genomes"; date
done

python scripts/compare_variants.py \
  --elements "$ELEMENTS" \
  --variant off0=$OUT/off0 --variant off200=$OUT/off200 --variant off400=$OUT/off400 \
  --variant off600=$OUT/off600 --variant off800=$OUT/off800 \
  --group_by_phase --out $OUT/results.json

echo "=== job end ==="; date
