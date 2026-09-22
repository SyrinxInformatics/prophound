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
OUT=analysis/grid_phase
mkdir -p analysis/_slurm_logs "$OUT"

# Five phases across the 1000 bp stride. Uses the three-checkpoint ensemble
# because that is the adopted default -- the question is whether the *shipped*
# configuration is phase-stable, not whether a single seed is.
CKPTS="runs/prophage_v2_pool_r1/checkpoints/best.pt runs/prophage_v2_pool_r2/checkpoints/best.pt runs/prophage_v2_pool_r3/checkpoints/best.pt"
CFG=runs/prophage_v2_pool_r1/config_snapshot.yaml

for OFF in 0 200 400 600 800; do
  while read -r G; do
    [ -z "$G" ] && continue
    D="$OUT/off$OFF/$G"
    [ -s "$D/regions.tsv" ] && continue
    python -m src.predict_genome \
      --input_fasta raw/genomes/$G.fna --outdir "$D" \
      --config "$CFG" --checkpoint $CKPTS --ensemble max \
      --mode default --device cuda \
      --grid-offset "$OFF" \
      >>"$OUT/off$OFF.log" 2>&1 || echo "FAILED off$OFF $G" >>"$OUT/failures.txt"
  done < "$LIST"
  echo "[done] offset $OFF: $(ls "$OUT/off$OFF" 2>/dev/null | wc -l) genomes"; date
done

python scripts/compare_variants.py \
  --elements scripts/eval_sizestrat_elements.csv \
  --variant off0=$OUT/off0 --variant off200=$OUT/off200 --variant off400=$OUT/off400 \
  --variant off600=$OUT/off600 --variant off800=$OUT/off800 \
  --group_by_phase --out $OUT/results.json

echo "=== job end ==="; date
