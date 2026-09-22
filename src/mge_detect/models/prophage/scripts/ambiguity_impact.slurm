#!/bin/bash
# docs/scientific-validity.md "1a": does fixing the ambiguity-code encoding change
# the headline numbers? Scores the held-out size-stratified sample under both
# --ambiguity-policy values, three seeds each, and compares.
# Submit from /work/$USER/mge_explorer/data_refactored:
#   sbatch scripts/ambiguity_impact.slurm
#SBATCH -J ambiguity_impact
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
OUT=analysis/ambiguity_impact
mkdir -p analysis/_slurm_logs "$OUT"

# default policy first: these runs must reproduce the committed numbers exactly,
# which is the regression check that the new flags changed nothing by default
for POL in n-channel uniform; do
  for R in r1 r2 r3; do
    CKPT=runs/prophage_v2_pool_$R/checkpoints/best.pt
    CFG=runs/prophage_v2_pool_$R/config_snapshot.yaml
    TAG="${POL}_${R}"
    while read -r G; do
      [ -z "$G" ] && continue
      D="$OUT/$TAG/$G"
      [ -s "$D/regions.tsv" ] && continue          # resumable: skip finished genomes
      python -m src.predict_genome \
        --input_fasta raw/genomes/$G.fna --outdir "$D" \
        --config "$CFG" --checkpoint "$CKPT" \
        --mode default --device cuda \
        --ambiguity-policy "$POL" \
        >>"$OUT/${TAG}.log" 2>&1 || echo "FAILED $TAG $G" >>"$OUT/failures.txt"
    done < "$LIST"
    echo "[done] $TAG: $(ls "$OUT/$TAG" 2>/dev/null | wc -l) genomes"; date
  done
done

python scripts/compare_variants.py \
  --elements scripts/eval_sizestrat_elements.csv \
  --variant n-channel_r1=$OUT/n-channel_r1 --variant uniform_r1=$OUT/uniform_r1 \
  --variant n-channel_r2=$OUT/n-channel_r2 --variant uniform_r2=$OUT/uniform_r2 \
  --variant n-channel_r3=$OUT/n-channel_r3 --variant uniform_r3=$OUT/uniform_r3 \
  --strata --out $OUT/results.json

echo "=== job end ==="; date
