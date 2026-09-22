#!/bin/bash
# docs/scientific-validity.md "1a": does fixing the ambiguity-code encoding change
# the headline numbers? Scores the held-out size-stratified sample under all three
# --ambiguity-policy values, three seeds each, and compares.
# 'impute' is the in-distribution fix and the one expected to ship; 'uniform' is
# kept as a control (fractional inputs are themselves out of distribution, so if
# the two disagree that is itself the finding).
# Submit from /work/$USER/mge_explorer/data_refactored:
#   sbatch scripts/ambiguity_impact.slurm
#SBATCH -J ambiguity_impact
#SBATCH -p gpu1v100
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -t 16:00:00
#SBATCH -o analysis/_slurm_logs/%x_%j.out
#SBATCH -e analysis/_slurm_logs/%x_%j.err

set -uo pipefail
echo "=== job start ==="; date; hostname; nvidia-smi --query-gpu=name --format=csv,noheader

module load miniconda
conda activate /home/$USER/.conda/envs/mge_explorer_gpu
cd /work/$USER/mge_explorer/data_refactored

LIST=scripts/eval_genome_list_sizestrat.txt
ELEMENTS=scripts/eval_sizestrat_elements.csv
OUT=analysis/ambiguity_impact
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

# default policy first: these runs must reproduce the committed numbers exactly,
# which is the regression check that the new flags changed nothing by default
for POL in n-channel impute uniform; do
  for R in r1 r2 r3; do
    CKPT=runs/prophage_v2_pool_$R/checkpoints/best.pt
    CFG=runs/prophage_v2_pool_$R/config_snapshot.yaml
    TAG="${POL}_${R}"
    while read -r G; do
      [ -z "$G" ] && continue
      D="$OUT/$TAG/$G"
      [ -s "$D/regions.tsv" ] && continue          # resumable: skip finished genomes
      python -m src.predict_genome \
        --input_fasta "$GENOMES/$G.fna" --outdir "$D" \
        --config "$CFG" --checkpoint "$CKPT" \
        --mode default --device cuda \
        --ambiguity-policy "$POL" \
        >>"$OUT/${TAG}.log" 2>&1 || echo "FAILED $TAG $G" >>"$OUT/failures.txt"
    done < "$LIST"
    echo "[done] $TAG: $(ls "$OUT/$TAG" 2>/dev/null | wc -l) genomes"; date
  done
done

ARGS=""
for POL in n-channel impute uniform; do
  for R in r1 r2 r3; do ARGS="$ARGS --variant ${POL}_${R}=$OUT/${POL}_${R}"; done
done
python scripts/compare_variants.py \
  --elements "$ELEMENTS" \
  $ARGS --strata --out $OUT/results.json

echo "=== job end ==="; date
