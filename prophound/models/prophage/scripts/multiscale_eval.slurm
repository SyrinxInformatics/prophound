#!/bin/bash
# docs/adaptive-multiscale-detection.md Parts 1+3, docs/downstream-postprocessing.md:
# does the fine refinement pass (split + edge placement) fix conflation and the
# small-element boundary problem without costing recall or undoing gap-filling?
# Scores the held-out size-stratified sample once per (seed | ensemble) x fine
# scale, writing both score tracks, then sweeps pass B's parameters offline.
# Submit from /work/$USER/mge_explorer/data_refactored:
#   sbatch scripts/multiscale_eval.slurm
#SBATCH -J multiscale_eval
#SBATCH -p gpu1v100
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -t 24:00:00
#SBATCH -o analysis/_slurm_logs/%x_%j.out
#SBATCH -e analysis/_slurm_logs/%x_%j.err

set -uo pipefail
echo "=== job start ==="; date; hostname; nvidia-smi --query-gpu=name --format=csv,noheader

module load miniconda
conda activate /home/$USER/.conda/envs/mge_explorer_gpu
cd /work/$USER/mge_explorer/data_refactored

LIST=scripts/eval_genome_list_sizestrat.txt
ELEMENTS=scripts/eval_sizestrat_elements.csv
OUT=analysis/multiscale
GENOMES=../raw/genomes        # sibling of data_refactored, not inside it (docs/environment.md)
SCALES="3000 5000"
mkdir -p analysis/_slurm_logs "$OUT"

# ---- preflight: 816440/816441 "COMPLETED" having done nothing -------------
fail=0
for f in "$LIST" "$ELEMENTS" src/multiscale.py scripts/refine_sweep.py scripts/compare_variants.py; do
  [ -s "$f" ] || { echo "PREFLIGHT FAIL: missing or empty $f" >&2; fail=1; }
done
grep -q "refine-window-len" src/predict_genome.py \
  || { echo "PREFLIGHT FAIL: src/predict_genome.py predates --refine-window-len (not pushed?)" >&2; fail=1; }
[ -d "$GENOMES" ] || { echo "PREFLIGHT FAIL: no genome dir $GENOMES" >&2; fail=1; }
for R in r1 r2 r3; do
  for f in runs/prophage_v2_pool_$R/checkpoints/best.pt runs/prophage_v2_pool_$R/config_snapshot.yaml; do
    [ -s "$f" ] || { echo "PREFLIGHT FAIL: missing $f" >&2; fail=1; }
  done
done
if [ -s "$LIST" ]; then
  miss=0
  while read -r G; do
    [ -z "$G" ] && continue
    [ -s "$GENOMES/$G.fna" ] || { miss=$((miss+1)); [ $miss -le 3 ] && echo "PREFLIGHT: no FASTA $GENOMES/$G.fna" >&2; }
  done < "$LIST"
  [ $miss -eq 0 ] || { echo "PREFLIGHT FAIL: $miss listed genomes have no FASTA" >&2; fail=1; }
  echo "[preflight] $(wc -l < "$LIST") genomes listed"
fi
[ $fail -eq 0 ] || { echo "PREFLIGHT FAILED -- nothing run" >&2; exit 1; }
echo "[preflight] all inputs present"
# ---------------------------------------------------------------------------

CKPT_ALL="runs/prophage_v2_pool_r1/checkpoints/best.pt runs/prophage_v2_pool_r2/checkpoints/best.pt runs/prophage_v2_pool_r3/checkpoints/best.pt"
RUN_ARGS=""
for W in $SCALES; do
  for R in r1 r2 r3 ens; do
    if [ "$R" = ens ]; then CKPT=$CKPT_ALL; CFG=runs/prophage_v2_pool_r1/config_snapshot.yaml
    else CKPT=runs/prophage_v2_pool_$R/checkpoints/best.pt; CFG=runs/prophage_v2_pool_$R/config_snapshot.yaml; fi
    TAG="${R}_${W}"
    while read -r G; do
      [ -z "$G" ] && continue
      D="$OUT/$TAG/$G"
      [ -s "$D/fine_windows.tsv" ] && continue          # resumable
      # --refine-threshold 0.95 only sets predict's own regions.tsv; the sweep
      # rebuilds everything from the saved score tracks at every threshold
      python -m src.predict_genome \
        --input_fasta "$GENOMES/$G.fna" --outdir "$D" \
        --config "$CFG" --checkpoint $CKPT \
        --mode default --device cuda \
        --refine-window-len "$W" --refine-stride 250 --refine-threshold 0.95 \
        >>"$OUT/${TAG}.log" 2>&1 || echo "FAILED $TAG $G" >>"$OUT/failures.txt"
    done < "$LIST"
    echo "[done] $TAG: $(ls "$OUT/$TAG" 2>/dev/null | wc -l) genomes"; date
    RUN_ARGS="$RUN_ARGS --run $TAG=$OUT/$TAG"
  done
done

python scripts/refine_sweep.py --elements "$ELEMENTS" $RUN_ARGS --out "$OUT/sweep.json"
rc=$?
if [ -s "$OUT/failures.txt" ]; then
  echo "WARNING: $(wc -l < "$OUT/failures.txt") genome runs failed -- see $OUT/failures.txt" >&2
  [ $rc -eq 0 ] && rc=2
fi
echo "=== job end (exit $rc) ==="; date
exit $rc
