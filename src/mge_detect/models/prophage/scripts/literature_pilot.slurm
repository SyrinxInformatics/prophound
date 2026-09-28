#!/bin/bash
# docs/final-benchmark.md, release checklist §3: literature-curated pilot.
# Scores the 7 genomes in scripts/literature_truth.tsv (4 test, 3 val) with each
# seed and the ensemble under shipped defaults, then scores every mode, the
# per-genome permutation test and geNomad against the literature truth.
# Submit from /work/$USER/mge_explorer/data_refactored:
#   sbatch scripts/literature_pilot.slurm
#SBATCH -J literature_pilot
#SBATCH -p gpu1v100
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -t 06:00:00
#SBATCH -o analysis/_slurm_logs/%x_%j.out
#SBATCH -e analysis/_slurm_logs/%x_%j.err

set -uo pipefail
echo "=== job start ==="; date; hostname; nvidia-smi --query-gpu=name --format=csv,noheader

module load miniconda
conda activate /home/$USER/.conda/envs/mge_explorer_gpu
cd /work/$USER/mge_explorer/data_refactored

TRUTH=scripts/literature_truth.tsv
GENOMES=../raw/genomes        # siblings of data_refactored/, not inside it
GENOMAD=../genomad_out
OUT=analysis/litpilot
mkdir -p analysis/_slurm_logs "$OUT"

# ---- preflight (docs/environment.md) ---------------------------------------
fail=0
for f in "$TRUTH" scripts/literature_pilot_eval.py scripts/mode_recheck.py scripts/compare_variants.py \
         scripts/refine_sweep.py src/predict_genome.py; do
  [ -s "$f" ] || { echo "PREFLIGHT FAIL: missing or empty $f" >&2; fail=1; }
done
grep -q SHORT_CONTIG_THRESHOLDS src/predict_genome.py \
  || { echo "PREFLIGHT FAIL: src/predict_genome.py predates calibrated short contigs (not pushed?)" >&2; fail=1; }
for d in "$GENOMES" "$GENOMAD"; do [ -d "$d" ] || { echo "PREFLIGHT FAIL: no directory $d" >&2; fail=1; }; done
for R in r1 r2 r3; do
  for f in runs/prophage_v2_pool_$R/checkpoints/best.pt runs/prophage_v2_pool_$R/config_snapshot.yaml; do
    [ -s "$f" ] || { echo "PREFLIGHT FAIL: missing $f" >&2; fail=1; }
  done
done
LIST=$(tail -n +2 "$TRUTH" | cut -f1 | sort -u)
for G in $LIST; do
  [ -s "$GENOMES/$G.fna" ] || { echo "PREFLIGHT FAIL: no FASTA $GENOMES/$G.fna" >&2; fail=1; }
done
echo "[preflight] $(echo $LIST | wc -w) genomes"
[ $fail -eq 0 ] || { echo "PREFLIGHT FAILED -- nothing run" >&2; exit 1; }
echo "[preflight] all inputs present"
# ---------------------------------------------------------------------------

CKPT_ALL="runs/prophage_v2_pool_r1/checkpoints/best.pt runs/prophage_v2_pool_r2/checkpoints/best.pt runs/prophage_v2_pool_r3/checkpoints/best.pt"
RUNS=""
for R in r1 r2 r3 ens; do
  if [ "$R" = ens ]; then CKPT=$CKPT_ALL; CFG=runs/prophage_v2_pool_r1/config_snapshot.yaml
  else CKPT=runs/prophage_v2_pool_$R/checkpoints/best.pt; CFG=runs/prophage_v2_pool_$R/config_snapshot.yaml; fi
  for G in $LIST; do
    D="$OUT/$R/$G"
    [ -s "$D/run_metadata.json" ] && continue          # resumable
    python -m src.predict_genome --input_fasta "$GENOMES/$G.fna" --outdir "$D" \
      --config "$CFG" --checkpoint $CKPT --mode default --device cuda \
      >>"$OUT/${R}.log" 2>&1 || echo "FAILED $R $G" >>"$OUT/failures.txt"
  done
  echo "[done] $R: $(ls "$OUT/$R" 2>/dev/null | wc -l) genomes"; date
  RUNS="$RUNS --run $R=$OUT/$R"
done

python scripts/literature_pilot_eval.py --truth "$TRUTH" --genomes_dir "$GENOMES" --genomad_dir "$GENOMAD" \
  $RUNS --out "$OUT/results.json"
rc=$?
if [ -s "$OUT/failures.txt" ]; then
  echo "WARNING: $(wc -l < "$OUT/failures.txt") genome runs failed -- see $OUT/failures.txt" >&2
  [ $rc -eq 0 ] && rc=2
fi
echo "=== job end (exit $rc) ==="; date
exit $rc
