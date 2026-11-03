#!/bin/bash
# SPDX-License-Identifier: Apache-2.0
#
# Copyright 2025 Jacob Alford <jalford0000@gmail.com>
#                and Benjamin Alford <ssolar5.dev@gmail.com>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# docs/adaptive-multiscale-detection.md Part 4: calibrate --short-contigs scan.
# Shreds the 42 held-out genomes into 3-10 kb contigs (a simulated fragmented
# assembly), scores them at 5 kb / 3 kb windows with each seed and the ensemble,
# then sweeps thresholds offline against the intact-genome 10 kb pipeline's
# false-positive load (817370's runs under analysis/multiscale/).
# Submit from /work/$USER/mge_explorer/data_refactored:
#   sbatch scripts/short_contig_calib.slurm
#SBATCH -J short_contig_calib
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
GENOMES=../raw/genomes        # sibling of data_refactored, not inside it (docs/environment.md)
OUT=analysis/short_contig_calib
REF=analysis/multiscale
mkdir -p analysis/_slurm_logs "$OUT"

# ---- preflight --------------------------------------------------------------
fail=0
for f in "$LIST" "$ELEMENTS" scripts/shred_genomes.py scripts/short_contig_sweep.py \
         scripts/refine_sweep.py scripts/compare_variants.py src/multiscale.py; do
  [ -s "$f" ] || { echo "PREFLIGHT FAIL: missing or empty $f" >&2; fail=1; }
done
grep -q "short-contigs" src/predict_genome.py \
  || { echo "PREFLIGHT FAIL: src/predict_genome.py predates --short-contigs (not pushed?)" >&2; fail=1; }
[ -d "$GENOMES" ] || { echo "PREFLIGHT FAIL: no genome dir $GENOMES" >&2; fail=1; }
for R in r1 r2 r3; do
  for f in runs/prophage_v2_pool_$R/checkpoints/best.pt runs/prophage_v2_pool_$R/config_snapshot.yaml; do
    [ -s "$f" ] || { echo "PREFLIGHT FAIL: missing $f" >&2; fail=1; }
  done
done
for R in r1 r2 r3 ens; do
  n=$(ls "$REF/${R}_3000"/*/windows.tsv 2>/dev/null | wc -l)
  [ "$n" -eq 42 ] || { echo "PREFLIGHT FAIL: reference $REF/${R}_3000 has $n/42 genomes" >&2; fail=1; }
done
if [ -s "$LIST" ]; then
  miss=0
  while read -r G; do
    [ -z "$G" ] && continue
    [ -s "$GENOMES/$G.fna" ] || { miss=$((miss+1)); [ $miss -le 3 ] && echo "PREFLIGHT: no FASTA $GENOMES/$G.fna" >&2; }
  done < "$LIST"
  [ $miss -eq 0 ] || { echo "PREFLIGHT FAIL: $miss listed genomes have no FASTA" >&2; fail=1; }
fi
[ $fail -eq 0 ] || { echo "PREFLIGHT FAILED -- nothing run" >&2; exit 1; }
echo "[preflight] all inputs present"
# ---------------------------------------------------------------------------

python scripts/shred_genomes.py --genome_list "$LIST" --genome_dir "$GENOMES" --out_dir "$OUT/fasta" \
  || { echo "shredding failed" >&2; exit 1; }

CKPT_ALL="runs/prophage_v2_pool_r1/checkpoints/best.pt runs/prophage_v2_pool_r2/checkpoints/best.pt runs/prophage_v2_pool_r3/checkpoints/best.pt"
SWEEP_ARGS=""
for R in r1 r2 r3 ens; do
  if [ "$R" = ens ]; then CKPT=$CKPT_ALL; CFG=runs/prophage_v2_pool_r1/config_snapshot.yaml
  else CKPT=runs/prophage_v2_pool_$R/checkpoints/best.pt; CFG=runs/prophage_v2_pool_$R/config_snapshot.yaml; fi
  while read -r G; do
    [ -z "$G" ] && continue
    D="$OUT/$R/$G"
    [ -s "$D/windows.tsv" ] && continue          # resumable
    # thresholds here only shape predict's own regions.tsv; the sweep re-applies
    # every threshold offline from windows.tsv
    python -m src.predict_genome \
      --input_fasta "$OUT/fasta/$G.fna" --outdir "$D" \
      --config "$CFG" --checkpoint $CKPT \
      --mode default --device cuda --short-contigs scan \
      >>"$OUT/${R}.log" 2>&1 || echo "FAILED $R $G" >>"$OUT/failures.txt"
  done < "$LIST"
  echo "[done] $R: $(ls "$OUT/$R" 2>/dev/null | wc -l) genomes"; date
  SWEEP_ARGS="$SWEEP_ARGS --run $R=$OUT/$R --reference $R=$REF/${R}_3000"
done

python scripts/short_contig_sweep.py --elements "$ELEMENTS" --fasta_dir "$OUT/fasta" \
  $SWEEP_ARGS --out "$OUT/sweep.json"
rc=$?
if [ -s "$OUT/failures.txt" ]; then
  echo "WARNING: $(wc -l < "$OUT/failures.txt") genome runs failed -- see $OUT/failures.txt" >&2
  [ $rc -eq 0 ] && rc=2
fi
echo "=== job end (exit $rc) ==="; date
exit $rc
