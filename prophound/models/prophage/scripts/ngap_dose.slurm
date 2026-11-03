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
# Long N gaps (scaffold gaps) under the shipped defaults: impute, and the
# 2026-09-29 seed-vote modes. Extends the earlier ambiguity dose-response, which tested only
# 100 bp N runs. Spikes the 42-genome size-stratified sample with gaps of
# 1 / 5 / 20 kb -- 2 per Mb at random plus one inside every element (capped at
# 25% of its length) -- and scores the ensemble with --write-seed-scores so each
# mode's vote is re-applied offline. Reports recall lost against the clean run,
# and calls lying mostly inside random gaps (calls on imputed sequence).
# Submit from /work/$USER/mge_explorer/data_refactored:
#   sbatch scripts/ngap_dose.slurm
#SBATCH -J ngap_dose
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
GENOMES=../raw/genomes        # sibling of data_refactored/, not inside it
OUT=analysis/ngap_dose
TMP=$OUT/_fasta_tmp
CONDITIONS="clean:0 ngap:1000 ngap:5000 ngap:20000"
P=impute
mkdir -p analysis/_slurm_logs "$OUT"

# ---- preflight -------------------------------------------------------------
fail=0
for f in "$LIST" "$ELEMENTS" scripts/ambiguity_dose.py scripts/mode_recheck.py scripts/compare_variants.py \
         scripts/refine_sweep.py src/predict_genome.py; do
  [ -s "$f" ] || { echo "PREFLIGHT FAIL: missing or empty $f" >&2; fail=1; }
done
grep -q -- '--write-seed-scores' src/predict_genome.py \
  || { echo "PREFLIGHT FAIL: src/predict_genome.py lacks --write-seed-scores (not pushed?)" >&2; fail=1; }
grep -q '"ngap"' scripts/ambiguity_dose.py \
  || { echo "PREFLIGHT FAIL: scripts/ambiguity_dose.py lacks the ngap pattern (not pushed?)" >&2; fail=1; }
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
  echo "[preflight] $(wc -l < "$LIST") genomes"
fi
[ $fail -eq 0 ] || { echo "PREFLIGHT FAILED -- nothing run" >&2; exit 1; }
echo "[preflight] all inputs present"
# ---------------------------------------------------------------------------

CKPT="runs/prophage_v2_pool_r1/checkpoints/best.pt runs/prophage_v2_pool_r2/checkpoints/best.pt runs/prophage_v2_pool_r3/checkpoints/best.pt"
CFG=runs/prophage_v2_pool_r1/config_snapshot.yaml
for C in $CONDITIONS; do
  PAT=${C%%:*}; RATE=${C##*:}
  mkdir -p "$TMP"
  while read -r G; do
    [ -z "$G" ] && continue
    D="$OUT/${P}__${PAT}__${RATE}/$G"
    [ -s "$D/run_metadata.json" ] && continue            # resumable
    if [ "$PAT" = clean ]; then SRC="$GENOMES/$G.fna"
    else
      SRC="$TMP/$G.fna"                                  # keeps the genome's name
      python scripts/ambiguity_dose.py spike --fasta "$GENOMES/$G.fna" --pattern "$PAT" --rate "$RATE" \
        --elements "$ELEMENTS" --gaps_out "$OUT/gaps/${PAT}__${RATE}/$G.tsv" \
        --out "$SRC" >>"$OUT/spike.log" 2>&1 || { echo "FAILED spike $C $G" >>"$OUT/failures.txt"; continue; }
    fi
    # --mode default: regions.tsv is the shipped Default output; eval rebuilds every mode
    python -m src.predict_genome --input_fasta "$SRC" --outdir "$D" --config "$CFG" --checkpoint $CKPT \
      --mode default --device cuda --ambiguity-policy "$P" --write-seed-scores \
      >>"$OUT/${P}__${PAT}__${RATE}.log" 2>&1 || echo "FAILED $P $C $G" >>"$OUT/failures.txt"
    [ "$PAT" = clean ] || rm -f "$SRC"
  done < "$LIST"
  echo "[done] $C: $(ls "$OUT/${P}__${PAT}__${RATE}" 2>/dev/null | wc -l) genomes"; date
done
rmdir "$TMP" 2>/dev/null

python scripts/ambiguity_dose.py eval --elements "$ELEMENTS" --root "$OUT" --gaps_dir "$OUT/gaps" \
  --out "$OUT/results.json"
rc=$?
if [ -s "$OUT/failures.txt" ]; then
  echo "WARNING: $(wc -l < "$OUT/failures.txt") runs failed -- see $OUT/failures.txt" >&2
  [ $rc -eq 0 ] && rc=2
fi
echo "=== job end (exit $rc) ==="; date
exit $rc
