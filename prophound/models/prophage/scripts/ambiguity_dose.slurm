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
# docs/release-scope-0.1.0.md A.3: ambiguity dose-response. Spikes IUPAC codes into the
# 42-genome size-stratified sample (substitution only; truth unchanged) and
# scores each condition with the shipped ensemble under n-channel (old) and
# impute (default since 2026-09-24):
#   scattered codes at 10 / 100 / 300 / 1000 / 3000 ppm (EDL933 is 1,182)
#   100 bp N runs at 0.5 / 2 / 10 per Mb (scaffold gaps)
# Rates are nested within a pattern. Then scripts/ambiguity_dose.py eval.
# Submit from /work/$USER/mge_explorer/data_refactored:
#   sbatch scripts/ambiguity_dose.slurm
#SBATCH -J ambiguity_dose
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
GENOMES=../raw/genomes        # sibling of data_refactored/, not inside it
OUT=analysis/ambiguity_dose
TMP=$OUT/_fasta_tmp
CONDITIONS="clean:0 scattered:10 scattered:100 scattered:300 scattered:1000 scattered:3000 nruns:0.5 nruns:2 nruns:10"
POLICIES="n-channel impute"
mkdir -p analysis/_slurm_logs "$OUT"

# ---- preflight (docs/environment.md) ---------------------------------------
fail=0
for f in "$LIST" "$ELEMENTS" scripts/ambiguity_dose.py scripts/mode_recheck.py scripts/compare_variants.py \
         scripts/refine_sweep.py src/predict_genome.py; do
  [ -s "$f" ] || { echo "PREFLIGHT FAIL: missing or empty $f" >&2; fail=1; }
done
grep -q '"impute"' src/predict_genome.py \
  || { echo "PREFLIGHT FAIL: src/predict_genome.py lacks --ambiguity-policy impute (not pushed?)" >&2; fail=1; }
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
    done_all=1
    for P in $POLICIES; do [ -s "$OUT/${P}__${PAT}__${RATE}/$G/run_metadata.json" ] || done_all=0; done
    [ $done_all -eq 1 ] && continue                       # resumable
    if [ "$PAT" = clean ]; then SRC="$GENOMES/$G.fna"
    else
      SRC="$TMP/$G.fna"                                   # keeps the genome's name
      python scripts/ambiguity_dose.py spike --fasta "$GENOMES/$G.fna" --pattern "$PAT" --rate "$RATE" \
        --out "$SRC" >>"$OUT/spike.log" 2>&1 || { echo "FAILED spike $C $G" >>"$OUT/failures.txt"; continue; }
    fi
    for P in $POLICIES; do
      D="$OUT/${P}__${PAT}__${RATE}/$G"
      [ -s "$D/run_metadata.json" ] && continue
      python -m src.predict_genome --input_fasta "$SRC" --outdir "$D" --config "$CFG" --checkpoint $CKPT \
        --mode default --device cuda --ambiguity-policy "$P" \
        >>"$OUT/${P}__${PAT}__${RATE}.log" 2>&1 || echo "FAILED $P $C $G" >>"$OUT/failures.txt"
    done
    [ "$PAT" = clean ] || rm -f "$SRC"
  done < "$LIST"
  echo "[done] $C"; date
done
rmdir "$TMP" 2>/dev/null

python scripts/ambiguity_dose.py eval --elements "$ELEMENTS" --root "$OUT" --out "$OUT/results.json"
rc=$?
if [ -s "$OUT/failures.txt" ]; then
  echo "WARNING: $(wc -l < "$OUT/failures.txt") runs failed -- see $OUT/failures.txt" >&2
  [ $rc -eq 0 ] && rc=2
fi
echo "=== job end (exit $rc) ==="; date
exit $rc
