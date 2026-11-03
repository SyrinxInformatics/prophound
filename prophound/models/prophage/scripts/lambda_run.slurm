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
# LAMBDA clean subset, step 2 of 2. GPU. Scores each fetched sequence with the
# 3-checkpoint ensemble under the shipped defaults, writing per-checkpoint scores
# (--write-seed-scores) so every mode's seed vote is re-applied offline by
# scripts/lambda_score.py. The checkpoints must match the SHA-256s recorded
# 2026-09-29 (docs/research-integrity.md C8), the ones the release will freeze.
# Submit after lambda_fetch.slurm, from /work/$USER/mge_explorer/data_refactored:
#   sbatch --dependency=afterok:<fetch job> scripts/lambda_run.slurm
#SBATCH -J lambda_run
#SBATCH -p gpu1v100
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -t 08:00:00
#SBATCH -o analysis/_slurm_logs/%x_%j.out
#SBATCH -e analysis/_slurm_logs/%x_%j.err

set -uo pipefail
echo "=== job start ==="; date; hostname; nvidia-smi --query-gpu=name --format=csv,noheader

module load miniconda
conda activate /home/$USER/.conda/envs/mge_explorer_gpu
cd /work/$USER/mge_explorer/data_refactored

L=analysis/external_benchmarks/lambda
FNA=../raw/lambda
OUT=analysis/lambda/ens
mkdir -p analysis/_slurm_logs "$OUT"

declare -A SHA=(
  [r1]=3259673877ced18e79a92227f898b2009f69cfe82b692c59cb355d74dc27655f
  [r2]=7a50701a0be0f9cca87f000b428341e67d55242f77807321b3237ee73f31e59c
  [r3]=541cd72f1508ca35e2829b13eaf6818727db9cc53e55399c2ad467bba1199484
)

# ---- preflight -------------------------------------------------------------
fail=0
for f in "$L/lambda_genomes.tsv" src/predict_genome.py src/io_checks.py; do
  [ -s "$f" ] || { echo "PREFLIGHT FAIL: missing or empty $f" >&2; fail=1; }
done
grep -q '"exploratory": 0.90' src/predict_genome.py \
  || { echo "PREFLIGHT FAIL: src/predict_genome.py lacks the final modes" >&2; fail=1; }
grep -q -- '--write-seed-scores' src/predict_genome.py \
  || { echo "PREFLIGHT FAIL: src/predict_genome.py lacks --write-seed-scores" >&2; fail=1; }
for R in r1 r2 r3; do
  C=runs/prophage_v2_pool_$R/checkpoints/best.pt
  [ -s "$C" ] || { echo "PREFLIGHT FAIL: missing $C" >&2; fail=1; continue; }
  got=$(sha256sum "$C" | cut -d' ' -f1)
  [ "$got" = "${SHA[$R]}" ] || { echo "PREFLIGHT FAIL: $C sha256 $got != recorded ${SHA[$R]}" >&2; fail=1; }
done
n=$(ls "$FNA"/*.fna 2>/dev/null | wc -l)
[ "$n" -eq 39 ] || { echo "PREFLIGHT FAIL: $FNA has $n/39 sequences (run lambda_fetch.slurm)" >&2; fail=1; }
[ $fail -eq 0 ] || { echo "PREFLIGHT FAILED -- nothing run" >&2; exit 1; }
echo "[preflight] all inputs present; checkpoint checksums match"
# the offline scorer does not apply the >=50%-N call filter (MAX_CALL_N_FRAC), which is a
# no-op only if the sequences carry no long N runs; record their N content
for f in "$FNA"/*.fna; do
  n=$(grep -v '^>' "$f" | tr -cd 'Nn' | wc -c)
  [ "$n" -gt 0 ] && echo "[N content] $(basename "$f" .fna): $n N"
done > "$OUT/../n_content.txt"
echo "[preflight] sequences with any N: $(wc -l < "$OUT/../n_content.txt") of 39 (analysis/lambda/n_content.txt)"
# ---------------------------------------------------------------------------

CKPT="runs/prophage_v2_pool_r1/checkpoints/best.pt runs/prophage_v2_pool_r2/checkpoints/best.pt runs/prophage_v2_pool_r3/checkpoints/best.pt"
CFG=runs/prophage_v2_pool_r1/config_snapshot.yaml
tail -n +2 "$L/lambda_genomes.tsv" | cut -f1 | while read -r ASM; do
  D="$OUT/$ASM"
  [ -s "$D/run_metadata.json" ] && continue                        # resumable
  python -m src.predict_genome --input_fasta "$FNA/$ASM.fna" --outdir "$D" --config "$CFG" --checkpoint $CKPT \
    --mode default --device cuda --write-seed-scores >>"$OUT/../run.log" 2>&1 \
    || echo "FAILED $ASM" >> "$OUT/../failures.txt"
done
echo "[done] $(ls -d "$OUT"/*/ 2>/dev/null | wc -l) of 39 genomes scored"; date

rc=0
if [ -s "$OUT/../failures.txt" ]; then
  echo "WARNING: $(wc -l < "$OUT/../failures.txt") genome runs failed -- see $OUT/../failures.txt" >&2
  rc=2
fi
echo "=== job end (exit $rc) ==="; date
exit $rc
