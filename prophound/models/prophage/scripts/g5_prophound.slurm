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
# G5 fragmented-input comparison, step 2 of 3. GPU. Scores the four shred
# conditions from g5_shred.slurm with the 3-checkpoint ensemble under the shipped
# defaults (--short-contigs scan, impute), writing per-checkpoint scores so every
# mode's seed vote is re-applied offline by scripts/g5_score.py, exactly as for
# LAMBDA (lambda_run.slurm). The complete-genome condition reuses 822275's tracks
# (analysis/lambda/ens): same checkpoints, same defaults.
# Submit with the other two steps (g5_shred.slurm has the line).
# Paths come from G5_* environment variables (defaults: G5 on LAMBDA). G5b runs the
# complete genomes here too (no earlier run exists), from the frozen fix (ad5b28b):
#   sbatch --export=ALL,G5_DIR=analysis/g5b,G5_SHR=../raw/g5b_shreds,G5_FNA=analysis/g5b/fasta scripts/g5_prophound.slurm
#SBATCH -J g5_prophound
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

G=${G5_DIR:-analysis/g5}
SHR=${G5_SHR:-../raw/g5_shreds}
FNA=${G5_FNA:-}                    # set: also run the complete genomes, from $FNA/<genome>.fna
OUT=$G/prophound
CONDS="mag_s0 mag_s1 mag_s2 draft_s0"
[ -n "$FNA" ] && CONDS="complete $CONDS"
N=$(grep -c . "$G/genome_list.txt")
mkdir -p analysis/_slurm_logs "$OUT"

declare -A SHA=(
  [r1]=3259673877ced18e79a92227f898b2009f69cfe82b692c59cb355d74dc27655f
  [r2]=7a50701a0be0f9cca87f000b428341e67d55242f77807321b3237ee73f31e59c
  [r3]=541cd72f1508ca35e2829b13eaf6818727db9cc53e55399c2ad467bba1199484
)

# ---- preflight -------------------------------------------------------------
fail=0
for f in "$G/genome_list.txt" src/predict_genome.py src/io_checks.py; do
  [ -s "$f" ] || { echo "PREFLIGHT FAIL: missing or empty $f" >&2; fail=1; }
done
grep -q '"exploratory": 0.90' src/predict_genome.py \
  || { echo "PREFLIGHT FAIL: src/predict_genome.py lacks the final modes" >&2; fail=1; }
grep -q -- '--write-seed-scores' src/predict_genome.py \
  || { echo "PREFLIGHT FAIL: src/predict_genome.py lacks --write-seed-scores" >&2; fail=1; }
if [ "$G" != analysis/g5 ]; then   # G5b and later: the frozen fragmented-input fix must be on ARC
  grep -q 'DEFAULT_EDGE_SNAP_FRAC = 0.5' src/predict_genome.py \
    || { echo "PREFLIGHT FAIL: src/predict_genome.py lacks the frozen edge snap (ad5b28b)" >&2; fail=1; }
fi
for R in r1 r2 r3; do
  C=runs/prophage_v2_pool_$R/checkpoints/best.pt
  [ -s "$C" ] || { echo "PREFLIGHT FAIL: missing $C" >&2; fail=1; continue; }
  got=$(sha256sum "$C" | cut -d' ' -f1)
  [ "$got" = "${SHA[$R]}" ] || { echo "PREFLIGHT FAIL: $C sha256 $got != recorded ${SHA[$R]}" >&2; fail=1; }
done
for c in $CONDS; do
  d=$SHR/$c; [ "$c" = complete ] && d=$FNA
  while read -r ASM; do
    [ -s "$d/$ASM.fna" ] || { echo "PREFLIGHT FAIL: missing $d/$ASM.fna" >&2; fail=1; }
  done < "$G/genome_list.txt"
done
[ $fail -eq 0 ] || { echo "PREFLIGHT FAILED -- nothing run" >&2; exit 1; }
echo "[preflight] all inputs present; checkpoint checksums match"
# ---------------------------------------------------------------------------

CKPT="runs/prophage_v2_pool_r1/checkpoints/best.pt runs/prophage_v2_pool_r2/checkpoints/best.pt runs/prophage_v2_pool_r3/checkpoints/best.pt"
CFG=runs/prophage_v2_pool_r1/config_snapshot.yaml
for c in $CONDS; do
  echo "[$c]"; date
  while read -r ASM; do
    D="$OUT/$c/$ASM"
    [ -s "$D/run_metadata.json" ] && continue                        # resumable
    s=$(date +%s)
    IN=$SHR/$c/$ASM.fna; [ "$c" = complete ] && IN=$FNA/$ASM.fna
    python -m src.predict_genome --input_fasta "$IN" --outdir "$D" --config "$CFG" --checkpoint $CKPT \
      --mode default --device cuda --write-seed-scores >>"$OUT/run.log" 2>&1 \
      && echo -e "$c\t$ASM\t$(( $(date +%s) - s ))" >> "$OUT/wall_seconds.tsv" \
      || echo "FAILED $c $ASM" >> "$OUT/failures.txt"
  done < "$G/genome_list.txt"
done
echo "[done] $(ls -d "$OUT"/*/*/ 2>/dev/null | wc -l) of $(( N * $(echo $CONDS | wc -w) )) genome runs"; date

rc=0
if [ -s "$OUT/failures.txt" ]; then
  echo "WARNING: $(wc -l < "$OUT/failures.txt") genome runs failed -- see $OUT/failures.txt" >&2
  rc=2
fi
echo "=== job end (exit $rc) ==="; date
exit $rc
