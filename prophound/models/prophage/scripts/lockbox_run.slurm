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
# The first lockbox (analysis/lockbox/): run the frozen prophound on the locked test set. GPU. Run ONCE.
# Each mode is a separate run of the shipped CLI (`--mode <m>`), so regions.tsv is
# exactly what a user gets: seed vote, short-contig scan, edge snap and N-gap
# filter included, with nothing rebuilt offline. Default also writes per-checkpoint
# scores, for diagnostics only. Scoring is local and separate, also run once:
# scripts/lockbox_score.py. The preflight refuses to run anything but the frozen
# code (amendment 2026-10-16), the C8 checkpoints and the
# preregistered genome list.
# Submit from /work/$USER/mge_explorer/data_refactored:
#   sbatch scripts/lockbox_run.slurm
#SBATCH -J lockbox_run
#SBATCH -p gpu1v100
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -t 1-00:00:00
#SBATCH -o analysis/_slurm_logs/%x_%j.out
#SBATCH -e analysis/_slurm_logs/%x_%j.err

set -uo pipefail
echo "=== job start ==="; date; hostname; nvidia-smi --query-gpu=name --format=csv,noheader

module load miniconda
conda activate /home/$USER/.conda/envs/mge_explorer_gpu
cd /work/$USER/mge_explorer/data_refactored

LB=analysis/lockbox
LIST=$LB/lockbox_genome_list.txt
FNA=../raw/lockbox
OUT=$LB/prophound
MODES="default exploratory conservative"
LIST_SHA=34f21b35b35188403c75a93c45939f5a0ae259bc31c9779a02b23ac3767d0eb7   # the preregistered list
mkdir -p analysis/_slurm_logs "$OUT"

declare -A SHA=(   # docs/research-integrity.md C8
  [r1]=3259673877ced18e79a92227f898b2009f69cfe82b692c59cb355d74dc27655f
  [r2]=7a50701a0be0f9cca87f000b428341e67d55242f77807321b3237ee73f31e59c
  [r3]=541cd72f1508ca35e2829b13eaf6818727db9cc53e55399c2ad467bba1199484
)
declare -A MD5=(   # frozen at ad5b28b; unchanged at e046579
  [src/predict_genome.py]=62a40e2c38c0a22dbb5d8fa45e82b492
  [src/io_checks.py]=6d4012633a75c621012b41652f150e9a
  [src/multiscale.py]=a6a045651ab0aad1dbaf6b41a15865c0
  [src/train_prophage.py]=2e10a611fbcec78ded2b6a825035762f
  [src/model_rc_cnn.py]=cb95b3c89abc59426480efa5ec7fb720
)

# ---- preflight -------------------------------------------------------------
fail=0
[ -s "$LIST" ] || { echo "PREFLIGHT FAIL: missing or empty $LIST" >&2; fail=1; }
got=$(sha256sum "$LIST" 2>/dev/null | cut -d' ' -f1)
[ "$got" = "$LIST_SHA" ] || { echo "PREFLIGHT FAIL: $LIST sha256 $got != preregistered $LIST_SHA" >&2; fail=1; }
for f in "${!MD5[@]}"; do
  got=$(md5sum "$f" 2>/dev/null | cut -d' ' -f1)
  [ "$got" = "${MD5[$f]}" ] || { echo "PREFLIGHT FAIL: $f md5 $got != frozen ${MD5[$f]}" >&2; fail=1; }
done
for R in r1 r2 r3; do
  C=runs/prophage_v2_pool_$R/checkpoints/best.pt
  [ -s "$C" ] || { echo "PREFLIGHT FAIL: missing $C" >&2; fail=1; continue; }
  got=$(sha256sum "$C" | cut -d' ' -f1)
  [ "$got" = "${SHA[$R]}" ] || { echo "PREFLIGHT FAIL: $C sha256 $got != recorded ${SHA[$R]}" >&2; fail=1; }
done
while read -r G; do
  [ -s "$FNA/$G.fna" ] || { echo "PREFLIGHT FAIL: missing $FNA/$G.fna" >&2; fail=1; }
done < "$LIST"
[ $fail -eq 0 ] || { echo "PREFLIGHT FAILED -- nothing run" >&2; exit 1; }
echo "[preflight] list, frozen code, checkpoints and all $(wc -l < "$LIST") FASTAs verified"
# ---------------------------------------------------------------------------

CKPT="runs/prophage_v2_pool_r1/checkpoints/best.pt runs/prophage_v2_pool_r2/checkpoints/best.pt runs/prophage_v2_pool_r3/checkpoints/best.pt"
CFG=runs/prophage_v2_pool_r1/config_snapshot.yaml
for m in $MODES; do
  echo "[$m]"; date
  EXTRA=""; [ "$m" = default ] && EXTRA="--write-seed-scores"
  while read -r G; do
    D="$OUT/$m/$G"
    [ -s "$D/run_metadata.json" ] && continue                        # resumable
    s=$(date +%s)
    python -m src.predict_genome --input_fasta "$FNA/$G.fna" --outdir "$D" --config "$CFG" --checkpoint $CKPT \
      --mode "$m" --device cuda $EXTRA >>"$OUT/run.log" 2>&1 \
      && echo -e "$m\t$G\t$(( $(date +%s) - s ))" >> "$OUT/wall_seconds.tsv" \
      || echo "FAILED $m $G" >> "$OUT/failures.txt"
  done < "$LIST"
done
N=$(wc -l < "$LIST")
echo "[done] $(ls "$OUT"/*/*/run_metadata.json 2>/dev/null | wc -l) of $(( N * 3 )) genome runs"; date

rc=0
if [ -s "$OUT/failures.txt" ]; then
  echo "WARNING: $(wc -l < "$OUT/failures.txt") genome runs failed -- see $OUT/failures.txt" >&2
  rc=2
fi
echo "=== job end (exit $rc) ==="; date
exit $rc
