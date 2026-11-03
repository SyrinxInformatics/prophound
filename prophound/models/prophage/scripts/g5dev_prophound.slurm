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
# Development data for the fragmented-input fix (docs/claims-and-evidence.md
# "G5b preregistration": the fix is built on corpus shreds only, never on LAMBDA
# or Philympics). GPU. The 42 size-stratified eval genomes (tuning) and the 18
# confirmation genomes, each complete and shredded exactly as for G5
# (shred_genomes.py --dist lognormal --sigma 1.0 --min_contig 1500; MAG N50 20 kb
# seeds 0-2, draft N50 100 kb seed 0), scored with the 3-checkpoint ensemble under
# shipped defaults with --write-seed-scores. Truth is geNomad (the eval element
# tables); only prophound runs, so its retention is measured against itself.
# Shred checksums go to analysis/g5dev/shred_md5.tsv. Submit from
# /work/$USER/mge_explorer/data_refactored:
#   sbatch scripts/g5dev_prophound.slurm
#SBATCH -J g5dev_prophound
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

D=analysis/g5dev
GENOMES=../raw/genomes           # sibling of data_refactored
SHR=../raw/g5dev_shreds
OUT=$D/prophound
LIST=$D/genome_list.txt
mkdir -p analysis/_slurm_logs "$OUT" "$SHR"

declare -A SHA=(
  [r1]=3259673877ced18e79a92227f898b2009f69cfe82b692c59cb355d74dc27655f
  [r2]=7a50701a0be0f9cca87f000b428341e67d55242f77807321b3237ee73f31e59c
  [r3]=541cd72f1508ca35e2829b13eaf6818727db9cc53e55399c2ad467bba1199484
)

# ---- preflight -------------------------------------------------------------
fail=0
for f in scripts/eval_genome_list_sizestrat.txt scripts/eval_genome_list.txt scripts/shred_genomes.py \
         src/predict_genome.py; do
  [ -s "$f" ] || { echo "PREFLIGHT FAIL: missing or empty $f" >&2; fail=1; }
done
grep -q -- '--dist' scripts/shred_genomes.py || { echo "PREFLIGHT FAIL: shred_genomes.py lacks --dist" >&2; fail=1; }
grep -q '"exploratory": 0.90' src/predict_genome.py || { echo "PREFLIGHT FAIL: predict_genome.py lacks the final modes" >&2; fail=1; }
grep -q -- '--write-seed-scores' src/predict_genome.py || { echo "PREFLIGHT FAIL: predict_genome.py lacks --write-seed-scores" >&2; fail=1; }
for R in r1 r2 r3; do
  C=runs/prophage_v2_pool_$R/checkpoints/best.pt
  [ -s "$C" ] || { echo "PREFLIGHT FAIL: missing $C" >&2; fail=1; continue; }
  got=$(sha256sum "$C" | cut -d' ' -f1)
  [ "$got" = "${SHA[$R]}" ] || { echo "PREFLIGHT FAIL: $C sha256 $got != ${SHA[$R]}" >&2; fail=1; }
done
sort -u scripts/eval_genome_list_sizestrat.txt scripts/eval_genome_list.txt | grep . > "$LIST"
miss=0; while read -r g; do [ -s "$GENOMES/$g.fna" ] || { echo "PREFLIGHT FAIL: no $GENOMES/$g.fna" >&2; miss=$((miss + 1)); }; done < "$LIST"
[ $miss -eq 0 ] || fail=1
[ $fail -eq 0 ] || { echo "PREFLIGHT FAILED -- nothing run" >&2; exit 1; }
echo "[preflight] $(wc -l < "$LIST") genomes (42 + 18, deduplicated); checkpoint checksums match"
# ---------------------------------------------------------------------------

shred() {  # condition, then shred_genomes.py arguments
  local c=$1; shift
  [ "$(ls "$SHR/$c"/*.fna 2>/dev/null | wc -l)" -eq "$(wc -l < "$LIST")" ] && return 0
  rm -rf "$SHR/$c.tmp"
  python scripts/shred_genomes.py --genome_list "$LIST" --genome_dir "$GENOMES" --out_dir "$SHR/$c.tmp" \
    --dist lognormal --sigma 1.0 --min_contig 1500 "$@" && rm -rf "$SHR/$c" && mv "$SHR/$c.tmp" "$SHR/$c" \
    || { echo "FAILED shred $c" >&2; exit 1; }
}
shred mag_s0 --n50 20000 --seed 0
shred mag_s1 --n50 20000 --seed 1
shred mag_s2 --n50 20000 --seed 2
shred draft_s0 --n50 100000 --seed 0
( cd "$SHR" && md5sum */*.fna ) | awk '{print $2"\t"$1}' > "$D/shred_md5.tsv"

CKPT="runs/prophage_v2_pool_r1/checkpoints/best.pt runs/prophage_v2_pool_r2/checkpoints/best.pt runs/prophage_v2_pool_r3/checkpoints/best.pt"
CFG=runs/prophage_v2_pool_r1/config_snapshot.yaml
for c in complete mag_s0 mag_s1 mag_s2 draft_s0; do
  src=$([ "$c" = complete ] && echo "$GENOMES" || echo "$SHR/$c")
  echo "[$c]"; date
  while read -r g; do
    d="$OUT/$c/$g"
    [ -s "$d/run_metadata.json" ] && continue                        # resumable
    s=$(date +%s)
    python -m src.predict_genome --input_fasta "$src/$g.fna" --outdir "$d" --config "$CFG" --checkpoint $CKPT \
      --mode default --device cuda --write-seed-scores >>"$OUT/run.log" 2>&1 \
      && echo -e "$c\t$g\t$(( $(date +%s) - s ))" >> "$OUT/wall_seconds.tsv" \
      || echo "FAILED $c $g" >> "$OUT/failures.txt"
  done < "$LIST"
done
n=$(( $(wc -l < "$LIST") * 5 ))
echo "[done] $(ls -d "$OUT"/*/*/ 2>/dev/null | wc -l) of $n genome runs"; date

rc=0
if [ -s "$OUT/failures.txt" ]; then
  echo "WARNING: $(wc -l < "$OUT/failures.txt") genome runs failed -- see $OUT/failures.txt" >&2
  rc=2
fi
echo "=== job end (exit $rc) ==="; date
exit $rc
