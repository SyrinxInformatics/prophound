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
# G5 fragmented-input comparison, step 3 of 3. CPU. geNomad end-to-end, default
# settings, on the 39 complete LAMBDA chromosomes and the four shred conditions
# from g5_shred.slurm (195 runs), in the labelling env: geNomad 1.9.4 with DB 1.9
# (docs/research-integrity.md C8). end-to-end, not the labels' annotate +
# find-proviruses, because on fragments an element that fills its contig is
# reported by the classification steps as a whole-contig virus, not as a provirus.
# The find-proviruses table is kept too, as a secondary. Then collects every call
# into one small table, analysis/g5/genomad_calls.tsv, which is all that needs pulling.
# Submit with the other two steps (g5_shred.slurm has the line).
# Paths come from G5_* environment variables (defaults: G5 on LAMBDA). G5b on Philympics:
#   sbatch --export=ALL,G5_DIR=analysis/g5b,G5_FNA=analysis/g5b/fasta,G5_SHR=../raw/g5b_shreds,G5_OUT=../g5b_<tool> <script>
#SBATCH -J g5_genomad
#SBATCH -p compute1
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -t 2-00:00:00
#SBATCH -o analysis/_slurm_logs/%x_%j.out
#SBATCH -e analysis/_slurm_logs/%x_%j.err

set -uo pipefail
echo "=== job start ==="; date; hostname

module load miniconda
cd /work/$USER/mge_explorer/data_refactored

G=${G5_DIR:-analysis/g5}
FNA=${G5_FNA:-../raw/lambda}
SHR=${G5_SHR:-../raw/g5_shreds}
GMD=${G5_OUT:-../g5_genomad}    # raw outputs stay on ARC
DB=/work/$USER/mge_explorer/genomad_db/genomad_db
PYENV=/home/$USER/.conda/envs/mge_explorer_gpu
GMDENV=/home/$USER/.conda/envs/mge_explorer
CONDS="complete mag_s0 mag_s1 mag_s2 draft_s0"
mkdir -p analysis/_slurm_logs "$GMD"

# ---- preflight -------------------------------------------------------------
fail=0
for f in "$G/genome_list.txt" scripts/g5_score.py; do
  [ -s "$f" ] || { echo "PREFLIGHT FAIL: missing or empty $f" >&2; fail=1; }
done
[ -d "$DB" ] || { echo "PREFLIGHT FAIL: no geNomad database at $DB" >&2; fail=1; }
grep -qx 'v\{0,1\}1\.9' "$DB/version.txt" 2>/dev/null \
  || { echo "PREFLIGHT FAIL: geNomad DB version is '$(cat "$DB/version.txt" 2>/dev/null)', expected 1.9" >&2; fail=1; }
for e in "$PYENV" "$GMDENV"; do [ -d "$e" ] || { echo "PREFLIGHT FAIL: no conda env $e" >&2; fail=1; }; done
conda activate "$GMDENV" && command -v genomad >/dev/null \
  || { echo "PREFLIGHT FAIL: genomad not on PATH in $GMDENV" >&2; fail=1; }
v=$(genomad --version 2>&1 | grep -o '[0-9]\+\.[0-9]\+\.[0-9]\+' | head -1)
[ "$v" = "1.9.4" ] || { echo "PREFLIGHT FAIL: genomad --version gave '$(genomad --version 2>&1)', expected 1.9.4" >&2; fail=1; }
N=$(grep -c . "$G/genome_list.txt")
n=$(ls "$FNA"/*.fna 2>/dev/null | wc -l); [ "$n" -eq "$N" ] || { echo "PREFLIGHT FAIL: $FNA has $n/$N" >&2; fail=1; }
for c in mag_s0 mag_s1 mag_s2 draft_s0; do
  n=$(ls "$SHR/$c"/*.fna 2>/dev/null | wc -l)
  [ "$n" -eq "$N" ] || { echo "PREFLIGHT FAIL: $SHR/$c has $n/$N shredded genomes (g5_shred.slurm)" >&2; fail=1; }
done
[ $fail -eq 0 ] || { echo "PREFLIGHT FAILED -- nothing run" >&2; exit 1; }
echo "[preflight] geNomad $v, DB 1.9, all inputs present"
# ---------------------------------------------------------------------------

SCR="/scratch/$USER/${SLURM_JOB_ID:-job}"; [ -d "/scratch/$USER" ] || SCR="/work/$USER/tmp/${SLURM_JOB_ID:-job}"
mkdir -p "$SCR"; export TMPDIR="$SCR" MMSEQS_TMP="$SCR"
THREADS=4; P=$(( ${SLURM_CPUS_PER_TASK:-4} / THREADS )); [ "$P" -ge 1 ] || P=1
export GMD DB THREADS G
for c in $CONDS; do
  echo "[genomad $c]"; date
  src=$([ "$c" = complete ] && echo "$FNA" || echo "$SHR/$c")
  sed "s|^|$src/|; s|$|.fna|" "$G/genome_list.txt" | C=$c xargs -n1 -P"$P" -I{} bash -c '
    f="{}"; b=$(basename "$f" .fna); od="$GMD/$C/$b"
    [ -s "$od/${b}_summary/${b}_virus_summary.tsv" ] && exit 0      # resumable
    mkdir -p "$od"; s=$(date +%s)
    genomad end-to-end --cleanup --threads "$THREADS" "$f" "$od" "$DB" > "$od/genomad.log" 2>&1 \
      && [ -e "$od/${b}_summary/${b}_virus_summary.tsv" ] \
      && echo -e "$C\t$b\t$(( $(date +%s) - s ))" >> "$GMD/wall_seconds.tsv" \
      || echo "FAILED genomad $C $b" >> "$GMD/failures.txt"
  '
  echo "[genomad $c] $(ls "$GMD/$c"/*/*_summary/*_virus_summary.tsv 2>/dev/null | wc -l) of $N"; date
done

conda activate "$PYENV"
python scripts/g5_score.py collect-genomad --root "$GMD" --out "$G/genomad_calls.tsv"
rc=$?
cp "$GMD/wall_seconds.tsv" "$G/genomad_wall_seconds.tsv" 2>/dev/null
if [ -s "$GMD/failures.txt" ]; then
  cp "$GMD/failures.txt" "$G/genomad_failures.txt"
  echo "WARNING: $(wc -l < "$GMD/failures.txt") geNomad runs failed -- see $G/genomad_failures.txt" >&2
  [ $rc -eq 0 ] && rc=2
fi
echo "=== job end (exit $rc) ==="; date
exit $rc
