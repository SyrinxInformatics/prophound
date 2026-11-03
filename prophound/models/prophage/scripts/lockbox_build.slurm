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
# docs/lockbox.md: build the locked test set. CPU only; no model is run here.
#   1. fetch NCBI's current assembly_summary_refseq.txt (kept, with its SHA-256)
#   2. select 150 genus-disjoint complete RefSeq bacterial genomes (seeded; never
#      re-drawn once analysis/lockbox/lockbox_genome_list.txt exists)
#   3. download their FASTAs to ../raw/lockbox/
#   4. geNomad annotate + find-proviruses -> ../genomad_out_lockbox/
#   5. write the truth manifest analysis/lockbox/lockbox_elements.csv
# Scoring with prophound is a separate job, run once, after the defaults freeze.
# Submit from /work/$USER/mge_explorer/data_refactored:
#   sbatch scripts/lockbox_build.slurm
#SBATCH -J lockbox_build
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

OUT=analysis/lockbox
FNA=../raw/lockbox              # sibling of raw/genomes, kept apart from the corpus
GMD=../genomad_out_lockbox
DB=/work/$USER/mge_explorer/genomad_db/genomad_db
PYENV=/home/$USER/.conda/envs/mge_explorer_gpu
GMDENV=/home/$USER/.conda/envs/mge_explorer
N=150                           # ~220 geNomad elements at the corpus's 1.5 per genome (docs/lockbox.md)
SEED=20260929
mkdir -p analysis/_slurm_logs "$OUT" "$FNA" "$GMD"

# ---- preflight (docs/environment.md) ---------------------------------------
fail=0
for f in scripts/lockbox_select.py splits/train.txt splits/val.txt splits/test.txt; do
  [ -s "$f" ] || { echo "PREFLIGHT FAIL: missing or empty $f" >&2; fail=1; }
done
[ -d "$DB" ] || { echo "PREFLIGHT FAIL: no geNomad database at $DB" >&2; fail=1; }
for e in "$PYENV" "$GMDENV"; do [ -d "$e" ] || { echo "PREFLIGHT FAIL: no conda env $e" >&2; fail=1; }; done
conda activate "$GMDENV" && command -v genomad >/dev/null \
  || { echo "PREFLIGHT FAIL: genomad not on PATH in $GMDENV" >&2; fail=1; }
command -v curl >/dev/null || { echo "PREFLIGHT FAIL: no curl" >&2; fail=1; }
curl -fsI https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/ >/dev/null \
  || { echo "PREFLIGHT FAIL: no network access to ftp.ncbi.nlm.nih.gov from $(hostname)" >&2; fail=1; }
[ $fail -eq 0 ] || { echo "PREFLIGHT FAILED -- nothing run" >&2; exit 1; }
echo "[preflight] all inputs present"
# ---------------------------------------------------------------------------

conda activate "$PYENV"
LIST="$OUT/lockbox_genome_list.txt"
if [ ! -s "$LIST" ]; then
  SUM="$OUT/assembly_summary_refseq.$(date +%Y%m%d).txt"
  [ -s "$SUM" ] || curl -fsSL -o "$SUM" https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt \
    || { echo "FAILED: summary download" >&2; exit 1; }
  python scripts/lockbox_select.py select --summary "$SUM" \
    --corpus_lists splits/train.txt,splits/val.txt,splits/test.txt \
    --n $N --seed $SEED --out_dir "$OUT" || exit $?
else
  echo "[select] $LIST exists -- not re-drawing ($(wc -l < "$LIST") genomes)"
fi

echo "[download]"; date
tail -n +2 "$OUT/lockbox_selection.tsv" | while IFS=$'\t' read -r G ACC ORG GEN CAT REL URL; do
  [ -s "$FNA/$G.fna" ] && continue
  curl -fsSL --retry 3 -o "$FNA/$G.fna.gz" "$URL" && gunzip -f "$FNA/$G.fna.gz" \
    || { echo "FAILED download $G $URL" >>"$OUT/failures.txt"; rm -f "$FNA/$G.fna.gz"; }
done
echo "[download] $(ls "$FNA"/*.fna 2>/dev/null | wc -l) of $(wc -l < "$LIST") FASTAs present"

echo "[genomad]"; date
conda activate "$GMDENV"
SCR="/scratch/$USER/${SLURM_JOB_ID:-job}"; [ -d "/scratch/$USER" ] || SCR="/work/$USER/tmp/${SLURM_JOB_ID:-job}"
mkdir -p "$SCR"; export TMPDIR="$SCR" MMSEQS_TMP="$SCR"
THREADS=4; P=$(( ${SLURM_CPUS_PER_TASK:-4} / THREADS )); [ "$P" -ge 1 ] || P=1
export GMD DB THREADS OUT
sed "s|^|$FNA/|; s|$|.fna|" "$LIST" | xargs -n1 -P"$P" -I{} bash -c '
  f="{}"; b=$(basename "$f" .fna); od="$GMD/$b"
  [ -s "$od/${b}_find_proviruses/${b}_provirus.tsv" ] && exit 0
  [ -s "$f" ] || exit 0                                   # download failure already recorded
  mkdir -p "$od"
  genomad annotate --threads "$THREADS" "$f" "$od" "$DB" > "$od/genomad.log" 2>&1 \
    && genomad find-proviruses --threads "$THREADS" "$f" "$od" "$DB" >> "$od/genomad.log" 2>&1 \
    && [ -e "$od/${b}_find_proviruses/${b}_provirus.tsv" ] \
    || echo "FAILED genomad $b" >> "$OUT/failures.txt"
'
echo "[genomad] $(ls "$GMD"/*/*_find_proviruses/*_provirus.tsv 2>/dev/null | wc -l) of $(wc -l < "$LIST") done"; date

conda activate "$PYENV"
python scripts/lockbox_select.py manifest --genomad_dir "$GMD" --genome_list "$LIST" \
  --out "$OUT/lockbox_elements.csv"
rc=$?
if [ -s "$OUT/failures.txt" ]; then
  echo "WARNING: $(wc -l < "$OUT/failures.txt") failures -- see $OUT/failures.txt" >&2
  [ $rc -eq 0 ] && rc=2
fi
echo "=== job end (exit $rc) ==="; date
exit $rc
