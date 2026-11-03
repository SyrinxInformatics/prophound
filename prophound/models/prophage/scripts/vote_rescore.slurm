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
# Re-score existing ARC tracks under the 2026-09-29 modes (seed vote;
# predict_genome.MODE_THRESHOLDS / MODE_MIN_SEED_SUPPORT). CPU only, no scoring:
#   1. short-contig thresholds re-derived under the vote (818235's fragment tracks
#      against the 42-genome 10 kb reference) -> short_contig_calib/vote_sweep.json
#   2. literature pilot, shipped modes (821391's tracks) -> litpilot/results_shipped.json
#   3. sub-10 kb spike-in, shipped modes, refinement off (821410's tracks)
#      -> spikein_small/results_shipped.json
# Earlier results files are left untouched. Submit from
# /work/$USER/mge_explorer/data_refactored:
#   sbatch scripts/vote_rescore.slurm
#SBATCH -J vote_rescore
#SBATCH -p compute1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -t 06:00:00
#SBATCH -o analysis/_slurm_logs/%x_%j.out
#SBATCH -e analysis/_slurm_logs/%x_%j.err

set -uo pipefail
echo "=== job start ==="; date; hostname

module load miniconda
conda activate /home/$USER/.conda/envs/mge_explorer_gpu
cd /work/$USER/mge_explorer/data_refactored

GENOMES=../raw/genomes        # siblings of data_refactored/, not inside it
GENOMAD=../genomad_out
SC=analysis/short_contig_calib
LP=analysis/litpilot
SP=analysis/spikein_small
REF=analysis/multiscale
mkdir -p analysis/_slurm_logs

# ---- preflight -------------------------------------------------------------
fail=0
for f in scripts/eval_sizestrat_elements.csv scripts/literature_truth.tsv scripts/mode_recheck.py \
         scripts/short_contig_vote_sweep.py scripts/short_contig_sweep.py scripts/shred_genomes.py \
         scripts/literature_pilot_eval.py scripts/small_element_spikein.py scripts/refine_sweep.py \
         scripts/compare_variants.py src/predict_genome.py "$SP/spike_truth.tsv" "$SP/insitu_truth.tsv"; do
  [ -s "$f" ] || { echo "PREFLIGHT FAIL: missing or empty $f" >&2; fail=1; }
done
grep -q MODE_MIN_SEED_SUPPORT src/predict_genome.py \
  || { echo "PREFLIGHT FAIL: src/predict_genome.py predates the seed vote (not pushed?)" >&2; fail=1; }
grep -q shipped_calls scripts/literature_pilot_eval.py \
  || { echo "PREFLIGHT FAIL: scripts/literature_pilot_eval.py predates shipped modes (not pushed?)" >&2; fail=1; }
grep -q -- '--refine' scripts/small_element_spikein.py \
  || { echo "PREFLIGHT FAIL: scripts/small_element_spikein.py predates --refine (not pushed?)" >&2; fail=1; }
for d in "$GENOMES" "$GENOMAD" "$SC/fasta"; do [ -d "$d" ] || { echo "PREFLIGHT FAIL: no directory $d" >&2; fail=1; }; done
for R in r1 r2 r3 ens; do
  n=$(ls "$REF/${R}_3000"/*/windows.tsv 2>/dev/null | wc -l)
  [ "$n" -eq 42 ] || { echo "PREFLIGHT FAIL: $REF/${R}_3000 has $n/42 genomes" >&2; fail=1; }
  n=$(ls "$SC/$R"/*/windows.tsv 2>/dev/null | wc -l)
  [ "$n" -eq 42 ] || { echo "PREFLIGHT FAIL: $SC/$R has $n/42 fragment runs" >&2; fail=1; }
  n=$(ls "$LP/$R"/*/windows.tsv 2>/dev/null | wc -l)
  [ "$n" -eq 7 ] || { echo "PREFLIGHT FAIL: $LP/$R has $n/7 genomes" >&2; fail=1; }
  n=$(ls "$SP/tracks_$R"/*/windows.tsv 2>/dev/null | wc -l)
  [ "$n" -eq 62 ] || { echo "PREFLIGHT FAIL: $SP/tracks_$R has $n/62 genomes" >&2; fail=1; }
done
[ $fail -eq 0 ] || { echo "PREFLIGHT FAILED -- nothing run" >&2; exit 1; }
echo "[preflight] all inputs present"
# ---------------------------------------------------------------------------

rc=0
echo "[1/3] short-contig thresholds under the vote"; date
python scripts/short_contig_vote_sweep.py --elements scripts/eval_sizestrat_elements.csv \
  --fasta_dir "$SC/fasta" --shred_root "$SC" --ref_root "$REF" --ref_suffix _3000 \
  --out "$SC/vote_sweep.json" || { r=$?; echo "FAILED step 1 (exit $r)" >&2; rc=$r; }

echo "[2/3] literature pilot, shipped modes"; date
python scripts/literature_pilot_eval.py --truth scripts/literature_truth.tsv --genomes_dir "$GENOMES" \
  --genomad_dir "$GENOMAD" --run r1=$LP/r1 --run r2=$LP/r2 --run r3=$LP/r3 --run ens=$LP/ens \
  --out "$LP/results_shipped.json" || { r=$?; echo "FAILED step 2 (exit $r)" >&2; rc=$r; }

echo "[3/3] sub-10 kb spike-in, shipped modes"; date
python scripts/small_element_spikein.py eval --spike_truth "$SP/spike_truth.tsv" \
  --insitu_truth "$SP/insitu_truth.tsv" --refine off \
  --run r1="$SP/tracks_r1" --run r2="$SP/tracks_r2" --run r3="$SP/tracks_r3" --run ens="$SP/tracks_ens" \
  --out "$SP/results_shipped.json" || { r=$?; echo "FAILED step 3 (exit $r)" >&2; rc=$r; }

echo "=== job end (exit $rc) ==="; date
exit $rc
