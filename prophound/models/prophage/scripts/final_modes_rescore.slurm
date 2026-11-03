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
# Re-score existing ARC tracks under the final 2026-09-29 modes (exploratory 0.90
# keep all; default 0.95 3-of-3; conservative 0.995 3-of-3; SHORT_CONTIG_THRESHOLDS
# 0.95 / 0.95 / 0.99). CPU only, no new scoring. Writes *_final files, so no
# earlier result is overwritten:
#   1. short-contig sweep under the final modes (confirms the thresholds)
#   2. literature pilot, revised truth classes      -> litpilot/results_final.json
#   3. sub-10 kb spike-in (821410 tracks)           -> spikein_small/results_final.json
#   4. spike-in, model-independent sites (822068)   -> spikein_siterule/results_final.json
#   5. long N gaps (822054 tracks)                  -> ngap_dose/results_final.json
#   6. C5 containment, fixed (docs/research-integrity.md) -> integrity/containment.tsv
# Steps 4 and 5 need 822068 and 822054 finished: submit with
#   sbatch --dependency=afterany:822054:822068:822071 scripts/final_modes_rescore.slurm
# from /work/$USER/mge_explorer/data_refactored.
#SBATCH -J final_modes_rescore
#SBATCH -p compute1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -t 12:00:00
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
SR=analysis/spikein_siterule
NG=analysis/ngap_dose
mkdir -p analysis/_slurm_logs analysis/integrity

# ---- preflight -------------------------------------------------------------
fail=0
grep -q '"exploratory": 0.90' src/predict_genome.py \
  || { echo "PREFLIGHT FAIL: src/predict_genome.py lacks the final modes (not pushed?)" >&2; fail=1; }
grep -q phage_derived scripts/literature_truth.tsv \
  || { echo "PREFLIGHT FAIL: scripts/literature_truth.tsv lacks the revised classes" >&2; fail=1; }
grep -q set_genomes_dir scripts/eval_containment.py \
  || { echo "PREFLIGHT FAIL: scripts/eval_containment.py lacks the C5 fix" >&2; fail=1; }
for f in scripts/short_contig_vote_sweep.py scripts/literature_pilot_eval.py scripts/small_element_spikein.py \
         scripts/ambiguity_dose.py scripts/mode_recheck.py scripts/eval_sizestrat_elements.csv \
         analysis/lockbox/lockbox_elements.csv "$SP/spike_truth.tsv" "$SP/insitu_truth.tsv"; do
  [ -s "$f" ] || { echo "PREFLIGHT FAIL: missing or empty $f" >&2; fail=1; }
done
[ $fail -eq 0 ] || { echo "PREFLIGHT FAILED -- nothing run" >&2; exit 1; }
echo "[preflight] required inputs present"
# ---------------------------------------------------------------------------

rc=0
step() { echo "[$1]"; date; }
fail_step() { echo "FAILED $1 (exit $2)" >&2; rc=$2; }

step "1 short-contig sweep"
python scripts/short_contig_vote_sweep.py --elements scripts/eval_sizestrat_elements.csv \
  --fasta_dir "$SC/fasta" --shred_root "$SC" --ref_root analysis/multiscale --ref_suffix _3000 \
  --out "$SC/vote_sweep_final.json" || fail_step 1 $?

step "2 literature pilot"
python scripts/literature_pilot_eval.py --truth scripts/literature_truth.tsv --genomes_dir "$GENOMES" \
  --genomad_dir "$GENOMAD" --run r1=$LP/r1 --run r2=$LP/r2 --run r3=$LP/r3 --run ens=$LP/ens \
  --out "$LP/results_final.json" || fail_step 2 $?

step "3 spike-in"
python scripts/small_element_spikein.py eval --spike_truth "$SP/spike_truth.tsv" --insitu_truth "$SP/insitu_truth.tsv" \
  --refine off --run r1="$SP/tracks_r1" --run r2="$SP/tracks_r2" --run r3="$SP/tracks_r3" --run ens="$SP/tracks_ens" \
  --out "$SP/results_final.json" || fail_step 3 $?

step "4 spike-in, model-independent sites"
if [ -s "$SR/spike_truth.tsv" ] && [ "$(ls "$SR/tracks_ens" 2>/dev/null | wc -l)" -gt 0 ]; then
  python scripts/small_element_spikein.py eval --spike_truth "$SR/spike_truth.tsv" --insitu_truth "$SR/insitu_truth.tsv" \
    --refine off --run r1="$SR/tracks_r1" --run r2="$SR/tracks_r2" --run r3="$SR/tracks_r3" --run ens="$SR/tracks_ens" \
    --out "$SR/results_final.json" || fail_step 4 $?
else
  echo "SKIPPED 4: $SR not built (822068 not finished?)" >&2; rc=2
fi

step "5 long N gaps"
if [ "$(ls -d "$NG"/impute__* 2>/dev/null | wc -l)" -ge 4 ]; then
  python scripts/ambiguity_dose.py eval --elements scripts/eval_sizestrat_elements.csv --root "$NG" \
    --gaps_dir "$NG/gaps" --out "$NG/results_final.json" || fail_step 5 $?
else
  echo "SKIPPED 5: $NG incomplete (822054 not finished?)" >&2; rc=2
fi

step "6 C5 containment"
python scripts/eval_containment.py --beds_dir beds_pos --sequences_dir sequences --train_split splits/train.txt \
  --genomes_dir "$GENOMES" \
  --set sizestrat42=scripts/eval_sizestrat_elements.csv --set confirm18=scripts/eval_confirm18_elements.csv \
  --set litpilot=scripts/literature_truth.tsv --set insitu="$SP/insitu_truth.tsv" \
  --set lockbox=analysis/lockbox/lockbox_elements.csv --set_genomes_dir lockbox=../raw/lockbox \
  --out analysis/integrity/containment.tsv || fail_step 6 $?

echo "=== job end (exit $rc) ==="; date
exit $rc
