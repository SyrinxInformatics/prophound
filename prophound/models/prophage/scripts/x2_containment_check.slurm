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
# Shared-code check: re-run eval_containment.py with the evalkit-based scripts
# and compare its output byte for byte with the recorded analysis/integrity/containment.tsv.
# The current src/ and the two scripts are pushed to x2_check/ (beside the frozen tree, which
# stays untouched); the data paths are the original job's, so it runs from data_refactored/:
#   sbatch scripts/x2_containment_check.slurm      (after extracting x2_check.tar there)
#SBATCH -J x2_containment
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

X=x2_check
SP=analysis/spikein_small
fail=0
for f in $X/src/calls.py $X/src/predict_genome.py $X/scripts/evalkit.py $X/scripts/eval_containment.py \
         beds_pos sequences splits/train.txt ../raw/genomes ../raw/lockbox analysis/integrity/containment.tsv \
         scripts/eval_sizestrat_elements.csv scripts/eval_confirm18_elements.csv scripts/literature_truth.tsv \
         "$SP/insitu_truth.tsv" analysis/lockbox/lockbox_elements.csv; do
  [ -e "$f" ] || { echo "PREFLIGHT FAIL: missing $f" >&2; fail=1; }
done
[ $fail -eq 0 ] || exit 1

python $X/scripts/eval_containment.py --beds_dir beds_pos --sequences_dir sequences --train_split splits/train.txt \
  --genomes_dir ../raw/genomes \
  --set sizestrat42=scripts/eval_sizestrat_elements.csv --set confirm18=scripts/eval_confirm18_elements.csv \
  --set litpilot=scripts/literature_truth.tsv --set insitu="$SP/insitu_truth.tsv" \
  --set lockbox=analysis/lockbox/lockbox_elements.csv --set_genomes_dir lockbox=../raw/lockbox \
  --out $X/containment.tsv
rc=$?
if [ $rc -eq 0 ]; then
  if cmp -s $X/containment.tsv analysis/integrity/containment.tsv; then
    echo "X2 CONTAINMENT: IDENTICAL to analysis/integrity/containment.tsv"
  else
    echo "X2 CONTAINMENT: DIFFERS from analysis/integrity/containment.tsv" >&2
    rc=3
  fi
fi
echo "=== job end (rc $rc) ==="; date
exit $rc
