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
# Composition-only baseline (analysis/composition_baseline/): can GC/k-mer composition alone match the
# CNN's window-level AUPRC of 0.63-0.65? CPU only.
# Submit from /work/$USER/mge_explorer/data_refactored:
#   sbatch scripts/composition_baseline.slurm
#SBATCH -J composition_baseline
# Walltime kept short on purpose: compute1 showed 0 idle nodes at submission
# time, and a short request backfills far sooner than a padded one.
# ~150k+50k windows featurized then two sklearn fits: measured locally at well under an hour of compute, 4h is generous.
#SBATCH -p compute1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -t 04:00:00
#SBATCH -o analysis/_slurm_logs/%x_%j.out
#SBATCH -e analysis/_slurm_logs/%x_%j.err

set -euo pipefail
echo "=== job start ==="; date; hostname

module load miniconda
conda activate /home/$USER/.conda/envs/mge_explorer_gpu
cd /work/$USER/mge_explorer/data_refactored
mkdir -p analysis/_slurm_logs analysis/composition_baseline

echo "=== inputs ==="
ls -lh data/manifest_windows.filtered.csv splits/train.txt splits/test.txt

python scripts/composition_baseline.py \
  --manifest data/manifest_windows.filtered.csv \
  --splits_dir splits \
  --out_dir analysis/composition_baseline \
  --n_train 150000 --n_test 50000 --k 4

echo "=== job end ==="; date
