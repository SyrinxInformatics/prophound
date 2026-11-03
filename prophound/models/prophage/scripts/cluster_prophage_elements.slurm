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
# docs/scientific-validity.md item 2a, stage 1: cluster geNomad-called elements
# into clades so a clade-held-out split can be built. CPU only.
# NOTHING IS RETRAINED HERE -- this job exists so the clade statistics can be
# inspected before committing GPU days.
# Submit from /work/$USER/mge_explorer/data_refactored:
#   sbatch scripts/cluster_prophage_elements.slurm
#SBATCH -J cluster_elements
# Walltime kept short on purpose: compute1 showed 0 idle nodes at submission
# time, and a short request backfills far sooner than a padded one.
# 815766 (2026-09-21) took 1h04 for ~4k elements / 7.9M pairs; --thresholds scores
# every pair (no union short-circuit) so budget the same, 4h stays generous.
#SBATCH -p compute1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -t 04:00:00
#SBATCH -o analysis/_slurm_logs/%x_%j.out
#SBATCH -e analysis/_slurm_logs/%x_%j.err

set -euo pipefail
echo "=== job start ==="; date; hostname

module load miniconda
conda activate /home/$USER/.conda/envs/mge_explorer_gpu
cd /work/$USER/mge_explorer/data_refactored
# analysis/clades/ holds 815766's single-threshold output; the sweep writes
# beside it so that file stays as pulled/analysed.
mkdir -p analysis/_slurm_logs analysis/clades_sweep

echo "=== inputs ==="
echo "beds_pos: $(ls beds_pos | wc -l) files"; ls -d sequences

# --fasta_suffix must match what is actually in sequences/ -- the preflight
# prints a sample; change it here if they are .fasta rather than .fna.
python scripts/cluster_prophage_elements.py \
  --beds_dir beds_pos \
  --sequences_dir sequences \
  --out_dir analysis/clades_sweep \
  --k 21 --sketch 2000 --min_jaccard 0.2 \
  --thresholds 0.005,0.01,0.02,0.05,0.1,0.15
# k=21 makes Jaccard a steep function of divergence (J ~ p^21/(2-p^21)):
# 0.2 ~ 3.5% divergent, 0.1 ~ 7%, 0.05 ~ 10%, 0.02 ~ 14%, 0.01 ~ 17%, 0.005 ~ 20%.
# pairs_jaccard.csv keeps every pair >= 0.005, so lower thresholds than that
# need a re-run; anything above is free locally.

echo "=== job end ==="; date
