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

#SBATCH -J add_boundary_distance
#SBATCH -p compute1
#SBATCH -n 1
#SBATCH -c 4
# No explicit --mem: this partition's DefMemPerNode/MaxMemPerNode are both
# UNLIMITED, and an earlier attempt with --mem=32G failed with "Memory
# specification can not be satisfied" despite nodes reporting 300+ GB free
# (sinfo's raw MEMORY field also looked bogus, e.g. "1") -- looks like a
# misconfigured/miscalibrated per-node memory value in Slurm's node table
# that any explicit --mem trips over. Omitting it relies on the unlimited
# default instead of that broken table entry.
#SBATCH -t 0-04:00:00
#SBATCH -o logs/add_boundary_distance_%j.out
#SBATCH -e logs/add_boundary_distance_%j.err

set -euo pipefail

echo "=== job start ==="
date
hostname
pwd

module load miniconda
conda activate /home/$USER/.conda/envs/mge_explorer_gpu

cd /work/$USER/mge_explorer/data_refactored
mkdir -p logs

python 01_scripts/add_boundary_distance.py \
  --manifest_csv data/manifest_windows.filtered.csv \
  --positives_dir /work/$USER/mge_explorer/positives/by_genome \
  --out_csv data/manifest_windows.filtered.with_boundary.csv

echo "=== job end ==="
date
