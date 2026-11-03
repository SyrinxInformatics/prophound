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
# docs/cli-and-packaging.md "Input robustness": how much sequence, and how much
# known truth, sits in contigs too short to be scored? CPU only, no inference.
# Submit from /work/$USER/mge_explorer/data_refactored:
#   sbatch scripts/short_contig_audit.slurm
#SBATCH -J short_contigs
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
mkdir -p analysis/_slurm_logs analysis/short_contigs

echo "=== inputs ==="
ls -d sequences beds_pos

python scripts/short_contig_audit.py \
  --sequences_dir sequences \
  --beds_dir beds_pos \
  --out_dir analysis/short_contigs \
  --window_len 10000

echo "=== job end ==="; date
