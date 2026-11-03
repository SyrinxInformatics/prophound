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

#SBATCH -J prophage_v1_positional_boundary
#SBATCH -p gpu2v100
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --gres=gpu:1
#SBATCH -t 2-00:00:00
#SBATCH -o runs/prophage_v1_positional_boundary/logs/slurm_%j.out
#SBATCH -e runs/prophage_v1_positional_boundary/logs/slurm_%j.err

set -euo pipefail

echo "=== job start ==="
date
hostname
pwd

module load miniconda
conda activate /home/$USER/.conda/envs/mge_explorer_gpu

cd /work/$USER/mge_explorer/data_refactored
CONFIG_PATH="${CONFIG_PATH:-runs/prophage_v1_positional_boundary/config_snapshot.yaml}"

echo "=== after cd ==="
pwd
which python
python -V

echo "=== files ==="
ls -lh "$CONFIG_PATH"
ls -lh runs/prophage_v0/checkpoints/best.pt

echo "=== boundary manifest check ==="
python - << PY
import yaml
p = "$CONFIG_PATH"
cfg = yaml.safe_load(open(p))
manifest = cfg["data"]["manifest_csv"]
print("manifest_csv:", manifest)
import csv
with open(manifest) as f:
    header = next(csv.reader(f))
if "boundary_dist_bp" not in header:
    raise SystemExit(
        f"ERROR: {manifest} has no boundary_dist_bp column -- run "
        "01_scripts/add_boundary_distance.py first (see "
        "01_scripts/add_boundary_distance.slurm)."
    )
print("OK: boundary_dist_bp column present")
PY

echo "=== config preflight ==="
python - << PY
import yaml
cfg = yaml.safe_load(open("$CONFIG_PATH"))
t, m, l = cfg.get("training", {}), cfg.get("model", {}), cfg.get("loss", {})
print("run_name:", cfg.get("run_name"))
print("epochs:", t.get("epochs"), "amp:", t.get("amp"), "init_mode:", t.get("init_mode"))
print("pooling:", m.get("pooling"), "pool_bins:", m.get("pool_bins"))
print("predict_boundary:", m.get("predict_boundary"))
print("pos_weight:", l.get("pos_weight"))
print("boundary_weight:", l.get("boundary_weight"), "boundary_max_dist_bp:", l.get("boundary_max_dist_bp"))
if not m.get("predict_boundary"):
    raise SystemExit("ERROR: this launcher is for the boundary-head experiment but model.predict_boundary is not true in the config")
PY

echo "=== training ==="
python -m src.train_prophage --config "$CONFIG_PATH"

echo "=== job end ==="
date
