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

#SBATCH -J prophage_v1_r2
#SBATCH -p gpu1v100
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --gres=gpu:1
#SBATCH -t 2-20:00:00
#SBATCH -o runs/prophage_v1_r2/logs/slurm_%j.out
#SBATCH -e runs/prophage_v1_r2/logs/slurm_%j.err

set -euo pipefail

echo "=== job start ==="
date
hostname
pwd

module load miniconda
conda activate /home/$USER/.conda/envs/mge_explorer_gpu

cd /work/$USER/mge_explorer/data_refactored
CONFIG_PATH="${CONFIG_PATH:-runs/prophage_v1_r2/config_snapshot.yaml}"

echo "=== after cd ==="
pwd
which python
python -V
echo "=== ulimit (nofile) ==="
ulimit -n || true
ulimit -n 65535 || true
echo "nofile now: $(ulimit -n || true)"

echo "=== GPU ==="
nvidia-smi

echo "=== files ==="
ls -lh data/manifest_windows.filtered.csv
ls -lh splits/train.txt splits/val.txt splits/test.txt
ls -lh "$CONFIG_PATH"

echo "=== torch ==="
python - << 'PY'
import torch
print("torch:", torch.__version__)
print("cuda available:", torch.cuda.is_available())
if torch.cuda.is_available():
    print("gpu:", torch.cuda.get_device_name(0))
PY

echo "=== config preflight ==="
echo "CONFIG_PATH=$CONFIG_PATH"
python - << PY
import yaml
p = "$CONFIG_PATH"
cfg = yaml.safe_load(open(p))
t = cfg.get("training", {})
d = cfg.get("data", {})
print("run_name:", cfg.get("run_name"))
print("epochs:", t.get("epochs"))
print("auto_resume:", t.get("auto_resume"))
print("amp:", t.get("amp"))
print("lr:", t.get("lr"))
print("train neg types:", d.get("use_neg_types_in_train"))
print("manifest_csv:", d.get("manifest_csv"))
PY

echo "=== training ==="
python -m src.train_prophage --config "$CONFIG_PATH"

echo "=== job end ==="
date
