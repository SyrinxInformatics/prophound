#!/bin/bash
#SBATCH -J prophage_v0
#SBATCH -p gpu1v100
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --gres=gpu:1
#SBATCH -t 2-00:00:00
#SBATCH -o runs/prophage_v0/logs/slurm_%j.out
#SBATCH -e runs/prophage_v0/logs/slurm_%j.err

set -euo pipefail

echo "=== job start ==="
date
hostname
pwd

module load miniconda
conda activate /home/$USER/.conda/envs/mge_explorer_gpu

cd /work/$USER/mge_explorer/data_refactored
CONFIG_PATH="runs/prophage_v0/config_snapshot.yaml"

echo "=== after cd ==="
pwd
which python
python -V

echo "=== GPU ==="
nvidia-smi

echo "=== files ==="
ls -lh data/manifest_windows.csv
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
print("epochs:", t.get("epochs"))
print("auto_resume:", t.get("auto_resume"))
print("amp:", t.get("amp"))
print("manifest_csv:", d.get("manifest_csv"))
PY

echo "=== training ==="
python -m src.train_prophage --config "$CONFIG_PATH"

echo "=== job end ==="
date
