#!/bin/bash
#SBATCH -J prophage_v1_r_val_test
#SBATCH -p gpu1v100
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --gres=gpu:1
#SBATCH -t 2-00:00:00
#SBATCH -o runs/prophage_v1_r_val_test/logs/slurm_%j.out
#SBATCH -e runs/prophage_v1_r_val_test/logs/slurm_%j.err

set -euo pipefail

echo "=== job start ==="
date
hostname
pwd

module load miniconda
conda activate /home/$USER/.conda/envs/mge_explorer_gpu

cd /work/$USER/mge_explorer/data_refactored

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

echo "=== torch ==="
python - << 'PY'
import torch
print("torch:", torch.__version__)
print("cuda available:", torch.cuda.is_available())
if torch.cuda.is_available():
    print("gpu:", torch.cuda.get_device_name(0))
PY

echo "=== split test ==="
python -m src.train_prophage \
  --config runs/prophage_v1/config_snapshot.yaml \
  --eval-only \
  --split test \
  --checkpoint runs/prophage_v1/checkpoints/best.pt

python -m src.train_prophage \
  --config runs/prophage_v1_r2/config_snapshot.yaml \
  --eval-only \
  --split test \
  --checkpoint runs/prophage_v1_r2/checkpoints/best.pt

python -m src.train_prophage \
  --config runs/prophage_v1_r3/config_snapshot.yaml \
  --eval-only \
  --split test \
  --checkpoint runs/prophage_v1_r3/checkpoints/best.pt

echo "=== job end ==="
date
