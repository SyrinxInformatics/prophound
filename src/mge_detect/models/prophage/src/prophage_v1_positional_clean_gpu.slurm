#!/bin/bash
#SBATCH -J prophage_v1_positional_clean
#SBATCH -p gpu1v100
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --gres=gpu:1
#SBATCH -t 2-00:00:00
#SBATCH -o runs/prophage_v1_positional_clean/logs/slurm_%j.out
#SBATCH -e runs/prophage_v1_positional_clean/logs/slurm_%j.err

set -euo pipefail

echo "=== job start ==="
date
hostname
pwd

module load miniconda
conda activate /home/$USER/.conda/envs/mge_explorer_gpu

cd /work/$USER/mge_explorer/data_refactored
CONFIG_PATH="${CONFIG_PATH:-runs/prophage_v1_positional_clean/config_snapshot.yaml}"

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
ls -lh runs/prophage_v0/checkpoints/best.pt

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
m = cfg.get("model", {})
l = cfg.get("loss", {})
print("run_name:", cfg.get("run_name"))
print("epochs:", t.get("epochs"))
print("amp:", t.get("amp"))
print("init_checkpoint:", t.get("init_checkpoint"))
print("init_mode:", t.get("init_mode"))
print("pos_weight:", l.get("pos_weight"))
print("pooling:", m.get("pooling"), "pool_bins:", m.get("pool_bins"))
# Sanity check this actually matches v1 except pooling/init_mode.
v1_cfg = yaml.safe_load(open("runs/prophage_v1/config_snapshot.yaml"))
mismatches = []
for section in ("training", "loss"):
    for key, v1_val in v1_cfg.get(section, {}).items():
        this_val = cfg.get(section, {}).get(key)
        if key == "init_mode":
            continue  # expected to differ, see config comment
        if this_val != v1_val:
            mismatches.append((section, key, v1_val, this_val))
if mismatches:
    print("WARNING: unexpected differences from v1's config:")
    for section, key, v1_val, this_val in mismatches:
        print(f"  {section}.{key}: v1={v1_val!r} this={this_val!r}")
else:
    print("OK: training/loss hyperparameters match v1 exactly (besides init_mode).")
PY

echo "=== training ==="
python -m src.train_prophage --config "$CONFIG_PATH"

echo "=== job end ==="
date
