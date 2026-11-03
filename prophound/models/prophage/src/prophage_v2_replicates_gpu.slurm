#!/bin/bash
# Job-array launcher for the 2026-09-17 replicate protocol: N replicates of one
# base config, each with its own run_name/seed/runs dir (docs/environment.md,
# docs/model-reengineering.md). One array task = one replicate.
#
# Submit from /work/$USER/mge_explorer/data_refactored, once per base config:
#   sbatch -J prophage_v2_pool      --export=ALL,BASE=prophage_v2_pool      src/prophage_v2_replicates_gpu.slurm
#   sbatch -J prophage_v2_pool_pw20 --export=ALL,BASE=prophage_v2_pool_pw20 src/prophage_v2_replicates_gpu.slurm
# Array task i trains runs/${BASE}_r${i}/config_snapshot.yaml (must already be
# pushed). Logs: runs/_slurm_logs/<jobname>_<arrayjob>_<task>.{out,err}.
# Re-running the same command after a TIMEOUT resumes each replicate from its
# own last.pt (auto_resume) -- that's the intended resume path.
#SBATCH -p gpu1v100
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --gres=gpu:1
#SBATCH -t 3-00:00:00
#SBATCH --array=1-3
#SBATCH -o runs/_slurm_logs/%x_%A_%a.out
#SBATCH -e runs/_slurm_logs/%x_%A_%a.err

set -euo pipefail

: "${BASE:?set BASE=<base run name>, e.g. --export=ALL,BASE=prophage_v2_pool}"
RUN_NAME="${BASE}_r${SLURM_ARRAY_TASK_ID}"
CONFIG_PATH="runs/${RUN_NAME}/config_snapshot.yaml"

echo "=== job start ==="
date
hostname
echo "array job ${SLURM_ARRAY_JOB_ID} task ${SLURM_ARRAY_TASK_ID} -> ${RUN_NAME}"

module load miniconda
conda activate /home/$USER/.conda/envs/mge_explorer_gpu

cd /work/$USER/mge_explorer/data_refactored
mkdir -p "runs/${RUN_NAME}/logs" "runs/${RUN_NAME}/checkpoints" "runs/${RUN_NAME}/metrics"

# train_prophage.py's walltime guard reads SLURM_TIMELIMIT, which Slurm does not
# export -- without this the guard never fires and the job is just SIGTERM'd at
# the limit (812125/812131 both were). See docs/environment.md.
export SLURM_TIMELIMIT="$(squeue -h -j "${SLURM_JOB_ID}" -o %l)"
echo "SLURM_TIMELIMIT=${SLURM_TIMELIMIT}"

echo "=== env ==="
which python; python -V
ulimit -n 65535 || true
nvidia-smi

echo "=== files ==="
ls -lh data/manifest_windows.filtered.csv splits/train.txt splits/val.txt splits/test.txt "$CONFIG_PATH" runs/prophage_v0/checkpoints/best.pt

echo "=== config preflight ==="
python - "$CONFIG_PATH" "$RUN_NAME" "$SLURM_ARRAY_TASK_ID" << 'PY'
import sys, yaml
from pathlib import Path
cfg_path, run_name, task = sys.argv[1], sys.argv[2], int(sys.argv[3])
cfg = yaml.safe_load(open(cfg_path))
t, l, m = cfg["training"], cfg.get("loss", {}), cfg["model"]
print("run_name:", cfg["run_name"], "| seed:", t["seed"], "| epochs:", t["epochs"], "| amp:", t["amp"],
      "| pos_weight:", l.get("pos_weight"), "| pooling:", m.get("pooling"), m.get("pool_bins"),
      "| init:", t.get("init_checkpoint"), t.get("init_mode"))
if cfg["run_name"] != run_name:
    raise SystemExit(f"ERROR: config run_name {cfg['run_name']!r} != expected {run_name!r} -- would train into the wrong runs/ dir")
if t["seed"] != 12 + task:
    raise SystemExit(f"ERROR: expected seed {12 + task} for replicate {task}, got {t['seed']}")
last = Path("runs") / run_name / "checkpoints" / "last.pt"
if last.exists():
    print(f"NOTE: {last} exists -- auto_resume will continue this replicate from it (expected only after a TIMEOUT resubmit)")
PY

echo "=== training ==="
python -m src.train_prophage --config "$CONFIG_PATH"

echo "=== job end ==="
date
