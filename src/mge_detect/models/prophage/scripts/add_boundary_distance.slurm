#!/bin/bash
#SBATCH -J add_boundary_distance
# NOTE: partition name below is a guess (this repo has no CPU-only SLURM
# launcher to copy from, only GPU ones like prophage_v1_gpu.slurm's
# gpu1v100) -- check `sinfo` on ARC and correct this before submitting.
#SBATCH -p compute1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=32G
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
