#!/bin/bash
# docs/scientific-validity.md item 1a: did the model learn "ambiguity code =>
# not a prophage" from the training labels? CPU only, I/O bound.
# Submit from /work/$USER/mge_explorer/data_refactored:
#   sbatch scripts/ambiguity_label_audit.slurm
#SBATCH -J ambiguity_labels
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
mkdir -p analysis/_slurm_logs analysis/ambiguity_labels

echo "=== inputs ==="
ls -lh data/manifest_windows.filtered.csv splits/train.txt

python scripts/ambiguity_label_audit.py \
  --manifest data/manifest_windows.filtered.csv \
  --splits_dir splits \
  --out_dir analysis/ambiguity_labels \
  --n_sample 400000

echo "=== job end ==="; date
