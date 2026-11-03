#!/bin/bash
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
