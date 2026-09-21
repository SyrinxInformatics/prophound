#!/bin/bash
# docs/scientific-validity.md item 2a, stage 1: cluster geNomad-called elements
# into clades so a clade-held-out split can be built. CPU only.
# NOTHING IS RETRAINED HERE -- this job exists so the clade statistics can be
# inspected before committing GPU days.
# Submit from /work/$USER/mge_explorer/data_refactored:
#   sbatch scripts/cluster_prophage_elements.slurm
#SBATCH -J cluster_elements
# Walltime kept short on purpose: compute1 showed 0 idle nodes at submission
# time, and a short request backfills far sooner than a padded one.
# ~4k elements sketched + ~8M sketch comparisons: under an hour expected, 4h is generous.
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
mkdir -p analysis/_slurm_logs analysis/clades

echo "=== inputs ==="
echo "beds_pos: $(ls beds_pos | wc -l) files"; ls -d sequences

# --fasta_suffix must match what is actually in sequences/ -- the preflight
# prints a sample; change it here if they are .fasta rather than .fna.
python scripts/cluster_prophage_elements.py \
  --beds_dir beds_pos \
  --sequences_dir sequences \
  --out_dir analysis/clades \
  --k 21 --sketch 2000 --min_jaccard 0.2

echo "=== job end ==="; date
