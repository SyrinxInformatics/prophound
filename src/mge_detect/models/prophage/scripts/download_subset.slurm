#!/bin/bash
#SBATCH -J dl_subset
#SBATCH -t 08:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p compute1
#SBATCH -o logs/dl_subset_%j.out

cd /work/pke835/mge_explorer
module load parallel
bash 01_scripts/download_subset.sh
