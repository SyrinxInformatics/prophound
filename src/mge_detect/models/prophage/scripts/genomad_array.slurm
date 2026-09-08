#!/bin/bash
#SBATCH -J genomad
#SBATCH -p compute1
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -t 2-00:00:00
#SBATCH -o logs/gmd_%A_%a.out

set -euo pipefail

module load miniconda
conda activate /home/$USER/.conda/envs/mge_explorer
cd /work/$USER/mge_explorer/plasmids/genome_chunks

DB="/work/$USER/mge_explorer/genomad_db/genomad_db"
CHUNK_FILE="${CHUNK_FILE:?CHUNK_FILE not set or empty}"
[ -r "$CHUNK_FILE" ] || { echo "[error] cannot read $CHUNK_FILE"; exit 1; }

mkdir -p logs
SCR="/scratch/$USER/${SLURM_JOB_ID:-job}"
[ -d "/scratch/$USER" ] || SCR="/work/$USER/tmp/${SLURM_JOB_ID:-job}"
mkdir -p "$SCR"
export TMPDIR="$SCR" MMSEQS_TMP="$SCR"

THREADS=4                 # threads per genome
P=$(( SLURM_CPUS_PER_TASK / THREADS ))   # concurrent genomes
[ "$P" -ge 1 ] || P=1

xargs -a "$CHUNK_FILE" -n1 -P"$P" -I{} bash -lc '
  f="{}"; b=$(basename "$f" .fna); od="genomad_out/$b"
  ad="$od/${b}_annotate"; pd="$od/${b}_find_proviruses"
  mkdir -p "$od"

  if [ -s "$pd/${b}_provirus.tsv" ]; then
    echo "[skip] $b already done"; exit 0
  fi

  echo "[annotate] $b"
  genomad annotate --threads '"$THREADS"' "$f" "$od" "'"$DB"'" > "logs/gmd_${b}.log" 2>&1

  if [ ! -s "$ad/${b}_genes.tsv" ] || [ ! -s "$ad/${b}_proteins.faa" ] || [ ! -s "$ad/${b}_mmseqs2.tsv" ]; then
    echo "[error] $b annotate outputs missing"; exit 2
  fi

  echo "[find-proviruses] $b"
  genomad find-proviruses --threads '"$THREADS"' "$f" "$od" "'"$DB"'" >> "logs/gmd_${b}.log" 2>&1
  echo "[ok] $b"
'
