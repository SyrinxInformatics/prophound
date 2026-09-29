#!/bin/bash
#SBATCH -J gmd_afp
#SBATCH -p compute1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -t 12:00:00
#SBATCH -o logs/gmd_afp_%A_%a.out
set -euo pipefail

module load miniconda
conda activate /home/pke835/.conda/envs/mge
cd /work/pke835/mge_explorer

: "${CHUNK_FILE:?CHUNK_FILE is not set}"
GENOME="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$CHUNK_FILE")" || exit 0
[ -n "$GENOME" ] || exit 0

BASE="$(basename "$GENOME" .fna)"
OUTDIR="genomad_out/$BASE"
ANNOT_DIR="$OUTDIR/${BASE}_annotate"
PROV_DIR="$OUTDIR/${BASE}_find_proviruses"

mkdir -p "$OUTDIR" logs

# per-task scratch for mmseqs2/Prokka/Prodigal
SCR_BASE="/scratch/$USER"
[ -d "$SCR_BASE" ] || SCR_BASE="/work/$USER/tmp"
TMPDIR="$SCR_BASE/gmd_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$TMPDIR"
export TMPDIR MMSEQS_TMP="$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT

# Skip if a provirus table already exists
if [ -f "$PROV_DIR/${BASE}_provirus.tsv" ]; then
  echo "[afp] provirus TSV exists → $BASE (skip)"
  exit 0
fi

# Ensure annotate outputs exist (genes/proteins live under *_annotate/)
if [ ! -f "$ANNOT_DIR/${BASE}_genes.tsv" ] || [ ! -f "$ANNOT_DIR/${BASE}_proteins.faa" ]; then
  echo "[afp] running annotate for $BASE"
  genomad annotate "$GENOME" "$OUTDIR" genomad_db
fi

# Now call proviruses (outputs go under *_find_proviruses/)
echo "[afp] running find-proviruses for $BASE"
genomad find-proviruses "$GENOME" "$OUTDIR" genomad_db

# Verify
if [ ! -f "$PROV_DIR/${BASE}_provirus.tsv" ]; then
  echo "[afp] find-proviruses finished but no ${BASE}_provirus.tsv" >&2
  exit 1
fi
