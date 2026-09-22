#!/bin/bash
# docs/v1-scope.md gating audit A.1: are small elements represented as positives
# in training, and at what dilution? Gates the retraining data design.
# Submit from /work/$USER/mge_explorer/data_refactored:
#   sbatch scripts/size_representation_audit.slurm
#SBATCH -J size_representation
#SBATCH -p compute1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -t 04:00:00
#SBATCH -o analysis/_slurm_logs/%x_%j.out
#SBATCH -e analysis/_slurm_logs/%x_%j.err

set -uo pipefail
echo "=== job start ==="; date; hostname

module load miniconda
conda activate /home/$USER/.conda/envs/mge_explorer_gpu
cd /work/$USER/mge_explorer/data_refactored

MANIFEST=data/manifest_windows.filtered.csv
BEDS=beds_pos
SPLITS=splits
mkdir -p analysis/_slurm_logs analysis/size_representation

# preflight: 816440/816441 "COMPLETED" in 3 seconds having done nothing
fail=0
[ -s "$MANIFEST" ] || { echo "PREFLIGHT FAIL: missing $MANIFEST" >&2; fail=1; }
[ -d "$BEDS" ] || { echo "PREFLIGHT FAIL: missing $BEDS" >&2; fail=1; }
[ -s "$SPLITS/train.txt" ] || { echo "PREFLIGHT FAIL: missing $SPLITS/train.txt" >&2; fail=1; }
n=$(ls "$BEDS" 2>/dev/null | wc -l); [ "$n" -gt 0 ] || { echo "PREFLIGHT FAIL: $BEDS is empty" >&2; fail=1; }
[ $fail -eq 0 ] || { echo "PREFLIGHT FAILED -- nothing run" >&2; exit 1; }
echo "[preflight] manifest, $n BED files, splits present"

python scripts/size_representation_audit.py \
  --manifest "$MANIFEST" \
  --beds_dir "$BEDS" \
  --splits_dir "$SPLITS" \
  --out_dir analysis/size_representation \
  --split train

echo "=== job end ==="; date
