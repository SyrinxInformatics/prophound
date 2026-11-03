#!/bin/bash
# docs/release-scope-0.1.0.md A.1 source check: where do inoviruses already sit in the
# corpus, and do current labels call any of them negative (label conflicts)?
# Needs Roux's pI HMMs at analysis/inovirus_overlap/Final_marker_morph.hmm
# (from bitbucket.org/srouxjgi/inovirus Inovirus_dbs.tar.gz); without
# hmmsearch on PATH it runs geNomad-taxonomy-only and says so.
# Submit from /work/$USER/mge_explorer/data_refactored:
#   sbatch scripts/inovirus_overlap_audit.slurm
#SBATCH -J inovirus_overlap
#SBATCH -p compute1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -t 08:00:00
#SBATCH -o analysis/_slurm_logs/%x_%j.out
#SBATCH -e analysis/_slurm_logs/%x_%j.err

set -uo pipefail
echo "=== job start ==="; date; hostname

module load miniconda
conda activate /home/$USER/.conda/envs/mge_explorer_gpu
cd /work/$USER/mge_explorer/data_refactored

GENOMAD=../genomad_out   # sibling of data_refactored/, not inside it
BEDS=beds_all
SPLITS=splits
OUT=analysis/inovirus_overlap
HMM=$OUT/Final_marker_morph.hmm
mkdir -p analysis/_slurm_logs "$OUT"

# preflight: 816440/816441 "COMPLETED" in 3 seconds having done nothing
fail=0
[ -d "$GENOMAD" ] || { echo "PREFLIGHT FAIL: missing $GENOMAD" >&2; fail=1; }
[ -d "$BEDS" ] || { echo "PREFLIGHT FAIL: missing $BEDS" >&2; fail=1; }
[ -s "$SPLITS/train.txt" ] || { echo "PREFLIGHT FAIL: missing $SPLITS/train.txt" >&2; fail=1; }
EVAL_ARGS=()
for L in scripts/eval_genome_list.txt scripts/eval_genome_list_sizestrat.txt; do
  [ -s "$L" ] || { echo "PREFLIGHT FAIL: missing $L" >&2; fail=1; }
  EVAL_ARGS+=(--eval_list "$L")
done
n=$(ls -d "$GENOMAD"/*/ 2>/dev/null | wc -l); [ "$n" -gt 0 ] || { echo "PREFLIGHT FAIL: $GENOMAD is empty" >&2; fail=1; }
na=$(ls "$GENOMAD"/*/*_annotate/*_genes.tsv 2>/dev/null | wc -l)
[ "$na" -gt 0 ] || { echo "PREFLIGHT FAIL: no *_annotate/*_genes.tsv under $GENOMAD" >&2; fail=1; }
[ $fail -eq 0 ] || { echo "PREFLIGHT FAILED -- nothing run" >&2; exit 1; }
echo "[preflight] $n genomad dirs, $na with annotate genes.tsv"

PI_ARGS=()
if [ -s "$HMM" ] && command -v hmmsearch >/dev/null; then
  PI_ARGS=(--pi_hmm "$HMM"); echo "[preflight] pI HMM + $(command -v hmmsearch)"
else
  echo "WARNING: pI evidence DISABLED (hmm present: $([ -s "$HMM" ] && echo y || echo n), hmmsearch: $(command -v hmmsearch || echo none)) -- tax-only run" >&2
fi

python scripts/inovirus_overlap_audit.py \
  --genomad_dir "$GENOMAD" --beds_dir "$BEDS" --splits_dir "$SPLITS" \
  "${EVAL_ARGS[@]}" "${PI_ARGS[@]}" --cpus "$SLURM_CPUS_PER_TASK" --out_dir "$OUT"
rc=$?

echo "=== job end (python exit $rc) ==="; date
# 817355 crashed mid-scan yet SLURM reported COMPLETED 0:0: without -e the job's
# status is the last command's, i.e. this echo. Propagate python's instead.
exit $rc
