#!/bin/bash
# SPDX-License-Identifier: Apache-2.0
#
# Copyright 2025 Jacob Alford <jalford0000@gmail.com>
#                and Benjamin Alford <ssolar5.dev@gmail.com>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# docs/research-integrity.md C5 + C6 (CPU; no model scoring):
#   C5  21-mer containment of every evaluation element (42, 18, literature pilot,
#       in-situ inovirus loci, lockbox) against training-split positives
#       -> analysis/integrity/containment.tsv
#   C6  off-target annotation re-run with model-independent controls, avoiding
#       geNomad elements only, and fully random; same calls42.csv as 821371
#       -> analysis/offtarget/annotation_ctrl_{elements_only,none}.{json,csv}
# Submit from /work/$USER/mge_explorer/data_refactored:
#   sbatch scripts/integrity_cpu.slurm
#SBATCH -J integrity_cpu
#SBATCH -p compute1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -t 12:00:00
#SBATCH -o analysis/_slurm_logs/%x_%j.out
#SBATCH -e analysis/_slurm_logs/%x_%j.err

set -uo pipefail
echo "=== job start ==="; date; hostname

module load miniconda
conda activate /home/$USER/.conda/envs/mge_explorer_gpu
cd /work/$USER/mge_explorer/data_refactored

GENOMAD=../genomad_out      # siblings of data_refactored/, not inside it
GENOMES=../raw/genomes
CALLS=analysis/offtarget/calls42.csv
INSITU=analysis/spikein_small/insitu_truth.tsv
mkdir -p analysis/_slurm_logs analysis/integrity analysis/offtarget

# ---- preflight (docs/environment.md) ---------------------------------------
fail=0
for f in scripts/eval_containment.py scripts/small_element_spikein.py scripts/offtarget_annotation.py \
         scripts/eval_sizestrat_elements.csv scripts/eval_confirm18_elements.csv scripts/literature_truth.tsv \
         "$INSITU" "$CALLS" splits/train.txt analysis/lockbox/lockbox_elements.csv; do
  [ -s "$f" ] || { echo "PREFLIGHT FAIL: missing or empty $f" >&2; fail=1; }
done
grep -q -- '--control_busy' scripts/offtarget_annotation.py \
  || { echo "PREFLIGHT FAIL: scripts/offtarget_annotation.py lacks --control_busy (not pushed?)" >&2; fail=1; }
for d in "$GENOMAD" "$GENOMES" ../raw/lockbox beds_pos sequences; do [ -d "$d" ] || { echo "PREFLIGHT FAIL: no directory $d" >&2; fail=1; }; done
[ $fail -eq 0 ] || { echo "PREFLIGHT FAILED -- nothing run" >&2; exit 1; }
echo "[preflight] all inputs present"
# ---------------------------------------------------------------------------

rc=0
echo "[C5] containment"; date
python scripts/eval_containment.py --beds_dir beds_pos --sequences_dir sequences --train_split splits/train.txt \
  --genomes_dir "$GENOMES" \
  --set sizestrat42=scripts/eval_sizestrat_elements.csv --set confirm18=scripts/eval_confirm18_elements.csv \
  --set litpilot=scripts/literature_truth.tsv --set insitu=$INSITU \
  --set lockbox=analysis/lockbox/lockbox_elements.csv --set_genomes_dir lockbox=../raw/lockbox \
  --out analysis/integrity/containment.tsv || { r=$?; echo "FAILED C5 (exit $r)" >&2; rc=$r; }

for B in elements_only none; do
  echo "[C6] controls avoiding: $B"; date
  python scripts/offtarget_annotation.py --calls "$CALLS" --elements scripts/eval_sizestrat_elements.csv \
    --genomad_dir "$GENOMAD" --genomes_dir "$GENOMES" --control_busy $B \
    --out analysis/offtarget/annotation_ctrl_$B.json || { r=$?; echo "FAILED C6 $B (exit $r)" >&2; rc=$r; }
done

echo "=== job end (exit $rc) ==="; date
exit $rc
