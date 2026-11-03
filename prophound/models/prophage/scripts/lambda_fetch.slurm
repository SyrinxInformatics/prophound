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
# LAMBDA clean subset (docs/claims-and-evidence.md "Amendment 2026-09-29"), step 1 of 2. CPU.
#   1. download the one sequence per genome that LAMBDA scores (the chromosome
#      carrying the curated elements; lambda_genomes.tsv) from NCBI nuccore,
#      header rewritten to the bare accession so it matches the truth table,
#      and check its length against LAMBDA's genome_length (rounded up to 1 kb)
#   2. 21-mer containment of every curated element against training positives
#      -> analysis/integrity/containment_lambda.tsv
# No model is run. Submit from /work/$USER/mge_explorer/data_refactored:
#   sbatch scripts/lambda_fetch.slurm
#SBATCH -J lambda_fetch
#SBATCH -p compute1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -t 06:00:00
#SBATCH -o analysis/_slurm_logs/%x_%j.out
#SBATCH -e analysis/_slurm_logs/%x_%j.err

set -uo pipefail
echo "=== job start ==="; date; hostname

module load miniconda
conda activate /home/$USER/.conda/envs/mge_explorer_gpu
cd /work/$USER/mge_explorer/data_refactored

L=analysis/external_benchmarks/lambda
FNA=../raw/lambda                # sibling of raw/genomes, kept apart from the corpus
mkdir -p analysis/_slurm_logs analysis/integrity "$FNA"

# ---- preflight -------------------------------------------------------------
fail=0
for f in "$L/lambda_genomes.tsv" "$L/lambda_truth_manifest.csv" scripts/eval_containment.py \
         scripts/small_element_spikein.py splits/train.txt; do
  [ -s "$f" ] || { echo "PREFLIGHT FAIL: missing or empty $f" >&2; fail=1; }
done
for d in beds_pos sequences; do [ -d "$d" ] || { echo "PREFLIGHT FAIL: no directory $d" >&2; fail=1; }; done
# a real GET: E-utilities answers HEAD (curl -I) with 405, which failed 822256's preflight
curl -fsS -o /dev/null "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi?db=nuccore" \
  || { echo "PREFLIGHT FAIL: no network access to NCBI from $(hostname)" >&2; fail=1; }
[ $fail -eq 0 ] || { echo "PREFLIGHT FAILED -- nothing run" >&2; exit 1; }
echo "[preflight] all inputs present"
# ---------------------------------------------------------------------------

rc=0
tail -n +2 "$L/lambda_genomes.tsv" | while IFS=$'\t' read -r ASM SEQ ORG GLEN STATUS GENUS; do
  out="$FNA/$ASM.fna"
  if [ ! -s "$out" ]; then
    curl -fsSL --retry 3 "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$SEQ&rettype=fasta&retmode=text" \
      | awk -v id="$SEQ" 'NR==1{print ">" id; next} {print}' > "$out.tmp" && mv "$out.tmp" "$out" \
      || { echo "FAILED download $ASM $SEQ" >> "$L/fetch_failures.txt"; rm -f "$out.tmp"; continue; }
    sleep 0.4                                   # NCBI allows 3 requests/s without a key
  fi
  len=$(grep -v '^>' "$out" | tr -d '\n' | wc -c)
  # LAMBDA's genome_length is the length rounded up to the next kb, plus 1 kb
  # (checked on all 39); 822260 used a <1 kb tolerance and flagged every genome
  want=$(( (len + 999) / 1000 * 1000 + 1000 ))
  if [ "$want" -ne "$GLEN" ]; then
    echo "FAILED length $ASM $SEQ: $len bp vs genome_length $GLEN" >> "$L/fetch_failures.txt"
  fi
done
echo "[fetch] $(ls "$FNA"/*.fna 2>/dev/null | wc -l) of $(($(wc -l < "$L/lambda_genomes.tsv") - 1)) sequences"; date
[ -s "$L/fetch_failures.txt" ] && { echo "WARNING: $(wc -l < "$L/fetch_failures.txt") fetch problems" >&2; cat "$L/fetch_failures.txt" >&2; rc=2; }

python scripts/eval_containment.py --beds_dir beds_pos --sequences_dir sequences --train_split splits/train.txt \
  --genomes_dir "$FNA" --set lambda="$L/lambda_truth_manifest.csv" \
  --out analysis/integrity/containment_lambda.tsv || { r=$?; echo "FAILED containment (exit $r)" >&2; rc=$r; }

echo "=== job end (exit $rc) ==="; date
exit $rc
