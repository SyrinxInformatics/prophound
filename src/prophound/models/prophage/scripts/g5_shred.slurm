#!/bin/bash
# G5 fragmented-input comparison (docs/claims-and-evidence.md "G5 preregistration"),
# step 1 of 3. CPU, minutes. Shreds the 39 LAMBDA chromosomes (../raw/lambda, from
# lambda_fetch.slurm) into the four preregistered conditions and checks every
# input and every shred byte-for-byte against analysis/g5/input_md5.tsv, which
# was written locally from the same fetch; the local Jaeger runs use those same
# files. Nothing is scored here. Submit all three steps from
# /work/$USER/mge_explorer/data_refactored:
#   S=$(sbatch --parsable scripts/g5_shred.slurm) && \
#   sbatch --dependency=afterok:$S scripts/g5_prophound.slurm && \
#   sbatch --dependency=afterok:$S scripts/g5_genomad.slurm
#SBATCH -J g5_shred
#SBATCH -p compute1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH -t 01:00:00
#SBATCH -o analysis/_slurm_logs/%x_%j.out
#SBATCH -e analysis/_slurm_logs/%x_%j.err

set -uo pipefail
echo "=== job start ==="; date; hostname

module load miniconda
conda activate /home/$USER/.conda/envs/mge_explorer_gpu
cd /work/$USER/mge_explorer/data_refactored

G=analysis/g5
FNA=../raw/lambda
SHR=../raw/g5_shreds             # sibling of raw/lambda, kept apart from the corpus
mkdir -p analysis/_slurm_logs "$SHR"

# ---- preflight (docs/environment.md) ---------------------------------------
fail=0
for f in "$G/genome_list.txt" "$G/input_md5.tsv" scripts/shred_genomes.py; do
  [ -s "$f" ] || { echo "PREFLIGHT FAIL: missing or empty $f" >&2; fail=1; }
done
grep -q -- '--dist' scripts/shred_genomes.py \
  || { echo "PREFLIGHT FAIL: scripts/shred_genomes.py lacks --dist lognormal (not pushed?)" >&2; fail=1; }
n=$(ls "$FNA"/*.fna 2>/dev/null | wc -l)
[ "$n" -eq 39 ] || { echo "PREFLIGHT FAIL: $FNA has $n/39 sequences" >&2; fail=1; }
# the complete chromosomes must be the ones the local Jaeger runs used
while IFS=$'\t' read -r path want; do
  case "$path" in fasta/*) got=$(md5sum "$FNA/${path#fasta/}" | cut -d' ' -f1)
    [ "$got" = "$want" ] || { echo "PREFLIGHT FAIL: $FNA/${path#fasta/} md5 $got != $want" >&2; fail=1; } ;; esac
done < "$G/input_md5.tsv"
[ $fail -eq 0 ] || { echo "PREFLIGHT FAILED -- nothing run" >&2; exit 1; }
echo "[preflight] 39 chromosomes match the local fetch"
# ---------------------------------------------------------------------------

rc=0
shred() {  # condition, then shred_genomes.py arguments
  local c=$1; shift
  rm -rf "$SHR/$c.tmp"
  python scripts/shred_genomes.py --genome_list "$G/genome_list.txt" --genome_dir "$FNA" --out_dir "$SHR/$c.tmp" \
    --dist lognormal --sigma 1.0 --min_contig 1500 "$@" && rm -rf "$SHR/$c" && mv "$SHR/$c.tmp" "$SHR/$c" \
    || { echo "FAILED shred $c" >&2; rc=1; }
}
shred mag_s0 --n50 20000 --seed 0
shred mag_s1 --n50 20000 --seed 1
shred mag_s2 --n50 20000 --seed 2
shred draft_s0 --n50 100000 --seed 0

bad=0
while IFS=$'\t' read -r path want; do
  case "$path" in shreds/*) got=$(md5sum "$SHR/${path#shreds/}" 2>/dev/null | cut -d' ' -f1)
    [ "$got" = "$want" ] || { echo "MISMATCH $path: $got != $want" >&2; bad=$((bad + 1)); } ;; esac
done < "$G/input_md5.tsv"
if [ $bad -gt 0 ]; then
  echo "FAILED: $bad shreds differ from the local ones -- the comparison would not be on identical input" >&2
  rc=1
else
  echo "[verify] all 156 shreds identical to the local ones"
fi
echo "=== job end (exit $rc) ==="; date
exit $rc
