#!/bin/bash
# G5 fragmented-input comparison (docs/claims-and-evidence.md "G5 preregistration"),
# Jaeger arm. CPU. Jaeger 1.26.5 (installed by the user in mge_explorer_gpu,
# 2026-09-30) on the 39 complete LAMBDA chromosomes and the four shred
# conditions from g5_shred.slurm, with both models (57341_1.5M primary, bundled
# `default` secondary) and the preregistered flags (`-p --lc 1500 --cpu`, plus
# `--fsize 3077` for 57341_1.5M, its training length: amendment 2026-09-30), then
# collects every call into analysis/g5/jaeger_calls.tsv, which is all that needs
# pulling. Replaces the local runner scripts/g5_jaeger.sh (same models, flags and
# output layout). Resumable. Submit after the shred job, from
# /work/$USER/mge_explorer/data_refactored:
#   sbatch --dependency=afterok:822964 scripts/g5_jaeger.slurm
# (or plain `sbatch scripts/g5_jaeger.slurm` once 822964 has finished cleanly).
#SBATCH -J g5_jaeger
#SBATCH -p compute1
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -t 1-00:00:00
#SBATCH -o analysis/_slurm_logs/%x_%j.out
#SBATCH -e analysis/_slurm_logs/%x_%j.err

set -uo pipefail
echo "=== job start ==="; date; hostname

module load miniconda
conda activate /home/$USER/.conda/envs/mge_explorer_gpu
cd /work/$USER/mge_explorer/data_refactored

G=analysis/g5
FNA=../raw/lambda
SHR=../raw/g5_shreds
RAW=../g5_jaeger                 # raw outputs stay on ARC
MODELDIR=/work/$USER/mge_explorer/jaeger_models
PRIMARY=jaeger_57341_1.5M_fragment
# label:model:extra flags (the primary model was trained on 3077 nt fragments; Jaeger warns against its 2000 default)
MODELS="57341_1.5M:$PRIMARY:--fsize=3077 default:default:"
CONDS="complete mag_s0 mag_s1 mag_s2 draft_s0"
WORKERS=8; P=$(( ${SLURM_CPUS_PER_TASK:-8} / WORKERS )); [ "$P" -ge 1 ] || P=1
mkdir -p analysis/_slurm_logs "$RAW"
rm -f "$RAW/failures.txt" "$G/jaeger_failures.txt"   # per job; a retry must not inherit old failures

# ---- preflight (docs/environment.md) ---------------------------------------
fail=0
for f in "$G/genome_list.txt" "$G/input_md5.tsv" scripts/g5_score.py; do
  [ -s "$f" ] || { echo "PREFLIGHT FAIL: missing or empty $f" >&2; fail=1; }
done
command -v jaeger >/dev/null || { echo "PREFLIGHT FAIL: jaeger not on PATH in mge_explorer_gpu" >&2; fail=1; }
# the prophage step imports pyrodigal_gv and, if it is missing, logs an ERROR but exits 0 (822991)
python -c "import pyrodigal_gv" 2>/dev/null \
  || { echo "PREFLIGHT FAIL: pyrodigal_gv not importable (pip install 'pyrodigal-gv>=0.3.0')" >&2; fail=1; }
v=$(jaeger --version 2>&1 | grep -o '[0-9]\+\.[0-9]\+\.[0-9]\+' | tail -1)
[ "$v" = "1.26.5" ] || { echo "PREFLIGHT FAIL: jaeger --version gave '$(jaeger --version 2>&1 | tail -1)', expected 1.26.5" >&2; fail=1; }
# every input, complete and shredded, must be byte-identical to the local files
bad=0
while IFS=$'\t' read -r path want; do
  case "$path" in fasta/*) f="$FNA/${path#fasta/}" ;; shreds/*) f="$SHR/${path#shreds/}" ;; *) continue ;; esac
  got=$(md5sum "$f" 2>/dev/null | cut -d' ' -f1)
  [ "$got" = "$want" ] || { echo "PREFLIGHT FAIL: $f md5 '$got' != $want" >&2; bad=$((bad + 1)); }
done < "$G/input_md5.tsv"
[ $bad -eq 0 ] || { echo "PREFLIGHT FAIL: $bad inputs differ (has g5_shred.slurm finished?)" >&2; fail=1; }
[ $fail -eq 0 ] || { echo "PREFLIGHT FAILED -- nothing run" >&2; exit 1; }
echo "[preflight] jaeger $v; all 195 inputs match input_md5.tsv"
# ---------------------------------------------------------------------------

# the primary model is a separate download (1.26 no longer bundles it). A one-genome
# smoke run shows whether it is installed and loads; download only if it does not.
T=$(head -1 "$G/genome_list.txt")
# a clean exit is not enough: Jaeger exits 0 when its prophage step fails, so require the
# prophage table and no ERROR line
smoke() { rm -rf "$RAW/_smoke"; jaeger predict -m "$PRIMARY" -i "$FNA/$T.fna" -o "$RAW/_smoke" -p --lc 1500 --cpu \
            --fsize 3077 --workers "$WORKERS" -f > "$RAW/_smoke.log" 2>&1 \
          && ls "$RAW"/_smoke/*/*prophages/prophages_jaeger.tsv >/dev/null 2>&1 \
          && ! grep -q ' ERROR ' "$RAW/_smoke.log"; }
if ! smoke; then
  echo "[model] $PRIMARY not usable yet; downloading to $MODELDIR"
  jaeger download --path "$MODELDIR" --model_name "$PRIMARY" \
    || { echo "FAILED: could not download $PRIMARY" >&2; exit 1; }
  smoke || { echo "FAILED: smoke run with $PRIMARY after download -- see $RAW/_smoke.log" >&2; exit 1; }
fi
rm -rf "$RAW/_smoke"
echo "[model] $PRIMARY loads and runs"

export RAW WORKERS
for m in $MODELS; do
  IFS=: read -r label model extra <<< "$m"; extra=${extra//=/ }
  for c in $CONDS; do
    src=$([ "$c" = complete ] && echo "$FNA" || echo "$SHR/$c")
    echo "[$label $c]"; date
    sed "s|^|$src/|; s|$|.fna|" "$G/genome_list.txt" | L=$label M=$model X="$extra" C=$c xargs -n1 -P"$P" -I{} bash -c '
      f="{}"; b=$(basename "$f" .fna); o="$RAW/$L/$C/$b"
      ls "$o"/*/*.log >/dev/null 2>&1 && ls "$o"/*/"$b"*.tsv >/dev/null 2>&1 && exit 0     # resumable
      rm -rf "$o"; mkdir -p "$(dirname "$o")"; s=$(date +%s)
      jaeger predict -m "$M" -i "$f" -o "$o" -p --lc 1500 --cpu $X --workers "$WORKERS" -f > "$o.log" 2>&1 \
        && ! grep -q " ERROR " "$o.log" \
        && echo -e "$L\t$C\t$b\t$(( $(date +%s) - s ))" >> "$RAW/wall_seconds.tsv" \
        || echo "FAILED $L $C $b" >> "$RAW/failures.txt"
    '
  done
done
echo "[jaeger] $(ls -d "$RAW"/*/*/*/ 2>/dev/null | wc -l) of 390 runs present"; date

python scripts/g5_score.py collect-jaeger --root "$RAW" --out "$G/jaeger_calls.tsv"
rc=$?
cp "$RAW/wall_seconds.tsv" "$G/jaeger_wall_seconds.tsv" 2>/dev/null
if [ -s "$RAW/failures.txt" ]; then
  cp "$RAW/failures.txt" "$G/jaeger_failures.txt"
  echo "WARNING: $(wc -l < "$RAW/failures.txt") Jaeger runs failed -- see $G/jaeger_failures.txt" >&2
  [ $rc -eq 0 ] && rc=2
fi
echo "=== job end (exit $rc) ==="; date
exit $rc
