#!/bin/bash
# docs/lockbox2.md: run the frozen prophound on the second locked test set. GPU. Run ONCE.
# Each mode is a separate run of the shipped command (`prophound predict --mode <m>`, the
# bundled model, every default), so regions.detailed.tsv is exactly what a user gets:
# seed vote, short-contig scan, edge snap, N-gap filter, refinement, att search,
# composition edges and edge ranges, with nothing rebuilt offline. Default also writes
# per-checkpoint scores, for diagnostics only. Scoring is local and separate, also run
# once: scripts/lockbox_score.py --boundaries.
#
# The code is a `git archive` of the freeze commit, extracted to $CODE (not this
# data_refactored tree, whose src/ predates the package layout). The preflight refuses
# to run unless every file of it matches analysis/lockbox2/frozen_md5.txt, that file's
# SHA-256 is the one below, the genome list is the preregistered one, every FASTA is
# present, and the example genome gives its expected calls on this GPU.
# Push first (docs/lockbox2.md Status): the archive, frozen_md5.txt and this launcher.
# Submit from /work/$USER/mge_explorer/data_refactored:
#   sbatch scripts/lockbox2_run.slurm
#SBATCH -J lockbox2_run
#SBATCH -p gpu1v100
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -t 2-00:00:00
#SBATCH -o analysis/_slurm_logs/%x_%j.out
#SBATCH -e analysis/_slurm_logs/%x_%j.err

set -uo pipefail
echo "=== job start ==="; date; hostname; nvidia-smi --query-gpu=name --format=csv,noheader

module load miniconda
conda activate /home/$USER/.conda/envs/mge_explorer_gpu
cd /work/$USER/mge_explorer/data_refactored

FREEZE=0d71e2f                                                        # the freeze commit (docs/lockbox2.md Status)
CODE=/work/$USER/mge_explorer/prophound_$FREEZE
MD5S=analysis/lockbox2/frozen_md5.txt
MD5S_SHA=21e7485299d65f75a47573c9e92a53c092f86756d0fb2703ae95c8e2396d7baf   # sha256 of frozen_md5.txt (745 files under src/)
LB=analysis/lockbox2
LIST=$LB/lockbox_genome_list.txt
LIST_SHA=c6f57869ceaa8ce08a6b4750859a8cf683dd15195de4a23f6d31d0398974a4d0   # docs/lockbox2.md Status
FNA=../raw/lockbox2
OUT=$LB/prophound
MODES="default exploratory conservative"
export PYTHONPATH="$CODE/src"
PREDICT=(python -c "from prophound.cli import main; main()" predict)
mkdir -p analysis/_slurm_logs "$OUT"

# ---- preflight (docs/environment.md) ---------------------------------------
fail=0
case "$FREEZE$MD5S_SHA" in *__*) echo "PREFLIGHT FAIL: freeze commit / manifest hash not filled in" >&2; exit 1;; esac
[ -s "$LIST" ] || { echo "PREFLIGHT FAIL: missing or empty $LIST" >&2; fail=1; }
got=$(sha256sum "$LIST" 2>/dev/null | cut -d' ' -f1)
[ "$got" = "$LIST_SHA" ] || { echo "PREFLIGHT FAIL: $LIST sha256 $got != preregistered $LIST_SHA" >&2; fail=1; }
got=$(sha256sum "$MD5S" 2>/dev/null | cut -d' ' -f1)
[ "$got" = "$MD5S_SHA" ] || { echo "PREFLIGHT FAIL: $MD5S sha256 $got != frozen $MD5S_SHA" >&2; fail=1; }
[ -d "$CODE/src/prophound" ] || { echo "PREFLIGHT FAIL: no frozen code at $CODE" >&2; fail=1; }
if [ $fail -eq 0 ]; then
  (cd "$CODE" && md5sum -c --quiet "$OLDPWD/$MD5S") \
    || { echo "PREFLIGHT FAIL: frozen code at $CODE differs from $MD5S" >&2; fail=1; }
  extra=$(cd "$CODE" && find src -name '*.py' | sort | comm -23 - <(awk '{print $2}' "$OLDPWD/$MD5S" | sort))
  [ -z "$extra" ] || { echo "PREFLIGHT FAIL: modules not in $MD5S: $extra" >&2; fail=1; }
fi
python -c "import prophound.cli, sys; assert prophound.cli.__file__.startswith('$CODE'), prophound.cli.__file__" \
  || { echo "PREFLIGHT FAIL: prophound does not import from $CODE" >&2; fail=1; }
while read -r G; do
  [ -s "$FNA/$G.fna" ] || { echo "PREFLIGHT FAIL: missing $FNA/$G.fna" >&2; fail=1; }
done < "$LIST"
[ $fail -eq 0 ] || { echo "PREFLIGHT FAILED -- nothing run" >&2; exit 1; }
# the example genome must give its expected calls on this device (calls, not scores: docs/benchmarks.md)
SMOKE="$OUT/_smoke"; rm -rf "$SMOKE"
"${PREDICT[@]}" --input "$CODE/examples/hin_rd_kw20_region.fna" --outdir "$SMOKE" --device cuda -q \
  && python - "$SMOKE/regions.detailed.tsv" "$CODE/examples/expected_calls.tsv" <<'EOF' \
  || { echo "PREFLIGHT FAILED -- example genome calls differ on this device; nothing run" >&2; exit 1; }
import sys, pandas as pd
cols = ["call_id", "contig_id", "start", "end", "boundary", "contig_edge", "seed_support"]
got, exp = (pd.read_csv(p, sep="\t")[cols] for p in sys.argv[1:])
pd.testing.assert_frame_equal(got, exp)
print(f"[preflight] example genome: {len(got)} calls as expected")
EOF
echo "[preflight] list, frozen code ($FREEZE), example genome and all $(wc -l < "$LIST") FASTAs verified"
# ---------------------------------------------------------------------------

for m in $MODES; do
  echo "[$m]"; date
  EXTRA=(); [ "$m" = default ] && EXTRA=(--write-seed-scores)
  while read -r G; do
    D="$OUT/$m/$G"
    [ -s "$D/run_metadata.json" ] && continue                        # resumable
    s=$(date +%s)
    "${PREDICT[@]}" --input "$FNA/$G.fna" --outdir "$D" --mode "$m" --device cuda -q "${EXTRA[@]}" \
      >>"$OUT/run.log" 2>&1 \
      && echo -e "$m\t$G\t$(( $(date +%s) - s ))" >> "$OUT/wall_seconds.tsv" \
      || echo "FAILED $m $G" >> "$OUT/failures.txt"
  done < "$LIST"
done
N=$(wc -l < "$LIST")
echo "[done] $(ls "$OUT"/{default,exploratory,conservative}/*/run_metadata.json 2>/dev/null | wc -l) of $(( N * 3 )) genome runs"; date

rc=0
if [ -s "$OUT/failures.txt" ]; then
  echo "WARNING: $(wc -l < "$OUT/failures.txt") genome runs failed -- see $OUT/failures.txt" >&2
  rc=2
fi
echo "=== job end (exit $rc) ==="; date
exit $rc
