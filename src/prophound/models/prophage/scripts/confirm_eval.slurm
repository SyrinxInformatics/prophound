#!/bin/bash
# the internal release checklist §3: independent confirmation of the shipped pipeline
# on samples that played no part in tuning.
#   - the original 18-genome / 52-element sample (eval_genome_list.txt): modes,
#     short contigs (scan, now default), impute (now default; EDL933 is in it),
#     and the fine-pass refinement at 0.6 / 0.8 / 0.95 -- the 0.8 result decides
#     whether refinement goes on by default (user, 2026-09-28)
#   - 20 test-split genomes with no geNomad provirus (eval_phagefree_list.txt):
#     calls per genome and per Mb on geNomad-negative sequence, per mode
# Scores once per (seed | ensemble) with 5 kb fine tracks; everything else is
# rebuilt offline from the saved tracks (refine_sweep.py, mode_recheck.py).
# Submit from /work/$USER/mge_explorer/data_refactored:
#   sbatch scripts/confirm_eval.slurm
#SBATCH -J confirm_eval
#SBATCH -p gpu1v100
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -t 12:00:00
#SBATCH -o analysis/_slurm_logs/%x_%j.out
#SBATCH -e analysis/_slurm_logs/%x_%j.err

set -uo pipefail
echo "=== job start ==="; date; hostname; nvidia-smi --query-gpu=name --format=csv,noheader

module load miniconda
conda activate /home/$USER/.conda/envs/mge_explorer_gpu
cd /work/$USER/mge_explorer/data_refactored

LIST18=scripts/eval_genome_list.txt
ELEMENTS18=scripts/eval_confirm18_elements.csv
LISTPF=scripts/eval_phagefree_list.txt
GENOMES=../raw/genomes        # sibling of data_refactored, not inside it (docs/environment.md)
mkdir -p analysis/_slurm_logs analysis/confirm18 analysis/phagefree

# ---- preflight (docs/environment.md: 816440/816441, 817354) ----------------
fail=0
for f in "$LIST18" "$ELEMENTS18" "$LISTPF" src/predict_genome.py src/multiscale.py \
         scripts/refine_sweep.py scripts/compare_variants.py scripts/mode_recheck.py; do
  [ -s "$f" ] || { echo "PREFLIGHT FAIL: missing or empty $f" >&2; fail=1; }
done
for pat in SHORT_CONTIG_THRESHOLDS DEFAULT_REFINE_THRESHOLD 'DEFAULT_SPLIT_MIN_GAP_BP = 4000'; do
  grep -q "$pat" src/predict_genome.py \
    || { echo "PREFLIGHT FAIL: src/predict_genome.py lacks '$pat' (not pushed?)" >&2; fail=1; }
done
grep -q "def ramp_alpha" src/multiscale.py \
  || { echo "PREFLIGHT FAIL: src/multiscale.py predates --edge-alpha-ramp (not pushed?)" >&2; fail=1; }
grep -qF 'int(meta["window_len"])]' scripts/refine_sweep.py \
  || { echo "PREFLIGHT FAIL: scripts/refine_sweep.py lacks the short-window filter (not pushed?)" >&2; fail=1; }
[ -d "$GENOMES" ] || { echo "PREFLIGHT FAIL: no genome dir $GENOMES" >&2; fail=1; }
for R in r1 r2 r3; do
  for f in runs/prophage_v2_pool_$R/checkpoints/best.pt runs/prophage_v2_pool_$R/config_snapshot.yaml; do
    [ -s "$f" ] || { echo "PREFLIGHT FAIL: missing $f" >&2; fail=1; }
  done
done
for L in "$LIST18" "$LISTPF"; do
  [ -s "$L" ] || continue
  miss=0
  while read -r G; do
    [ -z "$G" ] && continue
    [ -s "$GENOMES/$G.fna" ] || { miss=$((miss+1)); [ $miss -le 3 ] && echo "PREFLIGHT: no FASTA $GENOMES/$G.fna" >&2; }
  done < "$L"
  [ $miss -eq 0 ] || { echo "PREFLIGHT FAIL: $miss genomes in $L have no FASTA" >&2; fail=1; }
  echo "[preflight] $L: $(wc -l < "$L") genomes"
done
[ $fail -eq 0 ] || { echo "PREFLIGHT FAILED -- nothing run" >&2; exit 1; }
echo "[preflight] all inputs present"
# ---------------------------------------------------------------------------

CKPT_ALL="runs/prophage_v2_pool_r1/checkpoints/best.pt runs/prophage_v2_pool_r2/checkpoints/best.pt runs/prophage_v2_pool_r3/checkpoints/best.pt"
RUNS18=""; RUNS18_5K=""; RUNSPF=""
for R in r1 r2 r3 ens; do
  if [ "$R" = ens ]; then CKPT=$CKPT_ALL; CFG=runs/prophage_v2_pool_r1/config_snapshot.yaml
  else CKPT=runs/prophage_v2_pool_$R/checkpoints/best.pt; CFG=runs/prophage_v2_pool_$R/config_snapshot.yaml; fi
  for S in confirm18 phagefree; do
    if [ "$S" = confirm18 ]; then L=$LIST18; else L=$LISTPF; fi
    OUT=analysis/$S
    while read -r G; do
      [ -z "$G" ] && continue
      D="$OUT/$R/$G"
      # resumable; keyed on metadata since call-free genomes write no fine track
      [ -s "$D/run_metadata.json" ] && continue
      # shipped defaults (impute, scan, coverage+merge) plus the tight-boundaries
      # preset, so regions.tsv is that preset's output; the sweeps rebuild the rest
      python -m src.predict_genome \
        --input_fasta "$GENOMES/$G.fna" --outdir "$D" \
        --config "$CFG" --checkpoint $CKPT \
        --mode default --device cuda \
        --refine-window-len 5000 --refine-stride 250 \
        --edge-alpha 0.6 --edge-alpha-ramp 15000 60000 \
        >>"$OUT/${R}.log" 2>&1 || echo "FAILED $S $R $G" >>"$OUT/failures.txt"
    done < "$L"
    echo "[done] $S $R: $(ls "$OUT/$R" 2>/dev/null | wc -l) genomes"; date
  done
  RUNS18="$RUNS18 --run $R=analysis/confirm18/$R"
  RUNS18_5K="$RUNS18_5K --run ${R}_5000=analysis/confirm18/$R"
  RUNSPF="$RUNSPF --run $R=analysis/phagefree/$R"
done

rc=0
python scripts/mode_recheck.py --elements "$ELEMENTS18" --genomes_dir "$GENOMES" $RUNS18 \
  --out analysis/confirm18/modes.json || rc=$?
python scripts/mode_recheck.py --genomes_dir "$GENOMES" $RUNSPF \
  --out analysis/phagefree/modes.json || rc=$?
python scripts/refine_sweep.py --elements "$ELEMENTS18" \
  $RUNS18_5K \
  --refine_thresholds 0.6,0.8,0.95 --alphas 0.6 --split_gaps 4000 --split_max_probs 0.05 \
  --ramps 15000:60000 --out analysis/confirm18/refine_sweep.json || rc=$?
for S in confirm18 phagefree; do
  if [ -s "analysis/$S/failures.txt" ]; then
    echo "WARNING: $(wc -l < "analysis/$S/failures.txt") genome runs failed -- see analysis/$S/failures.txt" >&2
    [ $rc -eq 0 ] && rc=2
  fi
done
echo "=== job end (exit $rc) ==="; date
exit $rc
