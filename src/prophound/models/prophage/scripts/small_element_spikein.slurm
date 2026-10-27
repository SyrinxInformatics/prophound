#!/bin/bash
# the internal release checklist §3 / docs/release-scope-0.1.0.md: sub-10 kb recall with REAL small
# prophage-class elements, which geNomad truth cannot certify (its floor is 3.5 kb
# and the corpus holds 2 elements below 5 kb).
#   1. disjoint  drop candidates sharing 21-mers with any training positive or
#                training inovirus locus (the validation set must be held out)
#   2. build     implant the survivors into quiet sites of the 42-genome sample
#   3. score     every spiked genome, plus the val/test genomes carrying confident
#                inovirus loci (the in-situ arm), with r1-r3 and the ensemble,
#                keeping 5 kb fine tracks so refinement variants are free offline
#   4. eval      detection / covered50 / IoU per size bin, for every mode x seed
#                vote x refinement setting -- so it answers whichever way the
#                pending vote and refinement decisions go
# Submit from /work/$USER/mge_explorer/data_refactored:
#   sbatch scripts/small_element_spikein.slurm
#SBATCH -J small_spikein
#SBATCH -p gpu1v100
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -t 16:00:00
#SBATCH -o analysis/_slurm_logs/%x_%j.out
#SBATCH -e analysis/_slurm_logs/%x_%j.err

set -uo pipefail
echo "=== job start ==="; date; hostname; nvidia-smi --query-gpu=name --format=csv,noheader

module load miniconda
conda activate /home/$USER/.conda/envs/mge_explorer_gpu
cd /work/$USER/mge_explorer/data_refactored

CAND=scripts/spikein_small
HOSTS=scripts/eval_genome_list_sizestrat.txt
HOST_TRACKS=analysis/multiscale/ens_3000     # clean 10 kb ensemble tracks (813733/813736 era)
ELEMENTS=scripts/eval_sizestrat_elements.csv
LOCI=analysis/inovirus_overlap/loci.tsv
GENOMES=../raw/genomes                       # sibling of data_refactored/, not inside it
BEDS=beds_pos
SEQ=sequences
OUT=analysis/spikein_small
mkdir -p analysis/_slurm_logs "$OUT"

# ---- preflight (docs/environment.md: 816440/816441, 817354) ----------------
fail=0
for f in "$CAND/candidates.tsv" "$CAND/candidates.fna" "$HOSTS" "$ELEMENTS" "$LOCI" \
         scripts/small_element_spikein.py scripts/mode_recheck.py scripts/refine_sweep.py \
         scripts/compare_variants.py src/predict_genome.py splits/train.txt; do
  [ -s "$f" ] || { echo "PREFLIGHT FAIL: missing or empty $f" >&2; fail=1; }
done
grep -q SHORT_CONTIG_THRESHOLDS src/predict_genome.py \
  || { echo "PREFLIGHT FAIL: src/predict_genome.py not the current version (not pushed?)" >&2; fail=1; }
for d in "$GENOMES" "$BEDS" "$SEQ" "$HOST_TRACKS"; do
  [ -d "$d" ] || { echo "PREFLIGHT FAIL: no directory $d" >&2; fail=1; }
done
miss=0
while read -r G; do
  [ -z "$G" ] && continue
  [ -s "$GENOMES/$G.fna" ] || { miss=$((miss+1)); [ $miss -le 3 ] && echo "PREFLIGHT: no FASTA $GENOMES/$G.fna" >&2; }
  [ -s "$HOST_TRACKS/$G/windows.tsv" ] || { miss=$((miss+1)); [ $miss -le 3 ] && echo "PREFLIGHT: no track for $G" >&2; }
done < "$HOSTS"
[ $miss -eq 0 ] || { echo "PREFLIGHT FAIL: $miss host inputs missing" >&2; fail=1; }
for R in r1 r2 r3; do
  for f in runs/prophage_v2_pool_$R/checkpoints/best.pt runs/prophage_v2_pool_$R/config_snapshot.yaml; do
    [ -s "$f" ] || { echo "PREFLIGHT FAIL: missing $f" >&2; fail=1; }
  done
done
[ $fail -eq 0 ] || { echo "PREFLIGHT FAILED -- nothing run" >&2; exit 1; }
echo "[preflight] all inputs present"
# ---------------------------------------------------------------------------

# 1. hold-out check (CPU; ~2,700 training elements sketched)
if [ ! -s "$CAND/candidates_disjoint.tsv" ]; then
  python scripts/small_element_spikein.py disjoint --candidates_dir "$CAND" \
    --beds_dir "$BEDS" --sequences_dir "$SEQ" --train_split splits/train.txt \
    --inovirus_loci "$LOCI" --max_containment 0.1 || exit $?
fi

# 2. implant
if [ ! -s "$OUT/spike_truth.tsv" ]; then
  python scripts/small_element_spikein.py build --candidates_dir "$CAND" \
    --hosts "$HOSTS" --host_tracks "$HOST_TRACKS" --host_elements "$ELEMENTS" \
    --genomes_dir "$GENOMES" --inovirus_loci "$LOCI" --out_dir "$OUT" --replicates 2 || exit $?
fi

# 3. score spiked genomes and the in-situ genomes
CKPT_ALL="runs/prophage_v2_pool_r1/checkpoints/best.pt runs/prophage_v2_pool_r2/checkpoints/best.pt runs/prophage_v2_pool_r3/checkpoints/best.pt"
for R in r1 r2 r3 ens; do
  if [ "$R" = ens ]; then CKPT=$CKPT_ALL; CFG=runs/prophage_v2_pool_r1/config_snapshot.yaml
  else CKPT=runs/prophage_v2_pool_$R/checkpoints/best.pt; CFG=runs/prophage_v2_pool_$R/config_snapshot.yaml; fi
  # the ensemble carries the fine pass; single seeds only need the coarse track
  if [ "$R" = ens ]; then FINE="--refine-window-len 5000 --refine-stride 250"; else FINE=""; fi
  while read -r G; do
    [ -z "$G" ] && continue
    D="$OUT/tracks_$R/$G"
    [ -s "$D/run_metadata.json" ] && continue                 # resumable
    python -m src.predict_genome --input_fasta "$OUT/fasta/$G.fna" --outdir "$D" \
      --config "$CFG" --checkpoint $CKPT --mode default --device cuda $FINE \
      >>"$OUT/score_$R.log" 2>&1 || echo "FAILED spiked $R $G" >>"$OUT/failures.txt"
  done < "$OUT/spike_genomes.txt"
  while read -r G; do
    [ -z "$G" ] && continue
    D="$OUT/tracks_$R/$G"
    [ -s "$D/run_metadata.json" ] && continue
    python -m src.predict_genome --input_fasta "$GENOMES/$G.fna" --outdir "$D" \
      --config "$CFG" --checkpoint $CKPT --mode default --device cuda $FINE \
      >>"$OUT/score_$R.log" 2>&1 || echo "FAILED insitu $R $G" >>"$OUT/failures.txt"
  done < "$OUT/insitu_genomes.txt"
  echo "[done] scored $R: $(ls "$OUT/tracks_$R" 2>/dev/null | wc -l) genomes"; date
done

# 4. evaluate
python scripts/small_element_spikein.py eval --spike_truth "$OUT/spike_truth.tsv" \
  --insitu_truth "$OUT/insitu_truth.tsv" \
  --run r1="$OUT/tracks_r1" --run r2="$OUT/tracks_r2" --run r3="$OUT/tracks_r3" --run ens="$OUT/tracks_ens" \
  --out "$OUT/results.json"
rc=$?
if [ -s "$OUT/failures.txt" ]; then
  echo "WARNING: $(wc -l < "$OUT/failures.txt") runs failed -- see $OUT/failures.txt" >&2
  [ $rc -eq 0 ] && rc=2
fi
echo "=== job end (exit $rc) ==="; date
exit $rc
