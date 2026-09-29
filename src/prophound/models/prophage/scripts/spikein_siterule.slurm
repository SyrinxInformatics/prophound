#!/bin/bash
# docs/research-integrity.md C7: the sub-10 kb spike-in (821410) with implant sites
# chosen without model scores (--site_rule genomad: geNomad elements and spacing
# only). Same hold-out-checked elements, hosts, seed and spacing. Only the spiked
# genomes are re-scored; the in-situ arm reuses 821410's tracks by symlink.
# Scored with the shipped modes, refinement off.
# Submit from /work/$USER/mge_explorer/data_refactored:
#   sbatch scripts/spikein_siterule.slurm
#SBATCH -J spikein_siterule
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

CAND=scripts/spikein_small
HOSTS=scripts/eval_genome_list_sizestrat.txt
HOST_TRACKS=analysis/multiscale/ens_3000     # used for scanned-bin coverage only under --site_rule genomad
ELEMENTS=scripts/eval_sizestrat_elements.csv
LOCI=analysis/inovirus_overlap/loci.tsv
GENOMES=../raw/genomes                       # sibling of data_refactored/, not inside it
OLD=analysis/spikein_small                   # 821410
OUT=analysis/spikein_siterule
mkdir -p analysis/_slurm_logs "$OUT"

# ---- preflight (docs/environment.md) ---------------------------------------
fail=0
for f in "$CAND/candidates_disjoint.tsv" "$CAND/candidates.fna" "$HOSTS" "$ELEMENTS" "$LOCI" \
         "$OLD/insitu_genomes.txt" scripts/small_element_spikein.py scripts/mode_recheck.py src/predict_genome.py; do
  [ -s "$f" ] || { echo "PREFLIGHT FAIL: missing or empty $f" >&2; fail=1; }
done
grep -q -- '--site_rule' scripts/small_element_spikein.py \
  || { echo "PREFLIGHT FAIL: scripts/small_element_spikein.py lacks --site_rule (not pushed?)" >&2; fail=1; }
grep -q MODE_MIN_SEED_SUPPORT src/predict_genome.py \
  || { echo "PREFLIGHT FAIL: src/predict_genome.py predates the seed vote" >&2; fail=1; }
n=$(ls "$HOST_TRACKS"/*/windows.tsv 2>/dev/null | wc -l)
[ "$n" -eq 42 ] || { echo "PREFLIGHT FAIL: $HOST_TRACKS has $n/42 genomes" >&2; fail=1; }
for R in r1 r2 r3 ens; do
  n=$(ls -d "$OLD/tracks_$R"/*/ 2>/dev/null | wc -l)
  [ "$n" -ge 25 ] || { echo "PREFLIGHT FAIL: $OLD/tracks_$R has $n genomes (need the in-situ ones)" >&2; fail=1; }
done
for R in r1 r2 r3; do
  for f in runs/prophage_v2_pool_$R/checkpoints/best.pt runs/prophage_v2_pool_$R/config_snapshot.yaml; do
    [ -s "$f" ] || { echo "PREFLIGHT FAIL: missing $f" >&2; fail=1; }
  done
done
[ -d "$GENOMES" ] || { echo "PREFLIGHT FAIL: no genome dir $GENOMES" >&2; fail=1; }
[ $fail -eq 0 ] || { echo "PREFLIGHT FAILED -- nothing run" >&2; exit 1; }
echo "[preflight] all inputs present"
# ---------------------------------------------------------------------------

if [ ! -s "$OUT/spike_truth.tsv" ]; then
  python scripts/small_element_spikein.py build --candidates_dir "$CAND" \
    --hosts "$HOSTS" --host_tracks "$HOST_TRACKS" --host_elements "$ELEMENTS" \
    --genomes_dir "$GENOMES" --inovirus_loci "$LOCI" --out_dir "$OUT" --replicates 2 \
    --site_rule genomad || exit $?
fi

CKPT_ALL="runs/prophage_v2_pool_r1/checkpoints/best.pt runs/prophage_v2_pool_r2/checkpoints/best.pt runs/prophage_v2_pool_r3/checkpoints/best.pt"
for R in r1 r2 r3 ens; do
  if [ "$R" = ens ]; then CKPT=$CKPT_ALL; CFG=runs/prophage_v2_pool_r1/config_snapshot.yaml
  else CKPT=runs/prophage_v2_pool_$R/checkpoints/best.pt; CFG=runs/prophage_v2_pool_$R/config_snapshot.yaml; fi
  mkdir -p "$OUT/tracks_$R"
  while read -r G; do
    [ -z "$G" ] && continue
    D="$OUT/tracks_$R/$G"
    [ -s "$D/run_metadata.json" ] && continue                 # resumable
    python -m src.predict_genome --input_fasta "$OUT/fasta/$G.fna" --outdir "$D" \
      --config "$CFG" --checkpoint $CKPT --mode default --device cuda \
      >>"$OUT/score_$R.log" 2>&1 || echo "FAILED spiked $R $G" >>"$OUT/failures.txt"
  done < "$OUT/spike_genomes.txt"
  while read -r G; do                                          # in-situ: 821410's tracks, unchanged
    [ -z "$G" ] && continue
    [ -e "$OUT/tracks_$R/$G" ] || ln -s "$(readlink -f "$OLD/tracks_$R/$G")" "$OUT/tracks_$R/$G"
  done < "$OUT/insitu_genomes.txt"
  echo "[done] scored $R: $(ls "$OUT/tracks_$R" 2>/dev/null | wc -l) genomes"; date
done

python scripts/small_element_spikein.py eval --spike_truth "$OUT/spike_truth.tsv" \
  --insitu_truth "$OUT/insitu_truth.tsv" --refine off \
  --run r1="$OUT/tracks_r1" --run r2="$OUT/tracks_r2" --run r3="$OUT/tracks_r3" --run ens="$OUT/tracks_ens" \
  --out "$OUT/results.json"
rc=$?
if [ -s "$OUT/failures.txt" ]; then
  echo "WARNING: $(wc -l < "$OUT/failures.txt") runs failed -- see $OUT/failures.txt" >&2
  [ $rc -eq 0 ] && rc=2
fi
echo "=== job end (exit $rc) ==="; date
exit $rc
