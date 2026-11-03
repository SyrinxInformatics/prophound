# Archived prophage scripts

Early development scripts, kept for the record and not maintained: some import modules or read paths that no longer exist. Nothing in the package, the tests, the launchers or the other scripts uses them. Each lived in `prophound/models/prophage/<old dir>/`.

| file | old dir | what it was |
|---|---|---|
| `1_to_1.py` | `scripts/` | early positive/negative window balancing, run at import time |
| `prep_labels_and_splits_contigs.py` | `scripts/` | early contig-level labels and splits |
| `select_genomad_bias_examples.py` | `scripts/` | picks examples for an early audit of geNomad's labeling bias |
| `generate_chunks.py` | `scripts/` | splits a file list into job chunks (hard-coded HPCC path) |
| `analyze_postprocessing_disagreements.py` | `scripts/` | pre-v1 post-processing comparison |
| `compare_postprocessing_candidates.py` | `scripts/` | pre-v1 post-processing comparison |
| `analyze_genomad_region_bias.py` | `scripts/` | geNomad-bias audit, region level |
| `analyze_genomad_bias.py` | `scripts/` | geNomad-bias audit, window level |
| `analyze_operating_modes.py` | `scripts/` | pre-vote operating-mode table |
| `analyze_region_boundary_errors.py` | `scripts/` | pre-v1 boundary-error analysis |
| `analyze_region_coherence.py` | `scripts/` | pre-v1 region-coherence analysis |
| `analyze_threshold_calibration.py` | `scripts/` | pre-v1 threshold calibration |
| `analyze_prediction_protein_join.py` | `scripts/` | joins window predictions with protein corroboration (superset of the `src/` copy) |
| `make_splits.py` | `scripts/` | duplicate of `src/make_splits.py` (same logic, different variable names); the `src/` copy is kept |
| `protein_corroboration.py` | `src/` | per-window gene-overlap features from geNomad annotations |
| `prophage_protein_corroboration.slurm` | `src/` | launcher for the above |
| `analyze_prediction_protein_join.py` | `src/` | older copy of the `scripts/` file above |

The training-data builders that made the corpus stay in place: `src/make_manifest_from_bed.py`, `src/make_splits.py`, `scripts/bed_generator.py`, `scripts/generate_negatives.py`, `scripts/make_hard_negatives.py` and `scripts/make_manifest.py`.

## The boundary-regression head

An auxiliary head that predicted each window's distance to an element boundary was tried and shelved; no shipped checkpoint has it, and `load_cfg` rejects a config with `model.predict_boundary: true`. These files trained or fed that head:

| file | old dir | what it was |
|---|---|---|
| `add_boundary_distance.py` | `scripts/` | adds the `boundary_dist_bp` target column to the window manifest |
| `add_boundary_distance.slurm` | `scripts/` | launcher for the above |
| `prophage_v1_positional_boundary_gpu.slurm` | `src/` | training launcher, boundary head on |
| `prophage_v1_positional_boundary_lowweight_gpu.slurm` | `src/` | the same at a lower boundary-loss weight |
| `configs/prophage_v1_positional_boundary.yaml` | `src/configs/` | its config |
| `configs/prophage_v1_positional_boundary_lowweight.yaml` | `src/configs/` | its config |
