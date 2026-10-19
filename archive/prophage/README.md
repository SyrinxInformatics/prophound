# Archived prophage scripts

Moved here on 2026-10-19 by the technical debt audit (`docs/tech-debt-audit.md` X1, user decision 2). Nothing in `docs/`, `CLAUDE.md`, a launcher, a test or another script refers to these files. Where they are mentioned at all, it is in the frozen `src/prophound/models/prophage/notes/`, which cite them by their old paths. They are kept for that record and are not maintained: some import modules or read paths that no longer exist.

Old location: `src/prophound/models/prophage/<dir>/<file>`.

| file | old dir | what it was |
|---|---|---|
| `1_to_1.py` | `scripts/` | early positive/negative window balancing, run at import time |
| `prep_labels_and_splits_contigs.py` | `scripts/` | early contig-level labels and splits |
| `select_genomad_bias_examples.py` | `scripts/` | picks examples for the geNomad-bias audit (`notes/V1_GENOMAD_BIAS_AUDIT.md`) |
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

The training-data builders that made the corpus stay in place: `src/make_manifest_from_bed.py` (`notes/RUNBOOK.md` step 2), `src/make_splits.py`, `scripts/bed_generator.py`, `scripts/generate_negatives.py`, `scripts/make_hard_negatives.py` and `scripts/make_manifest.py`.
