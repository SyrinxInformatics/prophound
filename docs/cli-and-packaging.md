# CLI & Packaging

Living doc — update in place, don't append a history here (that's `SESSIONS.md`'s job).

## Current CLI shape

`mge-detect predict` (entrypoint `src/mge_detect/mge_detect.py`, dispatching into `src/mge_detect/models/prophage/src/predict_genome.py`):

```
mge-detect predict --input <fasta> --outdir <dir> --config <config.yaml> --checkpoint <best.pt> \
  [--mode exploratory|default|conservative] [--threshold F] [--seed-threshold F] [--extend-threshold F] \
  [--device auto|cpu|cuda|rocm|mps|xpu|intel] [--batch-size N] [--num-workers N] [--num-tiles N]
```

Outputs: `windows.tsv`, `regions.tsv`, `run_metadata.json` (+ `bins.tsv` when `--num-tiles > 1`). This matches `src/mge_detect/models/prophage/notes/V1_INFERENCE_CONTRACT.md`.

`--config`/`--checkpoint` are required, with no bundled default weights — model-weight bundling is an open packaging question upstream (`notes/V1_PACKAGE_BOUNDARY.md`) that hasn't been decided, so this doesn't guess at an answer.

`--num-tiles` (default 1 = original frozen v1 behavior, unchanged) controls the multi-tile boundary-refinement mode — see `docs/boundary-refinement.md` for what that does and why.

Console script: `[project.scripts] mge-detect = "mge_detect.mge_detect:main"` in `pyproject.toml` (added 2026-09-08, replacing a fragile `python3 src/mge_detect/mge_detect.py` pixi task that broke package-relative imports when run as a bare script).

## Package structure

`src/mge_detect/models/prophage/src/` is a real importable subpackage of `mge_detect` (`__init__.py` added at `models/`, `models/prophage/`, and `models/prophage/src/`; internal imports changed from `from src.X import Y` to relative `from .X import Y`).

**Divergence from the HPCC copy**: the HPCC's own copy of this code (under `/work/pke835/mge_explorer/data_refactored/src/`) still uses the old `from src.X import Y` style and is run there as `python -m src.train_prophage --config ...` from that repo's root. This local copy has diverged in import style to fit inside the `mge_detect` package. If code is ever synced back to the HPCC (e.g. an updated `train_prophage.py`), reconcile the import style deliberately rather than copying blind.

## Dependency fix (2026-09-08)

`pyproject.toml` previously declared `tensorflow` (nothing in the codebase imports it) and was missing `torch`, which every model file actually needs. Fixed: now declares `pytorch`, `pandas`, `pyfaidx`, `tqdm`, `pyyaml` (plus the pre-existing `biopython`, `scikit-learn`, `pytest`, `sphinx`, `mypy`).

**Not yet fixed**: `metadata/conda/mge_detect_05_11_2025.yaml` is a separate, older conda env snapshot that still lists `tensorflow`/`keras`, not `pytorch`. Needs the same correction whenever that file is next touched — not yet done as of 2026-09-08.

## Positional-pooling experiment (found, evaluated, not yet promoted)

`runs/prophage_v1_positional/` (config + checkpoint pulled from HPCC to local scratch, not into the repo) was scaffolded by a prior agent, actually trained to completion (12 epochs, early-stopped) and test-evaluated — but the results were never written up, so they'd been sitting unused. Evaluated locally 2026-09-08 by running `analyze_boundary_slices.py` against its `test_predictions.csv` at the Default threshold:

| Metric | v1 (best of 3 replicates) | v1_positional |
|---|---|---|
| Test AUPRC | 0.5707 | 0.6063 |
| Test AUROC | 0.9451 | 0.9519 |
| `positive_core` recall | 0.6677 | 0.7053 |
| `positive_edge` recall | 0.3803 | 0.3812 (unchanged) |

Real overall improvement, but it does **not** fix the edge/boundary weakness specifically — the architecture change (adaptive pooling) helps general discrimination, not boundary localization. Caveat: `pos_weight` and `amp` also differ between the v1 and v1_positional configs, so the gain isn't a clean single-variable ablation.

Checkpoint confirmed to exist on HPCC: `runs/prophage_v1_positional/checkpoints/best.pt` (16M). This is the checkpoint currently being used for all local smoke-testing (see `docs/boundary-refinement.md`).

**Not yet decided**: whether to promote this as the new default checkpoint/config. No `RESULTS.md` has been written for this run yet (unlike `prophage_v0`/`prophage_v1`, which have one under `runs/<name>/RESULTS.md` on the HPCC).
