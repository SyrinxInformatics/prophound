# CLI & Packaging

Living doc — update in place, don't append a history here (that's `SESSIONS.md`'s job).

## Current CLI shape

`mge-detect predict` (entrypoint `src/mge_detect/mge_detect.py`, dispatching into `src/mge_detect/models/prophage/src/predict_genome.py`):

```
mge-detect predict --input <fasta> --outdir <dir> --config <config.yaml> --checkpoint <best.pt> \
  [--mode exploratory|default|conservative] [--threshold F] [--seed-threshold F] [--extend-threshold F] \
  [--device auto|cpu|cuda|rocm|mps|xpu|intel] [--batch-size N] [--num-workers N] [--num-tiles N]
  [--region-builder strict|coverage] [--coverage-frac F]
```

Outputs: `windows.tsv`, `regions.tsv`, `run_metadata.json` (+ `bins.tsv` when `--num-tiles > 1`). This matches `src/mge_detect/models/prophage/notes/V1_INFERENCE_CONTRACT.md`.

`--config`/`--checkpoint` are required, with no bundled default weights — model-weight bundling is an open packaging question upstream (`notes/V1_PACKAGE_BOUNDARY.md`) that hasn't been decided, so this doesn't guess at an answer.

`--region-builder` (added 2026-09-17; default `strict` = the frozen v1 builder, byte-identical output to before the flag existed — verified) selects how window scores become regions. `coverage` uses `predict_genome.build_region_table_coverage`: a stride bin is inside if ≥ `--coverage-frac` (default 0.3) of its covering windows individually score ≥ the extend threshold, instead of the strict builder's "mean of all covering windows ≥ extend". On the 18-genome eval set that lifts region recall 76.9%→88.5% (frac 0.3) / 94.2% (frac ≤0.1) and median IoU 0.56→0.75, at 2–4× the calls per genome; `run_metadata.json` records `region_builder: coverage_fraction` and `coverage_frac`. Not combinable with `--num-tiles > 1` (hard error). **Still opt-in, not the shipped default**: frac was chosen on the same 52 regions it was scored on — pending a held-out, size-stratified check (`docs/multi-model-architecture.md` "Cheap checks"). The flags are mirrored in `mge_detect.py` rather than imported so `--help` stays torch-free.

`--num-tiles` (default 1 = original frozen v1 behavior, unchanged) controls the multi-tile boundary-refinement mode — **now marked experimental/not-recommended** in its own `--help` text after 3 failed fix attempts; see `docs/boundary-refinement.md` ("Status: paused") before using it for anything but further investigation.

Console script: `[project.scripts] mge-detect = "mge_detect.mge_detect:main"` in `pyproject.toml` (added 2026-09-08, replacing a fragile `python3 src/mge_detect/mge_detect.py` pixi task that broke package-relative imports when run as a bare script).

## Package structure

`src/mge_detect/models/prophage/src/` is a real importable subpackage of `mge_detect` (`__init__.py` added at `models/`, `models/prophage/`, and `models/prophage/src/`; internal imports changed from `from src.X import Y` to relative `from .X import Y`).

**HPCC copy reconciled (2026-09-08)**: originally the HPCC's copy used the old `from src.X import Y` style while this local copy used relative imports (`from .X import Y`) — a deliberate divergence. As of the model-reengineering work (`docs/model-reengineering.md`), this local (relative-import) version of `model_rc_cnn.py`/`train_prophage.py` was pushed *to* the HPCC, overwriting the old-style copy (timestamped backups made first). Confirmed relative imports work fine there too under `python -m src.train_prophage --config ...` (namespace packages support this without `__init__.py`). The two copies should now match — if they're ever found to differ, treat the HPCC's as authoritative only if you can confirm which one actually trained the checkpoint you're using.

## Dependency fix (2026-09-08)

`pyproject.toml` previously declared `tensorflow` (nothing in the codebase imports it) and was missing `torch`, which every model file actually needs. Fixed: now declares `pytorch`, `pandas`, `pyfaidx`, `tqdm`, `pyyaml` (plus the pre-existing `biopython`, `scikit-learn`, `pytest`, `sphinx`, `mypy`).

**Fixed 2026-09-15**: `metadata/conda/mge_detect_05_11_2025.yaml` (a fully build-string-pinned, linux-64-only export of the TF-era `mge_explorer` HPCC env, still listing `tensorflow`/`keras` plus unrelated tooling like geNomad/xgboost/mmseqs2) was replaced by `metadata/conda/environment.yaml` — a clean, loosely-pinned spec mirroring `[tool.pixi.dependencies]` (torch, pandas, pyfaidx, tqdm, pyyaml, biopython, scikit-learn, plus dev tools) with an editable `pip: -e ../..` install. Keep it in sync with `pyproject.toml` when deps change. **Status (user, 2026-09-15): pixi is the supported dev/local setup; the conda yaml is an optional backup only, stored but not actively solve-tested. Eventual release distribution will be via bioconda, so neither file is the release artifact — a bioconda recipe will be.** Nothing referenced the old file; the HPCC SLURM launchers activate a separate `mge_explorer_gpu` env that neither file ever described.

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

**Follow-on work**: this experiment's ambiguity (real gain, but confounded with `pos_weight`/`amp` changes) motivated an isolated pooling-only rerun, now training on the HPCC alongside a new boundary-regression head — see `docs/model-reengineering.md`.

## Defaults changed + seed ensembling added (2026-09-21)

**`predict` defaults are now the measured operating point** (`docs/model-quality-baseline.md` "Operating point"), not the frozen v1 behavior:

| flag | was | now |
|---|---|---|
| `--region-builder` | `strict` | `coverage` |
| `--coverage-frac` | 0.3 | 0.1 |
| `--max-gap-bp` | 0 (off) | 50000 |
| `--bridge-prob` | 0.0 | 0.1 |

Rationale (user, 2026-09-21): both builders stay first-class and are measured separately, but the coverage+merge path carries more weight because **`strict` largely reproduces what geNomad already knows**, while the discovery path is the product. `strict` remains fully supported for the conservative, database-agreement comparison and is verified byte-identical to the pre-change frozen path.

**`--checkpoint` now accepts multiple paths**, with `--ensemble {max,mean}` (default `max`) combining their per-window score tracks. Passing the three `prophage_v2_pool` replicates is the recommended inference path: it removes the dependence of recall on which replicate happened to be packaged, and lands at the top of the seed range rather than a random draw from it (.978 covered recall in Default vs .961-.972 single-seed; .994 vs .989 in Exploratory). Cost is 3x inference time and ~30% more calls.

`max` rather than `mean` is deliberate and measured: mean-ensembling dilutes borderline detections below the seed threshold (.950/.978 vs .978/.994), the same failure that sank multi-tile averaging (`docs/boundary-refinement.md`).

Verified: a CLI run with three checkpoints and no other flags reproduces the sweep harness's `ensmax` variant coordinates exactly; single-checkpoint runs still record `ensemble: single`; `--region-builder strict` is unchanged.
