# CLI & Packaging

Living doc — update in place, don't append a history here (that's `SESSIONS.md`'s job).

## Current CLI shape

`prophound predict` (entrypoint `src/prophound/cli.py`, dispatching into `src/prophound/models/prophage/src/predict_genome.py`). Every option, default and mode table is defined once, in the torch-free `src/prophound/models/prophage/src/options.py` (`add_predict_arguments`). `prophound predict` and `python -m src.predict_genome` both use it; the latter still accepts `--input_fasta` for the HPCC launchers. `prophound predict --help` lists the full set in basic, advanced and experimental sections. `-q` prints only warnings, errors and the summary; `-v` adds stage timings and full tracebacks. A failure the user can fix ends in one `error:` line (exit 1; 130 on Ctrl-C). `--input` takes FASTA, GenBank/RefSeq, EMBL or GFF3 with `##FASTA`, optionally gzip/bzip2/xz compressed; the format is read from the content (`src/inputs.py`), and reads, alignments and protein are refused with a reason. Several genomes run as a **batch** (repeat `--input`, give a directory, or `--input-list`): one subdirectory per genome plus combined calls and `batch_summary.tsv`, checkpoints loaded once (`docs/output-files.md` "Batches"). Everything after scoring (building calls, the seed vote, edge snap, the N filter, seed support) is in the torch-free `src/calls.py`, which `scripts/mode_recheck.shipped_calls` also runs, so offline scorers rebuild exactly the CLI's calls (with `genomes_dir`; without it, snap and the N filter are skipped with a warning). Calls removed by the vote or the N filter are kept aside with their reason, not deleted (`docs/tech-debt-audit.md` Phase 2b). The core:

```
prophound predict --input <fasta> --outdir <dir> --config <config.yaml> --checkpoint <best.pt> [<best.pt> ...] \
  [--mode exploratory|default|conservative] [--threshold F] [--seed-threshold F] [--extend-threshold F] \
  [--device auto|cpu|cuda|rocm|mps|xpu|intel] [--batch-size N] [--num-workers N] \
  [--coverage-frac F] [--max-gap-bp N] [--bridge-prob F] [--write-region-fasta] ...
```

Outputs: `regions.tsv`, `regions.filtered.tsv`, `regions.bed`, `regions.gff3`, `windows.tsv`, `run_metadata.json`, plus `regions.fna` (`--write-region-fasta`) and `fine_windows.tsv` (refinement on). **Every column is described in `docs/output-files.md`**, the output contract. The frozen `notes/V1_INFERENCE_CONTRACT.md` and `notes/V1_OUTPUT_FILES.md` predate it.

**Removed 2026-10-20** (`docs/tech-debt-audit.md`, decision 4): `--region-builder` (the `strict` builder; `coverage` is now the only one, and `run_metadata.json` still records `region_builder: coverage_fraction`), `--num-tiles` with the multi-tile builder and `bins.tsv`, and the boundary-regression head. Default outputs are byte-identical. The paragraphs below on those options are history. The strict builder survives as a private copy in `scripts/coverage_builder_checks.py`, the one script that compares against it.

**The model ships with the package (2026-10-20).** With no `--checkpoint`, `predict` uses the bundled model in `src/prophound/models/prophage/model/`: three weights-only checkpoints (5.5 MB each; the training files were 16.5 MB with optimizer and RNG state), `model.yaml`, and `manifest.json` with SHA-256s that are checked at every run. So `prophound predict --input genome.fna --outdir out` is the whole command. `--checkpoint`/`--config` use another model. Checkpoints load with `torch.load(weights_only=True)`. A full training checkpoint is refused, with guidance, unless `--trust-checkpoint` is given, because it loads by unpickling arbitrary objects. `scripts/strip_checkpoint.py` makes weights-only files. The wheel is 15.5 MB with the model. `reference_outputs.py` and `benchmark.py` use the bundled files by default; `--model training` gives the original checkpoints for code from before Phase 2. Old ARC launchers that pass training checkpoints now need `--trust-checkpoint`. Since 2026-10-20 the config needs only `data.window_len`, `data.stride` and the `model` section for inference. The training paths (`data.manifest_csv`, `data.splits_dir`) are optional, and unknown keys are an error (`src/config.py`). `--config`/`--checkpoint` are required, with no bundled default weights — model-weight bundling is an open packaging question upstream (`notes/V1_PACKAGE_BOUNDARY.md`) that hasn't been decided, so this doesn't guess at an answer.

`--region-builder` (added 2026-09-17; default `strict` = the frozen v1 builder, byte-identical output to before the flag existed — verified) selects how window scores become regions. `coverage` uses `predict_genome.build_region_table_coverage`: a stride bin is inside if ≥ `--coverage-frac` (default 0.3) of its covering windows individually score ≥ the extend threshold, instead of the strict builder's "mean of all covering windows ≥ extend". On the 18-genome eval set that lifts region recall 76.9%→88.5% (frac 0.3) / 94.2% (frac ≤0.1) and median IoU 0.56→0.75, at 2–4× the calls per genome; `run_metadata.json` records `region_builder: coverage_fraction` and `coverage_frac`. Not combinable with `--num-tiles > 1` (hard error). **Still opt-in, not the shipped default**: frac was chosen on the same 52 regions it was scored on — pending a held-out, size-stratified check (`docs/multi-model-architecture.md` "Cheap checks"). The flags are mirrored in `cli.py` rather than imported so `--help` stays torch-free.

`--num-tiles` (default 1 = original frozen v1 behavior, unchanged) controls the multi-tile boundary-refinement mode — **now marked experimental/not-recommended** in its own `--help` text after 3 failed fix attempts; see `docs/boundary-refinement.md` ("Status: paused") before using it for anything but further investigation.

Console script: `[project.scripts] prophound = "prophound.cli:main"` in `pyproject.toml` (added 2026-09-08, replacing a fragile `python3 src/prophound/cli.py` pixi task that broke package-relative imports when run as a bare script).

## Package structure

`src/prophound/models/prophage/src/` is a real importable subpackage of `prophound` (`__init__.py` added at `models/`, `models/prophage/`, and `models/prophage/src/`; internal imports changed from `from src.X import Y` to relative `from .X import Y`).

**HPCC copy reconciled (2026-09-08)**: originally the HPCC's copy used the old `from src.X import Y` style while this local copy used relative imports (`from .X import Y`) — a deliberate divergence. As of the model-reengineering work (`docs/model-reengineering.md`), this local (relative-import) version of `model_rc_cnn.py`/`train_prophage.py` was pushed *to* the HPCC, overwriting the old-style copy (timestamped backups made first). Confirmed relative imports work fine there too under `python -m src.train_prophage --config ...` (namespace packages support this without `__init__.py`). The two copies should now match — if they're ever found to differ, treat the HPCC's as authoritative only if you can confirm which one actually trained the checkpoint you're using.

## Dependency fix (2026-09-08)

`pyproject.toml` previously declared `tensorflow` (nothing in the codebase imports it) and was missing `torch`, which every model file actually needs. Fixed: now declares `pytorch`, `pandas`, `pyfaidx`, `tqdm`, `pyyaml` (plus the pre-existing `biopython`, `scikit-learn`, `pytest`, `sphinx`, `mypy`).

**Fixed 2026-09-15**: `metadata/conda/mge_detect_05_11_2025.yaml` (a fully build-string-pinned, linux-64-only export of the TF-era `mge_explorer` HPCC env, still listing `tensorflow`/`keras` plus unrelated tooling like geNomad/xgboost/mmseqs2) was replaced by `metadata/conda/environment.yaml` — a clean, loosely-pinned spec mirroring `[tool.pixi.dependencies]` (torch, pandas, pyfaidx, tqdm, pyyaml, biopython, scikit-learn, plus dev tools) with an editable `pip: -e ../..` install. Keep it in sync with `pyproject.toml` when deps change. **Status (user, 2026-09-15): pixi is the supported dev/local setup; the conda yaml is an optional backup only, stored but not actively solve-tested. Eventual release distribution will be via bioconda, so neither file is the release artifact — a bioconda recipe will be.** Nothing referenced the old file; the HPCC SLURM launchers activate a separate `mge_explorer_gpu` env that neither file ever described.

## Packaging fix (2026-10-19)

`docs/tech-debt-audit.md` P1, exclude-only (user decision).
- `[project] dependencies` now lists the runtime packages: torch, numpy, pandas, pyfaidx, tqdm, pyyaml and, until inference stops importing `train_prophage`, scikit-learn. Before this the wheel declared none.
- Hatch excludes the research tree from the wheel and the sdist. The wheel went from ~28 MB to 60 KB.
- The sdist holds only `src/`, `tests/`, README and licence files, so internal `docs/` stay out.
- `LICENSE` and `NOTICE` now ship: `license-files` had named a non-existent `LICENSE.txt`.
- `ruff` is a dev dependency with a `[tool.ruff]` lint target.

## Devices (2026-10-19)

- `--device auto` prefers CUDA, then MPS, then XPU, then CPU. On Apple Silicon it picks MPS.
- **MPS gave wrong scores until 2026-10-19:** adaptive max pooling (`docs/tech-debt-audit.md` S10). It is fixed in `model_rc_cnn.py`.
- Calls are the same on M4 CPU, M4 MPS and V100S CUDA. Scores agree to ~1e-5 but not byte for byte (`docs/benchmarks.md`).
- Speed on the M4: ~9 s/Mb on MPS (43.5 s per 5 Mb genome), ~9 min/Mb on CPU (45 min per 5 Mb genome) with the generic conda-forge torch build.

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

**Loose end**: `notes/V1_INFERENCE_CONTRACT.md`, `notes/V1_OUTPUT_FILES.md` and `notes/V1_OPERATING_MODES.md` describe the pre-2026-09-21 CLI — single `--checkpoint`, `strict` builder as default, no merging, no `ensemble` field in `run_metadata.json`. Those notes are a frozen historical snapshot and are deliberately not being edited (`CLAUDE.md`), but anyone reading them for the *current* contract will be misled. This doc is the current contract; if a release note or package description is ever generated from the `notes/` files, regenerate it from here instead.

## Input robustness (2026-09-22)

Measured under the standing RC/N check (`docs/scientific-validity.md`, "Input-robustness standing check"). Three findings that affect what the CLI does with real assemblies:

- **Contigs shorter than `window_len` (10 kb) are silently skipped** by `build_windows`, with no warning. A mixed assembly returns a normal-looking `regions.tsv` while whole contigs went unexamined; an all-short input fails loudly instead. Fix direction: report skipped contigs and skipped base count.
- **IUPAC ambiguity codes collapse window scores to zero** (two contiguous N bases take p=1.000 to 0.000). Affects draft assemblies, MAGs and environmental data. Two equivalent fixes are measured; neither is adopted yet because either changes every score.
- **Calls depend on window-grid phase** — shifting a contig's start by 777 bp turned 3 regions into 5 (marginal calls only; confident ones are stable).

Reverse-complement handling is *not* a concern: window scores and whole-contig region calls are both exactly invariant.

## Window size and fine refinement flags (added 2026-09-23, experimental, off by default)

- `--window-len` / `--stride`: override the checkpoint's trained 10 kb / 1 kb at inference. The pooling accepts any length, but the mode thresholds were calibrated at 10 kb, so a warning is printed whenever the length differs. `run_metadata.json` records both the used and the trained values.
- `--refine-window-len W` (0 = off), plus `--refine-stride`, `--refine-threshold`, `--refine-margin`, `--edge-alpha`, `--split-min-gap-bp` and `--split-max-prob`: the fine pass that splits conflated calls and places edges. It needs `--region-builder coverage`. It writes an extra `fine_windows.tsv`, and `regions.tsv` gains `refine_status`, `coarse_start`, `coarse_end` and `n_pieces`. Design: `docs/downstream-postprocessing.md` "One pipeline for all four".
- The flags are defined twice, in `predict_genome.add_multiscale_args` and in `cli.py`, because the top-level CLI avoids importing torch at parse time. Both parse identically (checked). Keep them in sync.

## Seed vote per mode (2026-09-29, user decision)

Each mode is now a seed threshold **plus** a call-level seed vote (`predict_genome.MODE_THRESHOLDS`, `MODE_MIN_SEED_SUPPORT`). An ensemble call is kept only if at least *k* checkpoints make an overlapping call of their own at the same thresholds; its extent still comes from the `max` track.

Sample roles (`docs/research-integrity.md` C1): the 42 is the **tuning** sample; the 18 is **confirmation after use**. Neither is held out.

| mode | seed | vote *k* | 42: covered50 / calls | 18: covered50 / calls | phage-free calls |
|---|---|---|---|---|---|
| exploratory | 0.90 | 1 (all kept, `seed_support` flagged) | 0.994 / 35.6 | 1.000 / 39.6 | 27.1 |
| default | 0.95 | 3 | 0.972 / 14.1 | 1.000 / 13.3 | 7.4 |
| conservative | 0.995 | 3 | 0.944 / 8.5 | 0.962 / 6.4 | 2.45 |

Short-contig thresholds (`SHORT_CONTIG_THRESHOLDS`) are 0.95 / 0.95 / 0.99, set by the "never exceed" rule (user). Each is the lowest threshold whose shredded false-call rate, with the mode's vote, stays at or below the intact 10 kb pipeline's (`822053`).

History, all on 2026-09-29:
- Before the vote, the modes were 0.95 / 0.995 / 0.999 with no vote.
- Default was first set to 0.98 with 3 of 3. That kept geNomad recall at 22–27% fewer calls, but lost small and geNomad-negative elements (spike-in 0.30 → 0.25, in situ 20 → 18 of 21, *B. subtilis* prophage 3). So it moved to 0.95 with 3 of 3, which beats the old Default on every measure at the same cost (`docs/error-analysis.md`).
- Exploratory then moved to 0.90 to keep a reach beyond Default. Recall saturates there.

- **`--min-seed-support K`** overrides the mode's *k*; `1` keeps every call. The pre-2026-09-29 outputs are reproduced with `--threshold 0.995 --min-seed-support 1` (Default) or `--threshold 0.999 --min-seed-support 1` (Conservative).
- **The vote needs several checkpoints** and the coverage builder with one tile. With fewer checkpoints than *k*, it votes at the number available, and with one checkpoint (or the strict builder) no vote applies. Either way a warning says so, because the thresholds were calibrated with the vote.
- **It is applied before refinement**, to coarse and short-contig calls separately, so the fine pass never works on a call the vote will drop. `run_metadata.json` `seed_support` records `min_seed_support`, `requested_min_seed_support` and `calls_dropped_by_vote`.
- **`--write-seed-scores`** (opt-in) adds `prob_ckpt1..N` to `windows.tsv`, so any vote can be re-applied offline. `scripts/mode_recheck.py` `shipped_calls()` does exactly what the CLI does, checked end to end: offline equals online on a planted fixture, and the offline rebuild reproduces the 42-genome and 18-genome frontier numbers to 3 dp.
- Per-seed calls are now built from the pre-builder window frame, so they no longer depend on the builder preserving row order (it did, but only by construction).
- **Short-contig thresholds**: re-derived under the vote (above). `final_modes_rescore.slurm` confirms them under the final modes.

## Usability pass (2026-09-29)

From the 2026-09-28 audit against Seemann's recommendations (`docs/claims-and-evidence.md`). New module: `src/prophound/models/prophage/src/io_checks.py`.

- **`prophound --version`**. The version comes from the installed distribution (`pyproject.toml`); run from an uninstalled source tree it reads "unknown (not installed)". `run_metadata.json` records `prophound_version` and the SHA-256 of every checkpoint used.
- **stdout is empty**. Every `[warn]`, `[done]` and `[summary]` line, and the `tqdm` progress bars, go to stderr.
- **Input validation** runs before anything touches the file (one streaming pass), with a one-line `error:` message and a nonzero exit. It rejects:
  - a missing or empty file;
  - gzip (pyfaidx would fail obscurely);
  - text that isn't FASTA, or a header with no ID;
  - records with no sequence;
  - duplicate contig IDs;
  - protein FASTA (below 90% nucleotide or IUPAC letters).

  The same pass counts ambiguity codes per code into `run_metadata.json` `input`.
- **Run summary** on stderr at the end of every run: input size and ambiguity load, short contigs scanned or never scanned, calls and the percentage of input called, calls dropped by the seed vote, and the line "calls are candidates for follow-up, not annotations".
- **Outputs**: `regions.bed` (BED6, 0-based, score = max_prob × 1000) and `regions.gff3` (1-based, SO type `prophage`, with `note=candidate`, `max_prob`, `mean_prob`, `seed_support`, `n_frac`, `contig_edge` and `boundary_kind`) are written every run. `--write-region-fasta` adds `regions.fna`. `run_metadata.json` `outputs` lists what was written.
- **`predict --help` is grouped**:
  - basic (7): input, outdir, mode, config, checkpoint, device, region FASTA;
  - advanced (13): thresholds, seed vote, ambiguity, merging, short contigs, batch size;
  - experimental (14): tiles, grid offset, window length and stride, refinement, the region builder and the ensemble rule.

  Only the layout changes; definitions and defaults don't, and `tests/test_cli_defaults.py` still checks the two parsers agree.
- **Edge snap (2026-10-05, on by default)**: `--edge-snap-frac` (default 0.5, so 5 kb at the 10 kb window; 0 = off). A call edge within that fraction of `--window-len` of a contig end is moved onto it. Short-contig calls use the same distance, so any call on a contig of 10 kb or less spans that contig. `regions.tsv` and the GFF3 gain `contig_edge` (none/start/end/both). It marks every edge on a contig end, snapped or not, meaning the reported boundary is the end of the sequence and the element may continue on another contig. Applied after the short-contig merge and before the N-gap filter. `run_metadata.json` `edge_snap`. The fix for fragmented input: `docs/claims-and-evidence.md` "Amendment 2026-10-05: freeze".
- **N-gap filter (2026-09-29)**: calls that are ≥50% N are dropped (`MAX_CALL_N_FRAC`). `n_frac` is in `regions.tsv` and the GFF3; the counts are in `run_metadata.json` `n_gap_filter` and in the summary (`docs/v1-scope.md` A.3).
- **ARC note**: `predict_genome.py` now imports `io_checks.py`. Any push of `src/predict_genome.py` must include `src/io_checks.py`.
- Tests: `tests/test_io_checks.py` covers every rejection and the BED/GFF3 coordinate conventions (31 tests in total). An end-to-end check with planted tracks confirmed empty stdout, BED, GFF3 and FASTA agreeing with `regions.tsv`, the metadata blocks, and a clear protein-input error.

## Decisions 2026-09-24 (user)

- **`--ambiguity-policy` defaults to `impute`** in both CLIs and in the library (`predict_genome.DEFAULT_AMBIGUITY_POLICY`). It is a no-op on ambiguity-free sequence: checked end to end byte-for-byte against the prior commit, and unit-tested in `tests/test_encoding.py` (816446 also showed covered recall identical to 3 dp on the eval sample). Numbers published before 2026-09-24 are reproduced with `--ambiguity-policy n-channel`. Note: `scripts/spikein_harness.py` keeps its own explicit `n-channel` default so its recorded curves stay reproducible.
- **Distribution**: the model will be bundled through **bioconda**; install and quickstart docs are written at the very end, not now.
- **Short contigs**: `--short-contigs {skip,scan}` built 2026-09-24. Seed support (2026-09-28): ensemble runs add a `seed_support` column to `regions.tsv` (checkpoints whose own calls overlap the call; `docs/error-analysis.md`), a ranking signal only — no call changes. Refinement (2026-09-28): `--refine-threshold` defaults to 0.95 (was the mode's extend threshold), and the split gates to 4000 / 0.05 (were 2000 / 0.1), so the opt-in tight preset runs as measured; `tests/test_cli_defaults.py` keeps the two CLIs' defaults in sync. **Default `scan` since 2026-09-28** (user decision; `skip` reproduces older outputs). Calibrated per mode: Exploratory 0.95, Default 0.99, Conservative 0.998, at both 5 kb and 3 kb (`SHORT_CONTIG_THRESHOLDS`). An explicit `--short-contig-threshold` overrides them and is recorded as `calibrated: false`. `run_metadata.json` `short_contigs.thresholds` is now per scale (it replaces `threshold_seed`/`threshold_extend`). A warning is always printed. Design and calibration: `docs/adaptive-multiscale-detection.md` Part 4.
- **Flagged for later evaluation** (not v1 blockers):
  1. Per-call confidence or ranking for wet-lab prioritisation (today only raw mean/max probability, which per-genome thresholds would make non-comparable across genomes).
  2. Recall claims must name their overlap rule. ≥1 bp credits grazes; see `docs/downstream-postprocessing.md` "Per-element dig" for covered50.
  3. The duplicated flag definitions in `cli.py` and `predict_genome.py` (kept so the top-level CLI parses without importing torch; checked equal, but a maintenance hazard).
- CLI surface as of today: one subcommand (`predict`), 29 options, most of them research knobs; `--help` does not yet separate user-facing from experimental options.

