# Second locked test set (lockbox 2)

Written 2026-10-21, **before any lockbox-2 genome was selected, downloaded, labelled or scored.** This doc is the preregistration, written the same way as `docs/lockbox.md`: the selection rule, what is frozen before scoring, the endpoints and the reporting rules are fixed here first. Any change after the list exists is recorded below as an amendment, with its date and reason, and the original text stays.

## Why

The first lockbox was scored once, on 2026-10-19, under the defaults frozen at `ad5b28b`. On 2026-10-21 the user decided that the first release ships with **base-exact boundaries where an att site is found, the refinement pass on by default, and an edge range on every call** (`docs/att-boundaries.md`). That changes call coordinates, and covered50 depends on them. Under `docs/lockbox.md` ("If a default changes after scoring … a second lockbox would be needed"), the first lockbox's figure no longer describes the shipped calls. It stays the record of the earlier defaults and is not re-scored.

The new stages were tuned on the 42-genome sample and checked on the 18-genome sample and the literature pilot, so none of those is held out for them. This lockbox is.

## Selection (fixed)

`scripts/lockbox_select.py select`, run on ARC by `scripts/lockbox2_build.slurm`:

- **Source:** the first lockbox's NCBI summary, `analysis/lockbox/assembly_summary_refseq.20260929.txt`, refused unless its SHA-256 is `3459db2b5ab91279a26cc109b2917abe4a11e44105f5bf2d444d40011651e931`. The same file means the same eligibility rule over the same assemblies.
- **Eligible:** bacterial RefSeq `Complete Genome` assemblies, `latest` version, whose genus is absent from the corpus (train, val and test) **and** from `analysis/lockbox2/exclude_genera.txt`: every genus of the first lockbox (150), of the LAMBDA clean subset (15) and of the G5b Philympics set (6), 171 in all. Those three sets informed or measured earlier decisions. "uncultured" and similar genus names are excluded, as before.
- **One per genus:** reference genome, then representative genome, then any, ties broken by accession.
- **Sample:** 150 genera drawn with `numpy.random.default_rng(20261021)`, or all of them if fewer are eligible. 338 genera were eligible for the first lockbox and 150 of them were drawn, so at most 188 remain before the LAMBDA and G5b exclusions. The list's SHA-256 is recorded in `analysis/lockbox2/lockbox_selection.json`.
- **Truth:** geNomad 1.9.4 `annotate` + `find-proviruses`, DB 1.9, the corpus's settings. Truth is geNomad agreement, as everywhere else.

**Sample size.** The first lockbox had 202 elements in 150 genomes. A similar yield gives a 95% Wilson interval of about ±3 points at 96%. Strata under 10 kb will be thin again and are reported with their *n*, not as claims.

**What the agent may see before scoring:** the genome list, the selection table, and element counts per stratum. No model output on these genomes, as for the first lockbox.

## Frozen before scoring

Scoring happens **once**, after all of the following are final and committed. The commit is recorded here, and `scripts/lockbox2_run.slurm` (written at the freeze) refuses to run any other code.

- the checkpoints: `prophage_v2_pool_r1..r3`, the SHA-256s in `docs/research-integrity.md` C8 (bundled as `model/manifest.json`);
- the `predict` defaults: the three modes as in `docs/lockbox.md`, coverage 0.1, merge 50 kb / 0.1, `impute`, `scan`, edge snap 0.5; **refinement on** (5 kb windows, 250 bp step, threshold 0.6, α 0.6, ramp 15,000–60,000, split 4,000 / 0.05); **att search on** with the rule in `src/att.py` (search 25 / 10 kb, accept 10 / 5 kb, core ≤300 bp, ≤2 copies, entropy ≥1.5 bits, E ≤ 1e-6); edge ranges −0.055 / +0.494 × call length;
- the MD5 of every module the CLI loads, including the new `src/att.py`;
- the scorer, `scripts/lockbox_score.py --boundaries`, and this doc's endpoints.

## Endpoints (fixed)

- **Primary:** ensemble, Default mode, covered50 over all lockbox-2 elements, with a 95% Wilson interval and a genome-level bootstrap interval. Same definition as the first lockbox.
- **Secondary:**
  - the same for Exploratory and Conservative;
  - calls per genome per mode; covered recall (≥1 bp); conflated calls and absorbed elements; covered50 per size stratum, with *n*;
  - the gap to the 42- and 18-genome estimates (their covered50 is unchanged by the new stages: 0.994 / 0.972 / 0.944 and 1.000 / 1.000 / 0.962);
  - **boundary endpoints** (per mode): median IoU per element; the share of calls with att boundaries; the coverage of model-edge ranges (the geNomad edge inside `start_lo..start_hi` or `end_lo..end_hi`, on calls overlapping one element, contig-end edges excluded); att edges' distance to geNomad's (median, ≤100 bp, ≤1 kb); bp called and bp off every element.
- **Reported whatever the outcome.** No genome or element is dropped after scoring. A genome that fails to download or label is excluded *before* scoring, and listed.

### Reporting rules (fixed)

- **Recall:** the release states "≥95%" only if the primary endpoint's Wilson lower bound is at least 0.95. Otherwise it reports the number and its interval, as it does for the first lockbox.
- **Edge ranges:** they are called 90% ranges only if the measured Default coverage's Wilson interval includes 0.90 or lies above it. Otherwise the release states the measured coverage instead.
- **Att edges** are compared with geNomad's, whose borders sit at gene edges, so this distance measures agreement, not att accuracy. Base accuracy is stated only from author-curated truth (the literature pilot, `docs/att-boundaries.md`), as confirmation after use.

## Status

- [x] 2026-10-21: preregistration written. `lockbox_select.py` gained `--exclude_genera` and `--summary_sha256` (checked on a synthetic summary), `lockbox2_build.slurm` written, `exclude_genera.txt` built. **Not yet submitted.**
- [x] Build on ARC (`sbatch scripts/lockbox2_build.slurm`; the user submits). **Submitted 2026-10-23 as `848536`; finished the same day (exit 0, empty `.err`, no `failures.txt`); pulled and read 2026-10-26.** 188 genera were eligible (as predicted, at most 188), 150 drawn, 150 of 150 downloaded and labelled by geNomad. List SHA-256 `c6f57869ceaa8ce08a6b4750859a8cf683dd15195de4a23f6d31d0398974a4d0`; exclusion list SHA-256 `87f47161…`; summary SHA-256 matches the preregistered one. Checked locally: 150 distinct genera, none in `exclude_genera.txt` or the first lockbox, the list equals the selection table, 87 reference genomes and 63 others. **Truth: 186 elements in 73 genomes** (77 genomes carry none): 5–8 kb 2, 8–10 kb 2, 10–15 kb 26, 15–30 kb 39, 30–60 kb 94, 60 kb+ 23. That is 16 fewer than the first lockbox's 202, and the sub-10 kb strata (4 elements) are reported with *n* only. At *n* = 186 the "≥95%" rule (Wilson lower bound ≥ 0.95) needs at most 3 misses (183/186, lower bound 0.954); 4 misses gives 0.946. No model output exists on these genomes.
- [ ] Freeze the defaults (commit), write `scripts/lockbox2_run.slurm` with that commit's MD5s, run once on ARC, score once locally. **The last CPU gates passed on 2026-10-27** (steps 1–3). **Frozen at `0d71e2f` (2026-10-27; steps 4–6 done). The one run was `851599` (2026-10-27, 78 min). **Scored once on 2026-10-27: Default covered50 0.989 (0.962–0.997), so ≥95% is established held out; edge ranges cover 0.78, so they are not called 90% ranges ("Result").**
  - *Freeze (2026-10-27).* `0d71e2f` = `382efba`'s code plus docs and three LAMBDA result files, with no code change. Its archive (`arc_data/prophound_0d71e2f.tar.gz`, 27.4 MB, `src` and `examples`) has the same `src/` as the `wt_final` snapshot the matrix ran, except two scripts written later and never imported by `predict` (`lockbox2_run.slurm`, `split_support_check.py`) and `summary.csv`'s line endings. `analysis/lockbox2/frozen_md5.txt` lists all 745 files under `src/` (sha256 `21e7485299d65f75a47573c9e92a53c092f86756d0fb2703ae95c8e2396d7baf`), and `md5sum -c` passes in the extracted tree. The launcher's preflight was run locally on the extracted archive: it imports from it, and the example genome gives its 2 expected calls (MPS).
  - *Result (2026-10-27).* Head against new: 27 identical, and the 3 that differ are exactly `complete_exploratory`, `gappy_exploratory` and `shred_exploratory`. New against support-fix: **30 identical**. No reference case has a split piece whose support changes, so the matrix shows that the fix changes nothing else; its evidence is `tests/test_calls.py` and `scripts/split_support_check.py`. Both pass the rule in step 3.
  1. *Done.* **The paused matrices.** Two CPU reference matrices were frozen in place with `SIGSTOP` at 16:00 (the bracketed first letter in the patterns keeps `pkill -f` from matching the shell that runs it) (if the session or the machine restarted since, they are gone; see 2). Resume them with `pkill -CONT -f "[f]rom prophound.cli import main"; pkill -CONT -f "[r]eference_outputs.py run"`. They run from snapshots, so later edits cannot leak in: `arc_data/head_ef4faac` (the baseline), `arc_data/wt_expl_comp` (Exploratory composition only) and `arc_data/wt_final` (`382efba`'s code: composition plus the support fix; only `calls.py` differs from `wt_expl_comp`). Outputs go to `src/prophound/models/prophage/analysis/refout/{expl_comp_head_cpu,expl_comp_new_cpu,support_fix_cpu}`; logs to `arc_data/{expl_comp_head,expl_comp_new,support_fix}.log`. One chain runs head and then new; the support-fix matrix runs beside it.
  2. **If they died** (no `manifest.json` in an output directory, and no process left): `reference_outputs.py run` cannot resume. Delete that output directory and rerun it, with the same command as the next item but `run --code <snapshot> --inputs src/prophound/models/prophage/analysis/refout/inputs --out <dir>`, from the repo root (`pixi run python src/prophound/models/prophage/scripts/reference_outputs.py ...`).
  3. **Compare.** `reference_outputs.py compare expl_comp_head_cpu expl_comp_new_cpu`: only the Exploratory cases may differ (`complete_exploratory`, `gappy_exploratory`, `shred_exploratory`, and any other case run in Exploratory). `compare expl_comp_new_cpu support_fix_cpu`: only `seed_support` may differ (its columns in `regions*.tsv` and the GFF3, and `seed_support.calls_by_support` in `run_metadata.json`), and only in cases with split or refined calls. Anything else stops the freeze. Record the result in `docs/release-checklist.md` §2 (Exploratory composition, item 1) and in `docs/att-boundaries.md` "Split-piece support".
  4. **Freeze commit** (the user commits): `382efba` plus whatever the comparisons needed; no code change after it. Record its hash here.
  5. **Fill the launcher**: `git archive --format=tar.gz -o arc_data/prophound_<sha>.tar.gz <sha> src examples`; extract it to a scratch directory and write `analysis/lockbox2/frozen_md5.txt` as `md5sum` lines (`<md5>  src/...`) for every file under `src/` (`md5sum -c` on ARC reads it from the extracted root); set `FREEZE=<sha>` and `MD5S_SHA=<sha256 of frozen_md5.txt>` in `scripts/lockbox2_run.slurm`; commit both (second commit; the code stays the freeze commit's).
  6. **Push** (the user, `docs/environment.md`: `scp`, then a plain `ssh`, then verify with `ls -la`): the archive into `/work/pke835/mge_explorer/`, extracted to `/work/pke835/mge_explorer/prophound_<sha>/`; `frozen_md5.txt` into `data_refactored/analysis/lockbox2/`; the launcher into `data_refactored/scripts/`. Then `sbatch scripts/lockbox2_run.slurm` from `data_refactored`. The preflight runs the example genome on the V100 first: it also catches a missing dependency in `mge_explorer_gpu` (the package needs `biopython>=1.86`, never needed by an ARC job before).
  7. *Done 2026-10-27.* Pull `analysis/lockbox2/prophound/` and the logs, then score once: `lockbox_score.py --boundaries` (endpoints above; "≥95%" needs at most 3 misses of 186).

## Result (scored once, 2026-10-27)

ARC `851599` ran the archive of the freeze on a V100 in 78 minutes (09:34–10:52 on 2026-10-27). Scored once locally with `scripts/lockbox_score.py --boundaries`, giving `analysis/lockbox2/results.{json,elements.csv}`.

**Pre-scoring checks, all passed.**
- The preflight passed on ARC: the genome list's SHA-256, `frozen_md5.txt` (745 files) and its SHA-256, no extra modules, `prophound` importing from the archive, all 150 FASTAs, and the example genome's 2 expected calls on CUDA.
- All 450 runs (150 genomes × 3 modes) completed. Their directories match the genome list exactly, there is no `failures.txt`, and the job exited 0.
- `.err` holds only the preflight's own example-genome summary, which goes to stderr by design.
- Every run used CUDA and `boundary_refinement` = `fine_pass_then_att_then_composition`, with one checkpoint-hash set that equals `model/manifest.json`.
- `prophound_version` reads "unknown (not installed)" in every run's metadata, because the archive ran from `PYTHONPATH`. The code's identity is the MD5 preflight.

The truth is 186 geNomad elements in 73 of the 150 genomes; the other 77 have none.

**Primary endpoint: Default covered50 0.989 (184 of 186; Wilson 0.962–0.997; genome bootstrap 0.972–1.000), at 12.4 calls per genome.** The Wilson lower bound is above 0.95, so **by the preregistered rule the release may state "≥95%"** for Default, held out. The first lockbox (0.960, lower bound 0.924, under the pre-boundary defaults) could not.

| | Exploratory | **Default** | Conservative |
|---|---|---|---|
| covered50 | 0.984 (183) | **0.989 (184)** | 0.903 (168) |
| Wilson 95% | 0.954–0.994 | **0.962–0.997** | 0.852–0.938 |
| genome bootstrap 95% | 0.956–1.000 | 0.972–1.000 | 0.858–0.949 |
| covered (≥1 bp) | 1.000 | 0.989 | 0.925 |
| calls per genome | 37.4 | 12.4 | 5.0 |
| … on the 77 genomes with no element (descriptive) | 31.6 | 9.2 | 2.7 |
| conflated calls / absorbed elements | 5 / 10 | 5 / 10 | 5 / 10 |
| gap: the 42's covered50 minus this | +0.010 | −0.017 | +0.041 |
| gap: the 18's covered50 minus this | +0.016 | +0.011 | +0.059 |

**By size, Default** (Wilson 95%):

| stratum | n | covered50 | Wilson 95% |
|---|---|---|---|
| 5–8 kb | 2 | 1 of 2 | — |
| 8–10 kb | 2 | 2 of 2 | — |
| 10–15 kb | 26 | 1.000 | 0.871–1.000 |
| 15–30 kb | 39 | 1.000 | 0.910–1.000 |
| 30–60 kb | 94 | 0.989 | 0.942–0.998 |
| 60 kb+ | 23 | 1.000 | 0.857–1.000 |

Below 10 kb n = 2 per stratum: reported as counts, not claims. Conservative's weakest stratum is 15–30 kb (0.821, 0.673–0.910).

**Boundary endpoints:**

| | Exploratory | Default | Conservative |
|---|---|---|---|
| median IoU per element | 0.752 | 0.742 | 0.746 |
| calls with att boundaries | 0.055 (306 of 5,611) | 0.084 (157 of 1,863) | 0.108 (81 of 749) |
| model-edge range coverage (Wilson 95%) | 0.765 (238/311; 0.715–0.809) | **0.782 (237/303; 0.732–0.825)** | 0.777 (216/278; 0.724–0.822) |
| att edges vs geNomad's: median, ≤100 bp, ≤1 kb | 1,087 bp, 0.06, 0.48 (n = 48) | 871 bp, 0.10, 0.52 (n = 50) | 1,087 bp, 0.10, 0.48 (n = 48) |
| Mb called / Mb off every element | 82.8 / 76.4 | 45.0 / 38.7 | 21.3 / 15.5 |

- **Edge ranges: by the preregistered rule they are not "90% ranges".** Default coverage is 0.782, and its Wilson interval (0.732–0.825) lies below 0.90. The release states the measured coverage, about 78% held out, instead. On the development samples it was 0.88–0.92 (`docs/att-boundaries.md`), so the ranges were fitted too narrow for new genera.
- **Att edges sit further from geNomad's here** (median 0.9–1.1 kb, against 0.5–0.7 kb on the development samples). As preregistered, this measures agreement with geNomad's gene-edge borders, not att accuracy.
- **Median IoU is lower than on the tuning sample** (0.742 against 0.812), but in line with the check samples and LAMBDA (0.72–0.82).

**Misses (diagnosis after scoring; nothing is re-scored).**
- **Default misses 2 elements, both of which Exploratory reaches:**
  - a 6.1 kb element in `GCF_052828705.1` (no Default call; Exploratory 0.85);
  - a 34.9 kb element in `GCF_060445945.1` (no Default call; Exploratory 0.32).
- **Exploratory misses 3,** two of them partial detections that Default covers just over half:
  - `GCF_005601135.1` 597,338–618,137: Default 0.55, Exploratory 0.45;
  - `GCF_060445945.1` 5,499,362–5,535,423: Default 0.53, Exploratory 0.40.

  In both, Exploratory's looser composition gate (0.005) moved the right edge inward, by 2.3 kb and 5.9 kb. Default's gate (0.01) left that edge on the fine pass. This is the cost measured when Exploratory composition was adopted (one tuning-sample element, 0.994 → 0.989; `docs/att-boundaries.md`). Held out it is two elements, so Exploratory covers one element fewer than Default. Without them it would read 0.995. That is not a figure the release may use; it is the evidence for revisiting the gate after 0.1.0.

## Amendments

### Amendment 2026-10-21: composition edges added to the frozen defaults (before any selection)

Written the same day, before the genome list existed. After this doc was written, the user decided to ship composition edges in Default and Conservative (`docs/att-boundaries.md` "Composition edges"). The frozen defaults therefore also include `src/composition.py` (divergence ≥ 0.01, inward cap 0.15 × call length) with `--composition-edges auto`, and the edge ranges are 0.143 / 0.463 × call length in Default and Conservative, 0.055 / 0.494 in Exploratory. The scorer's run check expects `boundary_refinement` = `fine_pass_then_att_then_composition` in Default and Conservative and `fine_pass_then_att` in Exploratory. Endpoints and reporting rules are unchanged; the edge-range rule applies to each mode's own ranges.

### Amendment 2026-10-21: composition inward ceiling (before any selection)

Written the same day, still before the genome list existed. The composition cap is min(0.15 × call length, 8 kb) (`composition.MAX_INWARD_BP`), chosen after the LAMBDA re-score showed a 22.8 kb move cutting off a curated element (`docs/att-boundaries.md` "LAMBDA under the 0.1.0 pipeline"). This value is part of the frozen defaults.

### Amendment 2026-10-26: Exploratory composition edges and split-piece support (after selection, before any model output)

Written after the genome list and truth existed (build `848536`, Status above) and before any prophound output on a lockbox-2 genome. Only the selection table, the genome list and per-stratum element counts have been seen. Selection, truth, endpoints and reporting rules are unchanged.

- **Exploratory gains composition edges** (user decision 2026-10-23, `docs/att-boundaries.md` "Composition everywhere"): `composition.MODES` is all three modes, with the profile-divergence gate 0.005 in Exploratory (`MIN_JSD_EXPLORATORY`) and 0.01 in Default and Conservative. Edge ranges are 0.143 / 0.463 × call length in every mode. The scorer's run check now expects `boundary_refinement` = `fine_pass_then_att_then_composition` in all three modes (`lockbox_score.py` `BOUNDARY_STAGES`). Decided on the tuning, check, literature and LAMBDA samples, none of which shares a genus with this lockbox.
- **`seed_support` is counted on the extent the vote saw** (user decision 2026-10-26, `docs/att-boundaries.md` "Split-piece support"). A reported column only: calls and coordinates are unchanged (12 of 12 sample × mode rebuilds identical). No endpoint uses `seed_support`.
- **Execution.** `scripts/lockbox2_run.slurm` runs `prophound predict` from a `git archive` of the freeze commit, with the bundled model, instead of the ARC tree's `src/` as the first lockbox did: that tree predates the package layout. Its preflight checks every file of the archive against `analysis/lockbox2/frozen_md5.txt` (whose SHA-256 is in the launcher), the genome-list SHA-256, every FASTA, and that the example genome gives its expected calls on the V100 before any lockbox genome runs.

### Amendment 2026-10-27: commit hashes after a history rewrite (no code change)

History rewrite, before the first push. To keep three internal files off GitHub, they were stripped from the unpushed commits with `git filter-repo`; the 32 commits already on GitHub are unchanged. Every rewritten commit got a new hash, and another when the commits were then re-signed, so the hashes cited above refer to the pre-rewrite history (kept in the local backup bundle). The map gives the final, signed hashes. Code is unchanged: `git diff` of `src/` and `examples/` between each old and new commit is empty, so `frozen_md5.txt` and every MD5 list verify both. SLURM launchers and ARC paths (for example `prophound_0d71e2f`) keep the old names, because they record what ran. References in this doc to `docs/release-checklist.md`, `SESSIONS.md` and `CLAUDE.md` point to internal files that are not published; the rules they stated are in `docs/project-principles.md`. Map for the hashes this doc cites:

| cited | now |
|---|---|
| `0d71e2f` | `24f6c44` |
| `382efba` | `26880c3` |
| `ad5b28b` | `09c7756` |

The freeze is `24f6c44` (was `0d71e2f`). Run `851599` used the archive of `0d71e2f`, whose `src/` is identical.

### Amendment 2026-10-27: help text and docstrings updated after scoring (no output change); the lockbox figure stands

After the result, two hash-checked modules changed in text only:
- `src/options.py`: the `--mode` help's call counts now read about 37 / 12 / 5 per genome, this lockbox's figures (they were about 32 / 11 / 5, the first lockbox's).
- `src/att.py`: two docstrings no longer call the edge ranges "90%", since they held 0.78 here.

Checks:
- Removing docstrings leaves `att.py`'s syntax tree identical to the scored code's, and `options.py`'s differs only in that help string.
- The fast tests (198) and the smoke test pass.
- CPU reference matrix against `support_fix_cpu`, which ran the scored code: **30 of 30 identical** (`analysis/refout/helptext_cpu`, 2026-10-27 12:52–14:52).

Nothing that produces a call changed, so the result above describes the shipped code.
