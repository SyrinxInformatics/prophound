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
- [ ] Freeze the defaults (commit), write `scripts/lockbox2_run.slurm` with that commit's MD5s, run once on ARC, score once locally. **The last CPU gates passed on 2026-10-27** (steps 1–3; all code is in `382efba`; the launcher is drafted with two placeholders). **Next is step 4.**
  - *Result (2026-10-27).* Head against new: 27 identical, and the 3 that differ are exactly `complete_exploratory`, `gappy_exploratory` and `shred_exploratory`. New against support-fix: **30 identical**. No reference case has a split piece whose support changes, so the matrix shows that the fix changes nothing else; its evidence is `tests/test_calls.py` and `scripts/split_support_check.py`. Both pass the rule in step 3.
  1. *Done.* **The paused matrices.** Two CPU reference matrices were frozen in place with `SIGSTOP` at 16:00 (the bracketed first letter in the patterns keeps `pkill -f` from matching the shell that runs it) (if the session or the machine restarted since, they are gone; see 2). Resume them with `pkill -CONT -f "[f]rom prophound.cli import main"; pkill -CONT -f "[r]eference_outputs.py run"`. They run from snapshots, so later edits cannot leak in: `arc_data/head_ef4faac` (the baseline), `arc_data/wt_expl_comp` (Exploratory composition only) and `arc_data/wt_final` (`382efba`'s code: composition plus the support fix; only `calls.py` differs from `wt_expl_comp`). Outputs go to `src/prophound/models/prophage/analysis/refout/{expl_comp_head_cpu,expl_comp_new_cpu,support_fix_cpu}`; logs to `arc_data/{expl_comp_head,expl_comp_new,support_fix}.log`. One chain runs head and then new; the support-fix matrix runs beside it.
  2. **If they died** (no `manifest.json` in an output directory, and no process left): `reference_outputs.py run` cannot resume. Delete that output directory and rerun it, with the same command as the next item but `run --code <snapshot> --inputs src/prophound/models/prophage/analysis/refout/inputs --out <dir>`, from the repo root (`pixi run python src/prophound/models/prophage/scripts/reference_outputs.py ...`).
  3. **Compare.** `reference_outputs.py compare expl_comp_head_cpu expl_comp_new_cpu`: only the Exploratory cases may differ (`complete_exploratory`, `gappy_exploratory`, `shred_exploratory`, and any other case run in Exploratory). `compare expl_comp_new_cpu support_fix_cpu`: only `seed_support` may differ (its columns in `regions*.tsv` and the GFF3, and `seed_support.calls_by_support` in `run_metadata.json`), and only in cases with split or refined calls. Anything else stops the freeze. Record the result in `docs/release-checklist.md` §2 (Exploratory composition, item 1) and in `docs/att-boundaries.md` "Split-piece support".
  4. **Freeze commit** (the user commits): `382efba` plus whatever the comparisons needed; no code change after it. Record its hash here.
  5. **Fill the launcher**: `git archive --format=tar.gz -o arc_data/prophound_<sha>.tar.gz <sha> src examples`; extract it to a scratch directory and write `analysis/lockbox2/frozen_md5.txt` as `md5sum` lines (`<md5>  src/...`) for every file under `src/` (`md5sum -c` on ARC reads it from the extracted root); set `FREEZE=<sha>` and `MD5S_SHA=<sha256 of frozen_md5.txt>` in `scripts/lockbox2_run.slurm`; commit both (second commit; the code stays the freeze commit's).
  6. **Push** (the user, `docs/environment.md`: `scp`, then a plain `ssh`, then verify with `ls -la`): the archive into `/work/pke835/mge_explorer/`, extracted to `/work/pke835/mge_explorer/prophound_<sha>/`; `frozen_md5.txt` into `data_refactored/analysis/lockbox2/`; the launcher into `data_refactored/scripts/`. Then `sbatch scripts/lockbox2_run.slurm` from `data_refactored`. The preflight runs the example genome on the V100 first: it also catches a missing dependency in `mge_explorer_gpu` (the package needs `biopython>=1.86`, never needed by an ARC job before).
  7. Pull `analysis/lockbox2/prophound/` and the logs, then score once: `lockbox_score.py --boundaries` (endpoints above; "≥95%" needs at most 3 misses of 186).

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
