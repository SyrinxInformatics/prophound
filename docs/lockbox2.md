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
- [ ] Build on ARC (`sbatch scripts/lockbox2_build.slurm`; the user submits).
- [ ] Freeze the defaults (commit), write `scripts/lockbox2_run.slurm` with that commit's MD5s, run once on ARC, score once locally.

## Amendments

### Amendment 2026-10-21: composition edges added to the frozen defaults (before any selection)

Written the same day, before the genome list existed. After this doc was written, the user decided to ship composition edges in Default and Conservative (`docs/att-boundaries.md` "Composition edges"). The frozen defaults therefore also include `src/composition.py` (divergence ≥ 0.01, inward cap 0.15 × call length) with `--composition-edges auto`, and the edge ranges are 0.143 / 0.463 × call length in Default and Conservative, 0.055 / 0.494 in Exploratory. The scorer's run check expects `boundary_refinement` = `fine_pass_then_att_then_composition` in Default and Conservative and `fine_pass_then_att` in Exploratory. Endpoints and reporting rules are unchanged; the edge-range rule applies to each mode's own ranges.
