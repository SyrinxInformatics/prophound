# The held-out test

The test behind prophound's headline figures. 150 complete bacterial genomes from genera that no earlier step had seen were chosen, labeled and scored **once**, with the rules for selection, scoring and reporting fixed in advance. Those rules are kept below exactly as written.

## In short

- **Default mode covered at least half of 184 of 186 known prophages: 0.989** (95% Wilson interval 0.962–0.997), at 12.4 candidates per genome.
- The rule fixed in advance allowed the release to state "at least 95%" only if the interval's lower bound reached 0.95. It did, so **at least 95% of known elements are covered on held-out genomes.**
- "Known" means geNomad's calls. This is agreement with geNomad, not accuracy ([Claims and evidence](claims-and-evidence.md)).
- **Edge ranges held the true end in about 78% of cases,** below the 90% the rule required to call them 90% ranges, so they are reported as measured.

## What was tested

- **Genomes.** 150 complete RefSeq bacterial genomes, one per genus, drawn at random from the 188 genera absent from prophound's training data and from every earlier test (the first held-out test, the LAMBDA benchmark and the Philympics genomes). 87 are NCBI reference genomes.
- **Truth.** geNomad 1.9.4 with database 1.9, run the same way as for the training labels: **186 elements in 73 genomes.** The other 77 genomes carry none.
- **Code.** The 0.1.0 pipeline, frozen before the run: the bundled three-checkpoint model and every default, including the boundary stages. The run's launcher refused to start unless every file of the frozen code matched its recorded checksum.
- **Blinding.** Before scoring, only the genome list, the selection table and the element counts per size class had been seen. No model output on these genomes existed.
- **One run.** All 450 runs (150 genomes × 3 modes) ran once on a cluster V100S GPU, and the results were scored once ([`analysis/lockbox2/`](../prophound/models/prophage/analysis/lockbox2/): `results.json`, `results.elements.csv`; scorer `scripts/lockbox_score.py --boundaries`).

**Why a second test.** An earlier held-out test, on 150 other genomes, scored the defaults as they stood before the boundary stages (att-site ends, the refinement pass and composition edges) were added: Default 0.960 (0.924–0.980), whose lower bound fell short of 0.95. The boundary stages change where calls start and end, and so what they cover, which meant that test no longer described the calls that ship. It stays as the record of the earlier defaults, and this test was set up for the new ones.

## Results

| | Exploratory | **Default** | Conservative |
|---|---|---|---|
| known elements covered at least half (covered50) | 0.984 (183) | **0.989 (184)** | 0.903 (168) |
| Wilson 95% interval | 0.954–0.994 | **0.962–0.997** | 0.852–0.938 |
| genome-bootstrap 95% interval | 0.956–1.000 | 0.972–1.000 | 0.858–0.949 |
| known elements touched at all | 1.000 | 0.989 | 0.925 |
| candidates per genome | 37.4 | 12.4 | 5.0 |
| … on the 77 genomes with no known element | 31.6 | 9.2 | 2.7 |
| calls holding more than one element / elements in them | 5 / 10 | 5 / 10 | 5 / 10 |

**By element size, Default:**

| size | elements | covered50 | Wilson 95% |
|---|---|---|---|
| 5–8 kb | 2 | 1 of 2 | |
| 8–10 kb | 2 | 2 of 2 | |
| 10–15 kb | 26 | 1.000 | 0.871–1.000 |
| 15–30 kb | 39 | 1.000 | 0.910–1.000 |
| 30–60 kb | 94 | 0.989 | 0.942–0.998 |
| over 60 kb | 23 | 1.000 | 0.857–1.000 |

Below 10 kb there are two elements per class, reported as counts, not as claims.

**Boundaries:**

| | Exploratory | Default | Conservative |
|---|---|---|---|
| median overlap with the element (IoU) | 0.752 | 0.742 | 0.746 |
| calls with att-site ends | 0.055 (306 of 5,611) | 0.084 (157 of 1,863) | 0.108 (81 of 749) |
| estimated ends whose range holds geNomad's end | 0.765 (0.715–0.809) | **0.782 (0.732–0.825)** | 0.777 (0.724–0.822) |
| att ends against geNomad's ends: median distance; within 100 bp; within 1 kb | 1,087 bp; 0.06; 0.48 | 871 bp; 0.10; 0.52 | 1,087 bp; 0.10; 0.48 |
| Mb called / Mb called off every known element | 82.8 / 76.4 | 45.0 / 38.7 | 21.3 / 15.5 |

## What the results mean

- **The 95% target holds held out for Default,** by a rule fixed before the test.
- **The ranges were fitted too narrow for new genera.** They held 0.88–0.92 of true ends on the development samples and 0.78 here. Widening them is a change for after 0.1.0.
- **Att ends sit about 0.9–1.1 kb from geNomad's ends here** (0.5–0.7 kb on development samples). geNomad's ends sit on gene edges, so this measures agreement with geNomad, not the accuracy of the att ends. Against authors' coordinates in the literature, att ends sat a median 22 bp away ([release text](release-text-0.1.0.md)).
- **Overlap (IoU) is lower than on the tuning sample** (0.742 against 0.812) and in line with the other samples and LAMBDA (0.72–0.82).
- **The misses.** Default misses two elements (6.1 kb and 34.9 kb), both of which Exploratory reaches. Exploratory misses three, two of them partial detections that Default covers just over half: Exploratory's looser composition gate moved one edge of each inward (by 2.3 and 5.9 kb). That is why Exploratory covers one element fewer than Default. It is evidence for revisiting the gate after 0.1.0, not a figure the release may use.

---

## Preregistration (as written)

Everything below this line is the preregistration as written before the test, followed by every amendment made after it, each dated and unedited. References in it to other docs (`docs/lockbox.md`, `docs/att-boundaries.md` and others) and to cluster job numbers point to internal records that are not published; the figures they support are in the plain account above and in the [release text](release-text-0.1.0.md).

Written 2026-10-21, **before any lockbox-2 genome was selected, downloaded, labelled or scored.** This doc is the preregistration, written the same way as `docs/lockbox.md`: the selection rule, what is frozen before scoring, the endpoints and the reporting rules are fixed here first. Any change after the list exists is recorded below as an amendment, with its date and reason, and the original text stays.

### Why

The first lockbox was scored once, on 2026-10-19, under the defaults frozen at `ad5b28b`. On 2026-10-21 the user decided that the first release ships with **base-exact boundaries where an att site is found, the refinement pass on by default, and an edge range on every call** (`docs/att-boundaries.md`). That changes call coordinates, and covered50 depends on them. Under `docs/lockbox.md` ("If a default changes after scoring … a second lockbox would be needed"), the first lockbox's figure no longer describes the shipped calls. It stays the record of the earlier defaults and is not re-scored.

The new stages were tuned on the 42-genome sample and checked on the 18-genome sample and the literature pilot, so none of those is held out for them. This lockbox is.

### Selection (fixed)

`scripts/lockbox_select.py select`, run on ARC by `scripts/lockbox2_build.slurm`:

- **Source:** the first lockbox's NCBI summary, `analysis/lockbox/assembly_summary_refseq.20260929.txt`, refused unless its SHA-256 is `3459db2b5ab91279a26cc109b2917abe4a11e44105f5bf2d444d40011651e931`. The same file means the same eligibility rule over the same assemblies.
- **Eligible:** bacterial RefSeq `Complete Genome` assemblies, `latest` version, whose genus is absent from the corpus (train, val and test) **and** from `analysis/lockbox2/exclude_genera.txt`: every genus of the first lockbox (150), of the LAMBDA clean subset (15) and of the G5b Philympics set (6), 171 in all. Those three sets informed or measured earlier decisions. "uncultured" and similar genus names are excluded, as before.
- **One per genus:** reference genome, then representative genome, then any, ties broken by accession.
- **Sample:** 150 genera drawn with `numpy.random.default_rng(20261021)`, or all of them if fewer are eligible. 338 genera were eligible for the first lockbox and 150 of them were drawn, so at most 188 remain before the LAMBDA and G5b exclusions. The list's SHA-256 is recorded in `analysis/lockbox2/lockbox_selection.json`.
- **Truth:** geNomad 1.9.4 `annotate` + `find-proviruses`, DB 1.9, the corpus's settings. Truth is geNomad agreement, as everywhere else.

**Sample size.** The first lockbox had 202 elements in 150 genomes. A similar yield gives a 95% Wilson interval of about ±3 points at 96%. Strata under 10 kb will be thin again and are reported with their *n*, not as claims.

**What the agent may see before scoring:** the genome list, the selection table, and element counts per stratum. No model output on these genomes, as for the first lockbox.

### Frozen before scoring

Scoring happens **once**, after all of the following are final and committed. The commit is recorded here, and `scripts/lockbox2_run.slurm` (written at the freeze) refuses to run any other code.

- the checkpoints: `prophage_v2_pool_r1..r3`, the SHA-256s in `docs/research-integrity.md` C8 (bundled as `model/manifest.json`);
- the `predict` defaults: the three modes as in `docs/lockbox.md`, coverage 0.1, merge 50 kb / 0.1, `impute`, `scan`, edge snap 0.5; **refinement on** (5 kb windows, 250 bp step, threshold 0.6, α 0.6, ramp 15,000–60,000, split 4,000 / 0.05); **att search on** with the rule in `src/att.py` (search 25 / 10 kb, accept 10 / 5 kb, core ≤300 bp, ≤2 copies, entropy ≥1.5 bits, E ≤ 1e-6); edge ranges −0.055 / +0.494 × call length;
- the MD5 of every module the CLI loads, including the new `src/att.py`;
- the scorer, `scripts/lockbox_score.py --boundaries`, and this doc's endpoints.

### Endpoints (fixed)

- **Primary:** ensemble, Default mode, covered50 over all lockbox-2 elements, with a 95% Wilson interval and a genome-level bootstrap interval. Same definition as the first lockbox.
- **Secondary:**
  - the same for Exploratory and Conservative;
  - calls per genome per mode; covered recall (≥1 bp); conflated calls and absorbed elements; covered50 per size stratum, with *n*;
  - the gap to the 42- and 18-genome estimates (their covered50 is unchanged by the new stages: 0.994 / 0.972 / 0.944 and 1.000 / 1.000 / 0.962);
  - **boundary endpoints** (per mode): median IoU per element; the share of calls with att boundaries; the coverage of model-edge ranges (the geNomad edge inside `start_lo..start_hi` or `end_lo..end_hi`, on calls overlapping one element, contig-end edges excluded); att edges' distance to geNomad's (median, ≤100 bp, ≤1 kb); bp called and bp off every element.
- **Reported whatever the outcome.** No genome or element is dropped after scoring. A genome that fails to download or label is excluded *before* scoring, and listed.

#### Reporting rules (fixed)

- **Recall:** the release states "≥95%" only if the primary endpoint's Wilson lower bound is at least 0.95. Otherwise it reports the number and its interval, as it does for the first lockbox.
- **Edge ranges:** they are called 90% ranges only if the measured Default coverage's Wilson interval includes 0.90 or lies above it. Otherwise the release states the measured coverage instead.
- **Att edges** are compared with geNomad's, whose borders sit at gene edges, so this distance measures agreement, not att accuracy. Base accuracy is stated only from author-curated truth (the literature pilot, `docs/att-boundaries.md`), as confirmation after use.

### Amendments

#### Amendment 2026-10-21: composition edges added to the frozen defaults (before any selection)

Written the same day, before the genome list existed. After this doc was written, the user decided to ship composition edges in Default and Conservative (`docs/att-boundaries.md` "Composition edges"). The frozen defaults therefore also include `src/composition.py` (divergence ≥ 0.01, inward cap 0.15 × call length) with `--composition-edges auto`, and the edge ranges are 0.143 / 0.463 × call length in Default and Conservative, 0.055 / 0.494 in Exploratory. The scorer's run check expects `boundary_refinement` = `fine_pass_then_att_then_composition` in Default and Conservative and `fine_pass_then_att` in Exploratory. Endpoints and reporting rules are unchanged; the edge-range rule applies to each mode's own ranges.

#### Amendment 2026-10-21: composition inward ceiling (before any selection)

Written the same day, still before the genome list existed. The composition cap is min(0.15 × call length, 8 kb) (`composition.MAX_INWARD_BP`), chosen after the LAMBDA re-score showed a 22.8 kb move cutting off a curated element (`docs/att-boundaries.md` "LAMBDA under the 0.1.0 pipeline"). This value is part of the frozen defaults.

#### Amendment 2026-10-26: Exploratory composition edges and split-piece support (after selection, before any model output)

Written after the genome list and truth existed (build `848536`, Status above) and before any prophound output on a lockbox-2 genome. Only the selection table, the genome list and per-stratum element counts have been seen. Selection, truth, endpoints and reporting rules are unchanged.

- **Exploratory gains composition edges** (user decision 2026-10-23, `docs/att-boundaries.md` "Composition everywhere"): `composition.MODES` is all three modes, with the profile-divergence gate 0.005 in Exploratory (`MIN_JSD_EXPLORATORY`) and 0.01 in Default and Conservative. Edge ranges are 0.143 / 0.463 × call length in every mode. The scorer's run check now expects `boundary_refinement` = `fine_pass_then_att_then_composition` in all three modes (`lockbox_score.py` `BOUNDARY_STAGES`). Decided on the tuning, check, literature and LAMBDA samples, none of which shares a genus with this lockbox.
- **`seed_support` is counted on the extent the vote saw** (user decision 2026-10-26, `docs/att-boundaries.md` "Split-piece support"). A reported column only: calls and coordinates are unchanged (12 of 12 sample × mode rebuilds identical). No endpoint uses `seed_support`.
- **Execution.** `scripts/lockbox2_run.slurm` runs `prophound predict` from a `git archive` of the freeze commit, with the bundled model, instead of the ARC tree's `src/` as the first lockbox did: that tree predates the package layout. Its preflight checks every file of the archive against `analysis/lockbox2/frozen_md5.txt` (whose SHA-256 is in the launcher), the genome-list SHA-256, every FASTA, and that the example genome gives its expected calls on the V100 before any lockbox genome runs.

#### Amendment 2026-10-27: commit hashes after a history rewrite (no code change)

History rewrite, before the first push. To keep three internal files off GitHub, they were stripped from the unpushed commits with `git filter-repo`; the 32 commits already on GitHub are unchanged. Every rewritten commit got a new hash, and another when the commits were then re-signed, so the hashes cited above refer to the pre-rewrite history (kept in the local backup bundle). The map gives the final, signed hashes. Code is unchanged: `git diff` of `src/` and `examples/` between each old and new commit is empty, so `frozen_md5.txt` and every MD5 list verify both. SLURM launchers and ARC paths (for example `prophound_0d71e2f`) keep the old names, because they record what ran. References in this doc to `docs/release-checklist.md`, `SESSIONS.md` and `CLAUDE.md` point to internal files that are not published; the rules they stated are in `docs/project-principles.md`. Map for the hashes this doc cites:

| cited | now |
|---|---|
| `0d71e2f` | `24f6c44` |
| `382efba` | `26880c3` |
| `ad5b28b` | `09c7756` |

The freeze is `24f6c44` (was `0d71e2f`). Run `851599` used the archive of `0d71e2f`, whose `src/` is identical.

#### Amendment 2026-10-27: help text and docstrings updated after scoring (no output change); the lockbox figure stands

After the result, two hash-checked modules changed in text only:
- `src/options.py`: the `--mode` help's call counts now read about 37 / 12 / 5 per genome, this lockbox's figures (they were about 32 / 11 / 5, the first lockbox's).
- `src/att.py`: two docstrings no longer call the edge ranges "90%", since they held 0.78 here.

Checks:
- Removing docstrings leaves `att.py`'s syntax tree identical to the scored code's, and `options.py`'s differs only in that help string.
- The fast tests (198) and the smoke test pass.
- CPU reference matrix against `support_fix_cpu`, which ran the scored code: **30 of 30 identical** (`analysis/refout/helptext_cpu`, 2026-10-27 12:52–14:52).

Nothing that produces a call changed, so the result above describes the shipped code.

#### Amendment 2026-11-02: package moved to the repository root (no shipped file changed); the lockbox figure stands

The repository layout now mirrors the user's GenomeProfiler repository: `src/prophound/` moved to `prophound/`, and the recipe to `bioconda/recipes/prophound/`. Paths above that start `src/prophound/` now start `prophound/`. The text above keeps them as written.

- **No shipped file changed** in the move: all 50 files of the package code and model are pure renames, so every hash-checked module has the same MD5 as before. Only dev scripts were edited, for their paths.
- The built wheel is byte-identical to the one before the move.
- The fast tests (198), the smoke test and the golden set (calls identical on the Apple GPU) pass.
- `scripts/lockbox2_run.slurm` and `analysis/lockbox2/frozen_md5.txt` still name `src/prophound/`, because they check the frozen code archive on ARC (`0d71e2f`), which has the old layout.

#### Amendment 2026-11-02: license headers and a pylint comment (comments only); the lockbox figure stands

Every tracked code file now carries the Apache-2.0 license header. Four files on this lockbox's MD5 list changed, in comments only:

- the three empty `__init__.py` files (`prophound/models/`, `models/prophage/`, `models/prophage/src/`) now hold the header and nothing else;
- `src/att.py` has `# pylint: disable=unsubscriptable-object` on one line, a pylint false positive (`best` is checked for `None` on the same line).

Checks:

- All four files' syntax trees are identical to the scored code's, so the compiled code is the same.
- The fast tests (198), the smoke test, the golden set (calls identical on the Apple GPU) and the example (calls identical) pass.
- The CPU reference matrix was not re-run: with identical syntax trees it cannot differ.

#### Amendment 2026-11-03: a warning for non-prokaryotic input (stderr only); the lockbox figure stands

After the result, the CLI gained a check for input that may not be prokaryotic (`docs/cli-and-packaging.md` "The non-prokaryote warning"; the panel is `analysis/robustness/domain_check.tsv`). It reads only contig lengths and contig ends. It prints a warning to stderr when eukaryotic telomere repeats sit at two or more contig ends or the input totals more than 17 Mb, and the run continues. `--no-domain-check` skips it.

- `src/domain.py` is new. `src/predict_genome.py` calls it before the windows are built, and `src/options.py` adds the flag. Nothing that produces a call changed.

Checks:

- CPU reference matrix against `helptext_cpu`: **30 of 30 identical**, output files and stderr (`analysis/refout/domain_cpu`, 2026-11-02 18:54 to 2026-11-03 14:17, paused overnight). The warning fires on none of the test inputs. Besides the usual keys, the compare ignored the model, config and input paths, which name `src/prophound/` in the reference (the checkpoint hashes are equal), and the OS release (macOS 27.0 to 27.2).
- **The matrix imported the main checkout, not the snapshot it was given.** Since the package moved to the repository root, `python -c` started from the root puts the checkout ahead of `PYTHONPATH`. The checkout held this change and differed from the snapshot only in comments and docstrings (the next amendment), so the result holds for this code. Earlier matrices are unaffected: before the move the root held no `prophound/`. `scripts/reference_outputs.py` and `scripts/benchmark.py` now run Python with `-P`, and the runner stops unless `prophound` imports from `--code`. Two cases re-run that way from the snapshot are identical.
- The new tests (`tests/test_domain.py`) and the fast tests (243, on the tree with the next amendment's change) pass, and so does the smoke test. The golden set passes on CPU, byte for byte, on the same tree.

#### Amendment 2026-11-03: citations of unpublished notes reworded (comments and docstrings only); the lockbox figure stands

Only the docs the README references are published. Comments and docstrings that cited the project's internal notes now cite a published doc or a result folder in `analysis/`, or nothing. Files on this lockbox's MD5 list changed as follows:

- **Modules the CLI loads:** `src/att.py`, `calls.py`, `composition.py`, `encoding.py`, `model_rc_cnn.py`, `multiscale.py`, `predict_genome.py` and `report.py` changed in comments and docstrings, and `src/options.py` in comments only. `src/domain.py`, added after scoring for the non-prokaryote warning, changed the same way.
- **The scorer,** `scripts/lockbox_score.py`, changed in its docstring and one comment.
- **Dev scripts, launchers and training configs** changed in comments, including the preflight banner of `scripts/lockbox2_run.slurm`. Three scripts also changed a string: a help text (`coverage_builder_checks.py`), a stderr notice (`evalkit.py`), and a descriptive value that `composition_baseline.py` writes to its results file. The tracked `analysis/composition_baseline/results.json` is unchanged.

Checks:

- With docstrings removed, every changed module's syntax tree is identical to the version before, and `options.py`'s is identical outright. So is the scorer's, so scoring is unchanged.
- 50 of the 53 changed dev scripts and three tests match the same way; the other three scripts differ only in the strings above. The 37 shell files are identical outside comments, and the 13 YAML configs load to the same values. (The commit guard and its test also changed, deliberately; neither is on the MD5 list.)
- The fast tests (243) pass.
- The smoke test and the golden set (CPU, byte for byte) pass. This is the same run as in the amendment above.
- The CPU reference matrix was not re-run for this change: with identical syntax trees it cannot differ.
