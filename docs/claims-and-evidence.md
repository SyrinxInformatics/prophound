# Claims and evidence

> **Sample roles** (`docs/research-integrity.md` C1, added 2026-09-29). In this doc, "held out" means held out from *training*. No sample here is held out from *tuning*. The 42-genome / 179-element sample is the **tuning sample**: thresholds, merging, the seed vote and the modes were chosen on it. The 18-genome / 52-element sample is **confirmation after use**: it chose the coverage builder, surfaced the ambiguity bug, and was in view for the 2026-09-29 mode decisions. Only the lockbox (`docs/lockbox.md`) is held out from both. Numbers below are as measured at the time and may predate the final modes (`docs/cli-and-packaging.md` "Seed vote per mode").

Living doc, updated in place. Built 2026-09-28. It says which novelty and improvement claims the v1 write-up can make, what evidence each needs, and what the tool itself must do so that its exploratory outputs can be trusted. The positioning argument and related tools are in `docs/related-tools.md`. Release tasks are in `docs/release-checklist.md`.

## The principle (user, 2026-09-28)

The tool should make a novelty or improvement claim confidently without overselling. It should **not rest on the claim alone.** It must be technically robust and easy to use, and its accuracy on *known* elements must be measured and published, so that its exploratory outputs can be considered with confidence. The order that follows from this: known accuracy first, then robustness and usability, then the improvement claims, and the concept last. The concept is only as credible as the measured layer under it.

## Release gate: new, or measurably better (user, 2026-09-28)

> "If this isn't doing anything new, it at least needs to be better."

The concept (rung 1 below) is new only as a *contract*, and the methods are not new. So v1 has to show it is **better than existing tools on at least one measured axis**, on independent truth or identical input, not only on agreement with geNomad. **Scope rule (user, 2026-09-28): an axis is in the gate only if measuring it will not cause major delay.** The axes are sorted by what they cost to measure.

### In the gate: all measured by one run, no new tools to install

One ensemble run of prophound (all three modes) on LAMBDA's 39 clean genomes, scored locally against the published predictions of the other tools. Prerequisites: a 21-mer containment check against training positives (reuse `small_element_spikein.py`'s `disjoint`) and confirming LAMBDA's coordinate frame. The cost is about one GPU job the size of the 42-genome runs (1–2 h). Scoring against the other tools' predictions is plain Python and **can run on ARC** in the same job or a CPU follow-up (user, 2026-09-28), so it adds no local turnaround. The LAMBDA JSONs are about 20 MB and are pushed alongside.

| # | axis | bar to clear (LAMBDA clean subset) | outlook |
|---|---|---|---|
| G1 | **Complement to geNomad** | geNomad ∪ prophound gains **more than +8** curated elements over geNomad alone (0.809, 123/152). +8 is the best existing pairing (PHASTER), at +3.9 calls/genome. | **Plausible**, most likely in Exploratory. It is the direct test of the concept. Report the extra calls/genome with it. |
| G2 | **Sequence-only, at matched cost** | top-*k* calls per genome (k matched to ~3.5) reach covered50 **above 0.727** (best genomic LM, NTv2; Evo2 0.724) | **Plausible.** It is trained on integrated elements in host context, whereas the LMs trained on free phage genomes. |
| G3 | **Overall curated recall** | covered50 **above 0.818** (PHASTER; PIDE 0.812, geNomad 0.809), with calls/genome | **Unlikely at Default**, because the model inherits geNomad's misses. Possible in Exploratory, at about 7× the calls/genome, and must be stated that way. |
| G4 | **Footprint** | install size versus geNomad's database: `genomad_db_v1.9.tar.gz` is **842 MB compressed** (Zenodo 14886553; measure the unpacked size). Our side is 3 checkpoints plus the package. | **Likely.** It counts only alongside a G1–G3 recall figure, since a small tool that finds less is not better. |

### Time-boxed: in the gate only if setup goes smoothly

| # | axis | bar | cost and fallback |
|---|---|---|---|
| G5 | **Fragmented input** | on identical MAG-like shreds of the 39 clean genomes, prophound retains more complete-genome recall than geNomad and Jaeger `-p` (`--lc` lowered, and stated) | **The most likely win**, but the only axis that needs installing other tools. ARC lacks MMseqs2 and Prodigal (preflight), so run the comparators locally, where both are conda/pip-installable. **Time-box: about one day for installs.** If it slips, fall back to a self-referential figure ("retains X% of its complete-genome recall"), which is measured already on eval genomes (covered 0.939) and needs no comparator. That is not a "better" claim. |
| G6 | **CPU runtime** | wall-clock per genome versus geNomad on the same machine | Comes free with G5's local geNomad install. Without it, use no runtime comparison; published figures come from other hardware. |

### G5 preregistration (2026-09-30, before any tool output on the shreds is scored)

Fixed before scoring, like the LAMBDA bars. Changes after this point are dated amendments below it, never edits.

**What was seen beforehand.** On one genome (GCA_000143885.2), Jaeger was run on the complete chromosome and on a single MAG-like shred. That was a format and runtime check: it produced 1 and 2 calls overlapping its one curated element. Nothing else was run or scored. prophound's complete-genome LAMBDA result (`822275`) was already known.

**Correction to the G5 row above.** ARC does have geNomad: the labelling env (`mge_explorer`: geNomad 1.9.4 with DB 1.9, which the lockbox build used) runs it. The "ARC lacks MMseqs2 and Prodigal" note referred to the prophound env. So geNomad runs on ARC in the exact label version, and only Jaeger runs locally.

**Input.** The 39 LAMBDA clean-subset chromosomes (149 curated elements), fetched as in `lambda_fetch.slurm`. MD5s for these and every shred are in `analysis/g5/input_md5.tsv`, and ARC refuses to run on anything that differs. Shreds come from `scripts/shred_genomes.py --dist lognormal --sigma 1.0 --min_contig 1500`. Fragment lengths are log-normal, with the length-weighted median set to the target N50. Fragments under 1.5 kb are dropped, as an assembly or binning cutoff would drop them. Breakpoints are uniform: real assemblies break preferentially at repeats, which this does not model, and the write-up says so.

| condition | N50 target (realized) | fragments | bp kept | role |
|---|---|---|---|---|
| `complete` | — | 39 | 1.000 | denominator |
| `mag_s0`, `mag_s1`, `mag_s2` | 20 kb (20.0–21.0 kb) | 14,570–14,971 | 0.995 | **primary**, pooled over the three shred seeds |
| `draft_s0` | 100 kb (100.0 kb) | 3,164 | 1.000 | secondary |

The N50 values are chosen as representative of short-read MAG and draft-isolate assemblies. They are not calibrated to any catalogue.

**Tools**, identical settings on every condition:
- **prophound:** the 3-checkpoint ensemble (SHA-256s as in C8) with shipped defaults (`--short-contigs scan`, `impute`). Modes are rebuilt offline with the shipped seed vote, as for LAMBDA. Contigs under 3 kb are never scanned; that is part of the tool. The complete condition reuses `822275`'s tracks.
- **geNomad 1.9.4, DB 1.9, `end-to-end`,** default settings (plus `--cleanup`, which does not change results). Calls are every `virus_summary` row: proviruses at their coordinates, and whole-contig viruses at full length. On fragments, an element that fills its contig is reported as a virus rather than a provirus, and scoring only proviruses would penalise geNomad for that. Secondary: the find-proviruses table alone (the labelling pipeline).
- **Jaeger 1.26.5,** CPU TensorFlow, `predict -p --lc 1500`. `--lc` is lowered from its 500 kb default to the shred minimum; otherwise it scans no fragment. This is stated with every Jaeger number. Calls are non-rejected prophage regions [sstart, eend] (raw_start/raw_end where no att site was found), plus contigs classified `phage` at full length. **Primary model: `jaeger_57341_1.5M_fragment`.** The bundled `default` model is deprecated in 1.26, and the model the README recommends (`jaeger_38341_1.4M_fragment`) is not on the download server. Secondary: `default`.

**Scoring** (`scripts/g5_score.py`). Calls are mapped back to chromosome coordinates. Covered50 is computed against LAMBDA truth on all 39 genomes. Retention = shredded covered50 / complete covered50, per tool. MAG is pooled as elements × 3 shreds. Intervals are 95% paired genome bootstrap (10,000 resamples, seed 20260930), with the same resamples for every tool.

**Primary endpoint:** prophound **Default**, MAG condition, retention against geNomad (virus_summary) and against Jaeger (57341_1.5M).

**Bar: G5 clears only if, against both comparators:**
1. prophound's retention is higher, with the paired interval of the difference above 0; and
2. prophound's absolute MAG covered50 is not lower (point estimate).

A higher retention on a lower absolute recall is not a win and is reported as such. Calls/genome and prophound at each comparator's per-genome call count (Exploratory calls ranked by max_prob) are always reported alongside, because the LAMBDA result showed that much of a raw margin can be cost.

**Secondary, reported whatever the outcome:** Exploratory and Conservative; the draft condition; geNomad find-proviruses; Jaeger `default`; the 34 not-in-corpus genomes; elements under 10 kb; the reachable ceiling (elements ≥50% present in the kept fragments); geNomad and prophound wall-clock per run (a start on G6; hardware differs, so this is not the G6 comparison).

**Wording set by outcome.** Cleared: *"on MAG-like fragmented input (N50 20 kb), prophound kept X of its complete-genome recall on LAMBDA's curated prophages, against Y for geNomad and Z for Jaeger, at N calls per genome"*. Not cleared: report the three retention figures plainly, with no comparative claim.

**Runs:** `g5_shred.slurm` (**`822964`**) → `g5_prophound.slurm` (GPU, **`822965`**) and `g5_genomad.slurm` (CPU, **`822966`**) on ARC, submitted 2026-09-30; Jaeger locally with `scripts/g5_jaeger.sh` in `metadata/g5-comparators/` (pixi).

**Amendment 2026-09-30 (same day, before any scoring): Jaeger moves to ARC.** The user installed Jaeger 1.26.5 in `mge_explorer_gpu`, so the Jaeger arm runs as `scripts/g5_jaeger.slurm` on `compute1` (**`822991`**). It uses the same models, the same flags (`-p --lc 1500 --cpu`), the same output layout and collector, and the same MD5-checked inputs. The local run (`scripts/g5_jaeger.sh`) was stopped after 17 of 195 primary-model runs, all on complete genomes. Its outputs were deleted unscored, and the local script stays only as a fallback. The design, endpoint and bar are unchanged. The TensorFlow build differs from the local pip wheel, which is not expected to change predictions.

**Amendment 2026-09-30, second (before any scoring): Jaeger's primary model runs at `--fsize 3077`.** `822991` failed without producing any output, for two reasons: a launcher bug (a missing log directory) and `pyrodigal-gv` missing from the ARC env. The smoke log carried Jaeger's own warning that `jaeger_57341_1.5M_fragment` was trained on 3,077 nt fragments and "may degrade or collapse" at the 2,000 nt default, and that `--fsize 3077` is preferred. Running a comparator against its own documented guidance would handicap it, so the primary model uses `--fsize 3077`. The secondary `default` model keeps its defaults. Contigs shorter than the window are left at Jaeger's default `--min-len` (= `--fsize`). The launcher now also rejects any run whose log contains an ERROR line, because Jaeger exits 0 when its prophage step fails. `pyrodigal-gv` 0.3.2 was installed, matching the local env, and the job was resubmitted as **`823003`**.

### G5 result (`822964`–`822966`, `823003`; scored once on 2026-10-01)

All runs were clean: every `.err` was empty and every job exited 0. That covers 156/156 prophound, 195/195 geNomad and 390/390 Jaeger runs, with no failures. All shreds were byte-identical across tools, and checkpoint SHA-256s match C8. Every element is reachable in every condition (≥50% of its bp lies in the kept fragments). Scored with `scripts/g5_score.py score`; outputs in `analysis/g5/results.{json,csv}`.

**Covered50 on LAMBDA curated truth, 39 genomes / 149 elements.** MAG is pooled over three shred seeds. Intervals are 95% paired genome bootstrap.

| tool | complete | MAG (N50 20 kb) | **retention** | MAG calls/genome | draft (N50 100 kb) | draft retention |
|---|---|---|---|---|---|---|
| prophound Exploratory | 0.946 | 0.908 | 0.960 | 40.4 | 0.946 | 1.000 |
| **prophound Default** | 0.933 | **0.855** | **0.916 (0.844–0.968)** | 17.6 | 0.933 | 1.000 |
| prophound Conservative | 0.899 | 0.781 | 0.868 | 11.0 | 0.893 | 0.993 |
| **geNomad** (end-to-end) | 0.812 | **0.830** | **1.022 (0.991–1.065)** | 14.4 | 0.819 | 1.008 |
| **Jaeger 57341_1.5M** | 0.859 | **0.848** | **0.987 (0.969–1.009)** | 22.6 | 0.852 | 0.992 |
| geNomad find-proviruses only | 0.819 | 0.273 | 0.333 | 3.2 | 0.745 | 0.910 |
| Jaeger `default` model | 0.832 | 0.617 | 0.742 | 11.4 | 0.262 | 0.315 |

**Verdict: G5 does not clear. Its claim is refuted, not just unsupported.**
- Retention, prophound Default minus each comparator: geNomad **−0.106** (−0.215 to −0.030), Jaeger **−0.071** (−0.148 to −0.013). Both intervals lie below zero.
- Absolute MAG covered50: +0.025 over geNomad (−0.036 to +0.071) and +0.007 over Jaeger (−0.048 to +0.052). These are ties.
- Matched cost (Exploratory cut to the comparator's per-genome call count): MAG 0.81–0.85 against geNomad's 0.82–0.84, and 0.83–0.86 against Jaeger's 0.84–0.85. Also ties. prophound's complete-genome lead over both tools does not survive fragmentation. On LAMBDA it was mostly cost already.

**Why (post-hoc, exploratory; not part of the preregistered analysis).** Split by how many MAG fragments overlap an element (elements × 3 seeds), Default vs geNomad vs Jaeger:
- 1 contig: 0.838 / 0.635 / 0.784 (n=74);
- 2: 0.886 / 0.829 / 0.814 (70);
- 3–4: 0.844 / 0.874 / 0.859 (135);
- 5+: 0.857 / 0.881 / 0.881 (168).

prophound leads when an element survives in one or two pieces. It trails once the element is in three or more, and on 60 kb+ elements (0.767 against 0.800 and 0.800). The comparators classify whole contigs: a fragment lying entirely inside an element is reported as a virus at full length. That is why geNomad's calls rise from 4.2 to 14.4 per genome and its retention exceeds 1. prophound has no contig-level call and never merges across contigs. Skipping contigs under 3 kb is not the cause: the 22 elements only prophound misses have a median of 3% of their bp in such pieces, against a median of 0.47 of the element in their largest piece. A contig-level "this whole fragment is phage" output is a plausible post-v1 fix, untested.

**What the release may say.**
- **Not supported:** any claim that prophound handles fragmented input better than geNomad or Jaeger.
- **Supported, self-referential:** *"on MAG-like fragmented input (N50 20 kb), prophound Default kept 0.92 of its complete-genome recall on LAMBDA's curated prophages (0.855 covered50), at 17.6 calls per genome; on draft-like input (N50 100 kb) it kept all of it"*, stated alongside geNomad (1.02) and Jaeger (0.99) on the same input.
- **Supported, design fact only:** prophound scans fragments by default, while Jaeger's prophage mode needs `--lc` lowered from 500 kb. This is a usability difference, not a performance one.
- The find-proviruses-only collapse (0.273) is a statement about geNomad's provirus step alone, not about geNomad as users run it, and is not used as a comparison. It does show why the comparator had to be end-to-end.
- The Jaeger `default` model's collapse on draft input (0.262) is unexplained and secondary. It is reported, not interpreted.

**Runtime, toward G6 (different hardware, so not the G6 comparison):** geNomad 340 s per complete genome (4 CPU threads); Jaeger 57341_1.5M 189 s and `default` 55 s (8 CPU workers); prophound about 22 s per genome on one V100 for the MAG shreds.

### G5b preregistration (2026-10-01, before any fix exists)

**Why.** G5 did not clear (above), and its post-hoc look used the LAMBDA shreds. The user decided on 2026-10-01 to fix the fragmented-input weakness and to confirm it on data no part of the work has touched. **The G5 result stands as reported.** A LAMBDA re-run with the fix is reported as post-hoc only and is never evidence for G5b.

**Protocol, in order:**
1. **Development on corpus shreds only.** The fix is built and tuned on the 42 size-stratified eval genomes, shredded exactly as in G5, with geNomad truth (`scripts/g5dev_prophound.slurm`, `analysis/g5dev/`). It is checked on the 18 confirmation genomes. Only prophound runs; its target is its own retention (shredded / complete covered50) and calls/genome. No LAMBDA or Philympics output is looked at during development.
2. **Freeze.** Before G5b is scored, an amendment below this section records the git commit, the fix's parameters and its default-on/off state. prophound runs on G5b only after that amendment.
3. **Comparators** may run on G5b before the freeze (they do not depend on the fix). Their output is collected but not looked at until prophound's G5b run exists.
3a. **Jobs** (submitted 2026-10-01): development `824148`; G5b shred `824149`, geNomad `824150`, Jaeger `824151`.
4. **Score once**, with `scripts/g5_score.py score --truth_tsv analysis/g5b/g5b_truth.tsv --genomes_tsv analysis/g5b/g5b_genomes.tsv`.

**Genome set, fixed now:** the 9 Philympics genomes that are outside our corpus *and* not among LAMBDA's 39. Philympics is pinned to `linsalrob/ProphagePredictionComparisons@dacfb153` (GenBank via Git LFS; `genbank/prophageCoords.tsv`). The 9 are *Achromobacter denitrificans* PR1, *Burkholderia cenocepacia* PC184, *Campylobacter jejuni* (GCA_002082995.2), *Cyanobacterium* sp. HL-69, and five *Paracoccus* (*aminophilus* JCM 7686, *denitrificans* PD1222, *sanguinis* 5503, sp. SCN 68-21, *yeei* TT13). The 4 genomes unmapped to our corpus are excluded by rule, because their corpus status cannot be verified cheaply. MD5s for the GenBank files, the coordinates file, the FASTAs and all shreds are in `analysis/g5b/input_md5.tsv`.

**Truth:** Philympics rows with ≥5 genes or ≥5 kb, giving **32 of 36 rows**. Coordinates are converted 1-based inclusive → 0-based half-open. One fuzzy bound (`>1753287`, *C. jejuni*) is taken at face value. Built by `scripts/g5b_prepare.py`. Every record of each GenBank file is kept, and headers are bare accessions. Two genomes are already drafts (*P. sanguinis*: 265 records; *Paracoccus* sp. SCN 68-21, a MAG: 87), so their "complete" condition is the deposited assembly.

**Conditions:** identical rule to G5 (log-normal, sigma 1.0, ≥1.5 kb; MAG N50 20 kb seeds 0–2, draft N50 100 kb seed 0). The drafts cap fragment length, so realized N50 is lower: MAG 18.5–19.4 kb, draft 75 kb. 3,145–3,247 MAG fragments.

**Tools:** exactly as in G5 with its amendments. geNomad 1.9.4 / DB 1.9 end-to-end (`virus_summary`). Jaeger 1.26.5, `57341_1.5M` at `--fsize 3077` (primary) and `default` (secondary), `-p --lc 1500 --cpu`. prophound is the frozen post-fix version.

**Endpoint and bar: unchanged from G5.** prophound Default, MAG pooled over three seeds. G5b clears only if prophound's retention exceeds both geNomad's and Jaeger's with paired genome-bootstrap intervals above 0, *and* its absolute MAG covered50 is not lower than either. Calls/genome and matched cost are reported alongside.

**Power, stated in advance.** There are 32 elements in 9 genomes, against LAMBDA's 149 in 39, so intervals will be roughly twice as wide. Clearing needs a retention lead of about 0.13 over both tools. On LAMBDA, prophound trailed by 0.07–0.11 before any fix. **An inconclusive G5b is the likely outcome** unless the fix is large. If it is inconclusive, the release gets the fix with self-referential retention figures and no comparative claim, as the G5 verdict already requires. Five of the 9 genomes are one genus, and the write-up says so.

### Measured by the same run, where the tool is probably not better (report plainly)

| axis | existing level (LAMBDA clean) | our reference point | reading |
|---|---|---|---|
| **Boundary precision** | median IoU of detected elements: geNomad **0.921**, PIDE 0.912, VIBRANT 0.885, PHASTER 0.831. IoU ≥ 0.5 on 0.72 / 0.76 / 0.68 / 0.70 of elements | against geNomad truth, our median IoU is about 0.69–0.74 (`docs/scientific-validity.md`) | **Probably worse.** Annotation tools snap edges to genes, tRNAs and integrases; ours come from 1 kb bins of the score track. State it as a known limitation. It is a real cost of skipping gene calling. |
| **Small elements (<10 kb)** | 18 of 154 clean LAMBDA elements | 5–8 kb covered50 0.857 on geNomad truth (n=7) | Known weakness, already in the limitations. |
| **Calls/genome at the same recall** | 3.4–4.5 for most tools | Default 14.5, Conservative 9.4 | Costlier by design. Frame it as follow-up budget, never hide it. |

### Not in the gate (would cause delay)

- **Ambiguity robustness versus other tools.** It needs spiked input run through every comparator. Our own dose-response (`821396`) is reported as a robustness fact, not a comparison.
- **Contributing prophound to the Philympics framework** (a Snakemake rule and a PR to an external repo, on their timeline). Post-v1. *Scoring* Philympics' 23 clean genomes ourselves can be added to the G1–G4 job if it costs less than a day, since the published CDS-level results are in its repo (`jupyter_notebooks/all_tptn.tsv`).
- **Anything requiring retraining** (small elements, boundary head): deferred past v1 already.

### Decision rule

- **If any of G1–G3 or G5 clears its bar**, the write-up leads with the axes that cleared: each with its bar, our number, the calls/genome and the n. G4 and G6 go alongside, never alone.
- **Every axis measured is reported**, including those where the tool is not better (boundary precision, small elements, cost per genome). Reporting the losses is what makes the wins credible.
- **If none clears**, v1 cannot claim to be new or better, and the release decision goes back to the user.

## Two layers of "known accuracy"

This is the finding that most changes how results should be written up.

| layer | what it measures | current number |
|---|---|---|
| **Agreement with geNomad** | how well the tool reproduces its label source on genomes held out from training | covered50 0.972 (0.949–0.994; ensemble, final Default, n=179). This is the **tuning sample**; the held-out figure comes from the lockbox |
| **Agreement with curated truth** | how well it finds prophages that experts annotated | **not yet measured**. geNomad itself reaches only **0.81** covered50 on LAMBDA's clean subset (below) |

geNomad covers about 81% of expert-curated prophages. A tool that matches geNomad at 96% will therefore land somewhere near 80% on curated truth. The 0.961 figure measures agreement, not accuracy. **Publishing it alone would oversell the tool, and the release must report both layers.** Curated truth is the number that tells a user how much of what is really there the tool finds.

## External curated benchmarks available now (checked 2026-09-28)

Two community gold standards exist. Their genome lists were matched against our corpus by assembly accession. Overlap tables: `src/prophound/models/prophage/analysis/external_benchmarks/`.

| benchmark | genomes | curated elements | not in our corpus | in our test split | in train/val (exclude) |
|---|---|---|---|---|---|
| **LAMBDA** (2026; Casjens gold standard, PHASTER manual annotations, Gauthier et al.) | 76 | 344 regions | 34 genomes (118 regions) | 5 | 37 |
| **Philympics** (2021; 78 manually curated genomes, CDS-level `is_phage` tags) | 73 files with coordinates (69 mapped) | 275 rows (includes 1-gene fragments) | 17 genomes (77 rows) | 6 | 46 |

**LAMBDA is the stronger instrument.** Its per-genome files include the curated truth, **every published tool's predicted regions** (geNomad, PHASTER, VIBRANT, VirSorter2, Phigaro, PhiSpy, PhageBoost, PIDE, and eight genomic language models) and 292 off-truth candidate regions that its authors classified from PHROG gene evidence. A head-to-head comparison therefore needs only an prophound run. The other tools are already scored on the same genomes.

### Amendment 2026-09-29, before any prophound run on LAMBDA: corrected baseline and bars

Rebuilt from pinned sources: `lambda-benchmark@d9cbf81` per-genome JSONs, and Zenodo record 21355242 (`LAMBDA_v1`, `metadata/prophage_reference_locations.csv`). Built by `scripts/lambda_prepare.py`, scored by `scripts/lambda_score.py`; outputs in `analysis/external_benchmarks/lambda/`. Three errors in the 2026-09-28 table below:
1. **Duplicate truth.** The dashboard JSON lists some elements twice (Photorhabdus TTO1: 25 entries, 20 unique). The clean subset has **149 curated elements, not 154**. The reference table is the deduplicated source, and the two agree as sets for all 39 genomes.
2. **"Ran" vs "called nothing".** Each tool is now scored only on genomes it was run on (listed in the JSON's `metrics`). geNomad ran on 37 genomes, not 38.
3. **Language models at one segment size.** Every model is now scored at 2k, 4k and 8k, and its best counts. That is conservative for us.

**Coordinate frame, confirmed for all 39:** LAMBDA scores the one sequence carrying the curated elements, the chromosome in every case. `genome_length` equals that sequence's length rounded up to 1 kb. prophound is therefore run on that sequence only, not the whole assembly.

Corrected covered50 (95% genome bootstrap), clean subset:

| tool | genomes | elements | covered50 | calls/genome |
|---|---|---|---|---|
| PHASTER | 38 | 148 | 0.824 (0.672–0.947) | 3.9 |
| PIDE | 38 | 148 | 0.824 (0.665–0.955) | 8.4 |
| geNomad | 37 | 146 | 0.822 (0.675–0.950) | 4.6 |
| VIBRANT | 37 | 146 | 0.788 (0.633–0.928) | 3.6 |
| EVO2 4k (best genomic LM) | 39 | 149 | 0.785 (0.626–0.920) | 4.8 |
| VirSorter2 | 18 | 101 | 0.762 | 3.9 |

The full table, all 37 tool labels, is in `baseline_scores.csv`. **Complement:** the best geNomad pairing gains **+5** elements, a tie among PHASTER, EVO2 4k, EVO2+SAE, Generanno 4k, PhageBoost and VirSorter2.

**Bars, corrected (supersede those below):**
- G1: complement **more than +5** (was +8).
- G2: matched-cost recall **above 0.785** against the best LM (was 0.727, NTv2 at one size).
- G3: overall covered50 **above 0.824** (was 0.818).

The intervals are wide because a few genomes carry many elements (Photorhabdus 20, EDL933 22). Gains are read through the paired comparison and the intervals, not from point estimates alone.

**Containment against training (`822260`, 2026-09-29).** Of the 149 elements:
- 36 have ≥10% of their canonical 21-mers in training-split positives;
- 29 have ≥50%;
- 10 have ≥90%;
- the median is 0.003.

About a quarter of LAMBDA's curated elements therefore closely resemble prophages the model was trained on, even though no LAMBDA genome is in the training split. Every LAMBDA figure is reported for the **low-containment subset (<0.1)** alongside the full set; `lambda_score.py --containment` does this. The coordinate note above holds: all 39 downloaded sequences have exactly the length NCBI reports, and LAMBDA's `genome_length` is that length rounded up to the next kb plus 1 kb.

### LAMBDA result (`822275`, 2026-09-29, final modes)

The ensemble ran on all 39 chromosomes (exit 0, checkpoint SHA-256s verified) and was scored locally with `scripts/lambda_score.py` (`analysis/external_benchmarks/lambda/scores.{json,csv}`). No call on the two N-bearing chromosomes is more than 2.9% N, so the offline scoring equals the current CLI's.

**Covered50**, with 95% genome bootstrap intervals:

| subset | prophound Exploratory | **Default** | Conservative | geNomad | PHASTER | best LM (EVO2 4k) |
|---|---|---|---|---|---|---|
| all 39 (149 elements) | 0.946 | **0.933 (0.874–0.983)** | 0.899 | 0.822 | 0.824 | 0.785 |
| **not in corpus, 34 (113)** | 0.938 | **0.920 (0.859–0.983)** | 0.885 | 0.791 | 0.795 | 0.779 |
| not in corpus, containment <0.1 (103) | 0.942 | **0.922** | 0.883 | 0.780 | 0.784 | 0.796 |
| test split, 5 (36) | 0.972 | 0.972 | 0.944 | 0.917 | 0.917 | 0.806 |
| under 10 kb, all (13) | 0.846 | 0.846 | 0.769 | 0.500 | 0.462 | 0.154 |
| calls/genome | 39.9 | 14.4 | 7.1 | 4.6 | 3.9 | 4.8 |

**Gates** (bars from the pre-run amendment):
- **G1 complement, bar >+5: cleared.** geNomad plus prophound Default gains **+16** elements (+18 Exploratory, +12 Conservative) at +14.6 extra calls/genome. The paired difference against geNomad, on the 37 genomes both ran, is **+16 / −0**: +0.110, interval 0.018–0.210.
- **G2 matched cost, bar >0.785: cleared narrowly.** Exploratory's calls, ranked by max_prob and cut to EVO2 4k's per-genome count, reach **0.805**, against its 0.785. At geNomad's cost it **ties** (0.822 vs 0.822). At PHASTER's lower cost it is **below** (0.797 vs 0.824).
- **G3 overall recall, bar >0.824: cleared** at every mode, but at 1.5–3.7× the competitors' call counts.

**Caveats that set the wording:**
- **Much of the raw margin is cost.** Default calls 9.8% of each chromosome, geNomad 3.5%. Rotating each tool's calls to random positions gives a chance covered50 of 0.144 for Default and 0.050 for geNomad. **Recall above chance is 0.787 against 0.772** (Conservative 0.808, PHASTER 0.771). Per unit of cost, prophound is level with the annotation tools. It does not beat them.
- **The gains are concentrated.** 12 of the 16 elements come from two genomes (Photorhabdus TTO1, 9; *P. aeruginosa* LESB61, 3), and only 5 genomes gain anything. Two are in the test-split EDL933.
- **The gains are not training look-alikes.** All 16 have containment ≤0.16 against training positives, 13 of them <0.01. Elements under 10 kb favour prophound, but n = 13 overall and 8 not in the corpus.
- The sample is 39 genomes, dominated by a few with many elements. Intervals overlap for every pairwise comparison except the paired one against geNomad.

**Reading for the release:**
- **Supported:** "run alongside geNomad, prophound recovered 16 additional curated prophages on LAMBDA's clean subset, at about 15 extra candidate calls per genome".
- **Supported:** it reaches the highest raw curated recall of the tools compared, *at a higher call cost*.
- **Not supported:** "more accurate than" any tool. At matched cost it ties geNomad and trails PHASTER.

**Caveat specific to prophound:** 5 of the 39 genomes are in our test split, including EDL933 (22 elements). EDL933 is also in the literature pilot, and was studied closely in the ambiguity work. The primary row is therefore the **34 genomes not in the corpus**.

### LAMBDA baseline on the clean subset (39 genomes, 154 curated regions, 18 under 10 kb) — superseded 2026-09-29 (see amendment above)

Computed 2026-09-28 from the published per-genome predictions, using our own metric (covered50: ≥50% of the element's bp under calls) so that every tool is scored the same way:

| tool | covered (≥1 bp) | **covered50** | calls/genome |
|---|---|---|---|
| PHASTER | 0.838 | **0.818** | 3.8 |
| PIDE | 0.844 | 0.812 | 8.5 |
| geNomad | 0.816 | 0.809 | 4.5 |
| VIBRANT | 0.789 | 0.783 | 3.6 |
| VirSorter2 (18 genomes only) | 0.764 | 0.755 | 3.9 |
| NTv2 (best genomic LM) | 0.747 | 0.727 | 3.4 |
| Evo2 | 0.730 | 0.724 | 3.7 |
| Phigaro | 0.766 | 0.701 | 4.4 |
| PhiSpy | 0.745 | 0.667 | 5.5 |
| PhageBoost | 0.753 | 0.644 | 8.7 |

**Complement baseline**, which tests "run alongside annotation tools". geNomad alone covers 123 of 152 elements (0.809; 38 genomes with geNomad output). Adding a second tool gives:

| geNomad + | covered50 | elements gained | extra calls/genome |
|---|---|---|---|
| PHASTER | 0.862 | +8 | +3.9 |
| VIBRANT | 0.842 | +5 | +3.6 |
| Evo2+SAE | 0.842 | +5 | +5.3 |
| NTv2 | 0.836 | +4 | +3.4 |
| PIDE | 0.822 | +2 | +8.7 |
| union of all 7 annotation-type tools | 0.908 | +15 | — |

These are the bars. **prophound's complement claim holds only if it gains more than +8 elements over geNomad** (the best existing pairing), and the gain is reported with its extra calls/genome.

### Caveats that apply to any number from these benchmarks

- **"Not in corpus" is genome-level only.** Our corpus is one RefSeq genome per genus, so a LAMBDA genome can share a genus, and possibly prophage relatives, with a training genome. Run the 21-mer containment check that `small_element_spikein.py` uses (`disjoint`) against training positives, and report the fully disjoint subset alongside.
- **Coordinate frame.** Each LAMBDA JSON gives a single `genome_length`. How multi-contig assemblies are laid out (concatenated or chromosome only) must be confirmed before scoring.
- Other tools' outputs are taken as LAMBDA published them. Versions and settings are theirs, and this should be stated.
- Philympics rows include 1-gene "prophages" of under 1 kb. Fix and state a filter (for example ≥5 genes or ≥5 kb) before scoring, and report Philympics' own CDS-level metric as well, since that is how its published tool results are stated.
- Our covered50 differs from LAMBDA's segment-level MCC and from Philympics' CDS-level F1. For our own table, rescore every tool with one rule, as above. Quote published figures only in their own units.

## The claim ladder

Each claim gives the evidence it needs, its status, and the strongest wording that evidence allows. Nothing moves up a rung without its evidence.

### Rung 1: claimable now

| claim | evidence | wording allowed |
|---|---|---|
| Design concept: recall-first screen, candidates, modes priced in calls/genome, runs alongside annotation tools | the design; survey in `docs/related-tools.md` | "designed as". Not "first". Being recall-leaning is not unique: Philympics characterizes VirSorter, VirSorter2 and DBSCAN-SWA that way, and geNomad has `--conservative`/`--relaxed` presets. What is distinct is stating the recall target explicitly and pricing each mode in follow-up cost. |
| No gene calling and no reference database at run time | code | plain statement |
| Agreement with geNomad on held-out genomes, stratified by size | `docs/model-quality-baseline.md` | "agrees with geNomad on 97% of elements in a sample held out from training but used for tuning (covered50, CI)", replaced by the lockbox figure once scored, **never** "96% accurate" |
| Robustness facts measured so far | RC exact; grid phase 1.7%; short-contig scan recovers covered 0.939 on shredded input; ambiguity `impute` | plain statements, with the numbers |

### Rung 2: improvement claims that need one new run

| claim | run needed | wording by outcome |
|---|---|---|
| **B1. Recall on independent curated truth** | prophound (ensemble, all three modes) on LAMBDA's clean subset, and on Philympics' clean subset | Above 0.818 covered50: "highest recall of the N tools compared, at X calls/genome versus Y". Within the CI of the leaders: "comparable to annotation-based tools without gene calling or a database". Below: report it plainly, since it is still known accuracy. |
| **B1'. Against sequence-only models at matched cost** | same run. Rank calls by score and keep the top *k* per genome, to match the genomic LMs' ~3.5 calls/genome | if higher: "exceeds published sequence-only (genomic language model) approaches at matched follow-up cost". **Matched cost is required.** Comparing Exploratory at 29 calls/genome with models at 3.5 would be an unfair comparison. |
| **B2. Complement to geNomad** | same run, geNomad ∪ prophound | above +8 elements: "adds more curated prophages to geNomad's calls than any other single tool tested, at +Z calls/genome". Otherwise report the gain against the table above. **This is the direct test of the concept.** |
| **B3. Fragmented input** | shred the clean LAMBDA genomes to a MAG-like contig-length distribution. Run prophound, Jaeger `-p` (with `--lc` lowered, and stated) and geNomad on the same shreds | "retains X% of its complete-genome recall on fragmented assemblies". A comparative claim is allowed only with the others run on identical input. It is a good candidate: DEPhT degrades on assemblies with more than 100 contigs, and Jaeger scans only contigs over 500 kb by default. |
| **B4. Robustness to ambiguity codes** | `821396` (in flight) | self-referential only ("recall within X at Y ppm"), unless other tools are tested on the same spiked input |
| **B5. Footprint** | measure the bundled weights (MB) and CPU time per genome. Measure the geNomad database size on ARC (`du -sh` of the DB directory). The survey did not find a documented figure. | plain statement, with both numbers |

### Rung 3: supports the exploratory framing, but can't be the headline

| evidence | source |
|---|---|
| Off-target candidates are enriched for phage features | `821371` (item 1a) |
| Where our off-truth calls land among LAMBDA's 292 classified candidates (LIKELY / POSSIBLE / WEAK phage, GI/ICE, UNLIKELY), clean subset | same run as B1, no extra compute |
| Recovery of geNomad-negative elements | `821410` in-situ arm, `821391` |
| Retrospective discovery | item 2d, not started |

### Not claimable in v1

Validated novel elements; mechanism (1b); protein-level novelty; "first"; "guarantee"; any accuracy figure that is really agreement with geNomad.

## What makes exploratory outputs trustworthy

The user's second point: the tool itself, not only its claims. Each item below is either verified today or points to a measurement.

### Usability, checked against the ten recommendations for bioinformatics CLIs (Seemann 2013)

| recommendation | status (checked 2026-09-28) |
|---|---|
| print usage with no arguments; `-h` | yes |
| **`--version`** | **missing**: `prophound --version` is rejected. It is needed for reproducibility and should also be written to `run_metadata.json`. |
| **messages to stderr, not stdout** | **no**: progress and warnings use `print()` to stdout |
| non-zero exit code on failure | argparse errors exit 2. Runtime failures need checking. |
| **validate parameters and input** | partial. Missing input, empty file, gzip, protein FASTA and duplicate contig IDs are untested and should fail with clear messages. |
| no hard-coded paths; bundled data found relative to the install | **open**: `--config`/`--checkpoint` are still required (checklist §5) |
| check dependencies | n/a beyond the Python environment |

Also:
- **38 flags on `predict`.** Group them into basic and advanced in `--help`, and put the experimental and not-recommended ones (`--num-tiles`, the edge-ramp and split parameters) under advanced.
- **Standard output formats.** Add BED and GFF3 for genome browsers, and optionally a FASTA of called regions for follow-up work.

### Per-call information, so a user can triage candidates

Today `regions.tsv` gives coordinates, `mean_prob`, `max_prob`, `boundary_kind`, `refine_status` and `n_pieces`. Checklist §4 lists evidence for triage as "flagged for later". Under this principle it is what makes exploratory calls usable. The minimum set:
- **Rank** within the genome.
- **Seed support** (3/3, 2/3, 1/3). Single-seed calls are 99% off-target (`docs/error-analysis.md`), so this is the most informative single column already available.
- **Flags**: found by a short-contig scan, contains ambiguity codes, touches a contig edge, merged from pieces or possibly conflated.
- **Size-class caveat**: element size classes with lower measured recall (sub-10 kb) are marked, with the recall figure in the docs.
- **Optional, calibrated**: an empirical "P(overlaps a known element)" per score/support bin, from data held out from tuning (the lockbox). This must be labelled as agreement with known elements, not as the probability the call is real.

### An optional known-element overlay (proposal, user decision)

A `--known <regions.bed|genomad_provirus.tsv>` flag that labels each call **known** (overlaps the supplied calls) or **candidate**. The tool stays database-free, because the user supplies another tool's output. It turns "runs alongside annotation tools" from a sentence in the README into a feature, and hands the user exactly the candidate list the concept promises.

### Run report

A short `summary.tsv` (or a block in `run_metadata.json`) with mode, number of calls, the candidate count if `--known` is used, contigs and bp skipped or short-scanned, ambiguity bases imputed, and the version and checkpoint checksums.

### Determinism

A golden `regions.tsv` on a fixed genome (checklist §6), plus a CPU / MPS / CUDA agreement check, so the published accuracy numbers can be reproduced on the user's hardware.

## Sources

- LAMBDA: [paper](https://academic.oup.com/nargab/article/8/3/lqag103/8787816), [data and visualizer files](https://github.com/leannmlindsey/lambda-benchmark) (per-genome JSONs in `public/data/`), [code](https://github.com/leannmlindsey/LAMBDA)
- Philympics: [paper](https://f1000research.com/articles/10-758), [framework and gold standard](https://github.com/linsalrob/ProphagePredictionComparisons) (`genbank/prophageCoords.tsv`; new tools are added as Snakemake rules, MIT licence)
- geNomad presets: [post-classification filtering](https://portal.nersc.gov/genomad/post_classification_filtering.html)
- Metagenomic phage-tool benchmark (contig-length effects; prophage tools not evaluated on fragments): [Frontiers 2023](https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2023.1078760/full)
- CLI usability: [Seemann 2013, *GigaScience*](https://pmc.ncbi.nlm.nih.gov/articles/PMC4076505/); [Bionitio](https://academic.oup.com/gigascience/article/8/9/giz109/5572530)
