# Claims and evidence

> **Sample roles** (`docs/research-integrity.md` C1, added 2026-09-29). In this doc, "held out" means held out from *training*. No sample here is held out from *tuning*. The 42-genome / 179-element sample is the **tuning sample**: thresholds, merging, the seed vote and the modes were chosen on it. The 18-genome / 52-element sample is **confirmation after use**: it chose the coverage builder, surfaced the ambiguity bug, and was in view for the 2026-09-29 mode decisions. Only the lockbox (`docs/lockbox.md`) is held out from both. Numbers below are as measured at the time and may predate the final modes (`docs/cli-and-packaging.md` "Seed vote per mode").

> **Commit hashes** (added 2026-10-27). Hashes cited in this doc before 2026-10-27 refer to the history before a rewrite that kept internal files off GitHub. The code is unchanged. `382efba` is now `26880c3` and `ad5b28b` (`ad5b28ba…`) is now `09c7756` (`docs/lockbox2.md` "Amendment 2026-10-27").


Living doc, updated in place. Built 2026-09-28. It says which novelty and improvement claims the 0.1.0 write-up can make, what evidence each needs, and what the tool itself must do so that its exploratory outputs can be trusted. The positioning argument and related tools are in `docs/related-tools.md`. Release tasks are tracked in an internal checklist.

## The principle (user, 2026-09-28)

The tool should make a novelty or improvement claim confidently without overselling. It should **not rest on the claim alone.** It must be technically robust and easy to use, and its accuracy on *known* elements must be measured and published, so that its exploratory outputs can be considered with confidence. The order that follows from this: known accuracy first, then robustness and usability, then the improvement claims, and the concept last. The concept is only as credible as the measured layer under it.

## Release gate: new, or measurably better (user, 2026-09-28)

> "If this isn't doing anything new, it at least needs to be better."

The concept (rung 1 below) is new only as a *contract*, and the methods are not new. So 0.1.0 has to show it is **better than existing tools on at least one measured axis**, on independent truth or identical input, not only on agreement with geNomad. **Scope rule (user, 2026-09-28): an axis is in the gate only if measuring it will not cause major delay.** The axes are sorted by what they cost to measure.

### In the gate: all measured by one run, no new tools to install

One ensemble run of prophound (all three modes) on LAMBDA's 39 clean genomes, scored locally against the published predictions of the other tools. Prerequisites: a 21-mer containment check against training positives (reuse `small_element_spikein.py`'s `disjoint`) and confirming LAMBDA's coordinate frame. The cost is about one GPU job the size of the 42-genome runs (1–2 h). Scoring against the other tools' predictions is plain Python and **can run on ARC** in the same job or a CPU follow-up (user, 2026-09-28), so it adds no local turnaround. The LAMBDA JSONs are about 20 MB and are pushed alongside.

| # | axis | bar to clear (LAMBDA clean subset) | outlook |
|---|---|---|---|
| G1 | **Complement to geNomad** | geNomad ∪ prophound gains **more than +8** curated elements over geNomad alone (0.809, 123/152). +8 is the best existing pairing (PHASTER), at +3.9 calls/genome. | **Plausible**, most likely in Exploratory. It is the direct test of the concept. Report the extra calls/genome with it. |
| G2 | **Sequence-only, at matched cost** | top-*k* calls per genome (k matched to ~3.5) reach covered50 **above 0.727** (best genomic LM, NTv2; Evo2 0.724) | **Plausible.** It is trained on integrated elements in host context, whereas the LMs trained on free phage genomes. |
| G3 | **Overall curated recall** | covered50 **above 0.818** (PHASTER; PIDE 0.812, geNomad 0.809), with calls/genome | **Unlikely at Default**, because the model inherits geNomad's misses. Possible in Exploratory, at about 7× the calls/genome, and must be stated that way. |
| G4 | **Footprint** (run-time data: the bundled weights-only model is 16.5 MB, against geNomad's DB at 1.47 GB unpacked, ~89×; the training checkpoints were 49.5 MB with optimizer state; see the release text, "Footprint") | install size versus geNomad's database: `genomad_db_v1.9.tar.gz` is **842 MB compressed** (Zenodo 14886553; measure the unpacked size). Our side is 3 checkpoints plus the package. | **Likely.** It counts only alongside a G1–G3 recall figure, since a small tool that finds less is not better. |

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

prophound leads when an element survives in one or two pieces. It trails once the element is in three or more, and on 60 kb+ elements (0.767 against 0.800 and 0.800). The comparators classify whole contigs: a fragment lying entirely inside an element is reported as a virus at full length. That is why geNomad's calls rise from 4.2 to 14.4 per genome and its retention exceeds 1. prophound has no contig-level call and never merges across contigs. Skipping contigs under 3 kb is not the cause: the 22 elements only prophound misses have a median of 3% of their bp in such pieces, against a median of 0.47 of the element in their largest piece. A contig-level "this whole fragment is phage" output is a plausible post-0.1.0 fix, untested.

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

### Amendment 2026-10-05: development result (not the freeze)

`824148` was clean: 300 of 300 runs (60 genomes × 5 conditions), no failures. `scripts/g5dev_eval.py` ran unchanged, using the selection rule from its docstring. Output: `analysis/g5dev/sweep.json`. No LAMBDA or Philympics output was read.

**No variant changes complete-genome recall** (Default 0.972 tune / 1.000 check for every variant), so every retention gain below is a gain in MAG recall, not an artifact of the ratio.

| mode | variant | set | MAG covered50 | retention | draft | MAG calls/genome | MAG false/Mb (cap) |
|---|---|---|---|---|---|---|---|
| Default | shipped | tune (42) | 0.808 | 0.831 | 0.955 | 17.2 | 2.22 |
| Default | ES 10 kb only | tune | 0.832 | 0.856 | 0.961 | 17.2 | 2.22 |
| Default | **ES 10 kb + WC(0.8, 0.7)**, selected | tune | **0.840** | **0.864** | 0.966 | 18.2 | 2.39 (2.46) |
| Default | shipped | check (18) | 0.878 | 0.878 | 0.981 | 15.1 | 1.93 |
| Default | ES 10 kb only | check | 0.904 | 0.904 | 0.981 | 15.1 | 1.92 |
| Default | selected | check | 0.904 | 0.904 | 0.981 | 15.9 | 2.06 |
| Exploratory | selected: ES 10 kb + WC(0.8, 0.7) | tune / check | 0.937 / 0.949 | 0.942 / 0.949 | — | 41.3 / 40.4 | 7.51 (7.63) |
| Conservative | selected: ES 10 kb only (every WC variant exceeds the cap) | tune / check | 0.732 / 0.808 | 0.775 / 0.840 | — | 10.4 / 8.7 | 0.91 (1.10) |

Reading:
- **Edge snapping does most of the work, and costs nothing** (calls and false/Mb unchanged). The whole-contig call adds +0.008 MAG recall on the tune set, about 1 call/genome, and **nothing on the check set**. Its Default false/Mb sits just under the cap (2.39 vs 2.46).
- **ES = 10 kb is the edge of the grid** (0 / 2 / 5 / 10 kb), so the optimum may lie beyond it. Extending the grid now would mean choosing after seeing results, and needs its own amendment if done.
- **Phage-free regression:** on the 20 geNomad-negative genomes the selected variants leave calls/genome unchanged in all three modes (Default 7.40, Exploratory 27.15, Conservative 2.45). WC fires on 2 contigs there and replaces the calls already on them.
- **Per-mode selection:** the rule was applied per mode, so Conservative gets ES only. The preregistered endpoint is Default only. Whether the other modes get the same fix or their own is a user decision taken at the freeze.
- **Expected G5b outcome:** a retention gain of +0.033 (tune) / +0.026 (check) covers about a third of the 0.07–0.11 deficit to the comparators seen on LAMBDA. Clearing G5b needs a lead of about 0.13, so **inconclusive or not clearing remains the expected result**, as the power note above says.

### Amendment 2026-10-05 (later): ES vs ES+WC, one fix for all modes

**User decisions, same day:** the fix applies to **all modes** (no per-mode variants). It is **on by default if the gain is meaningful.** ES and ES+WC were compared before choosing. Script: `scripts/g5dev_compare.py` (same post-processing code as the sweep) → `analysis/g5dev/compare.json`.

**Under one fix for all modes, the rule has to hold in every mode.** Applied jointly, the docstring rule (best Default MAG retention, with each mode's MAG false/Mb at or below its own complete-genome false/Mb) picks **ES 10 kb with no whole-contig call**. ES+WC(0.8, 0.7) is out because it raises Conservative's MAG false/Mb to 1.33, above its cap of 1.10. WC variants at the mode's own threshold tie with ES alone and lose the tie-break on simplicity.

MAG covered50, three shred seeds pooled. Paired genome bootstrap, 95% interval:

| set | mode | shipped | ES | ES+WC | ES − shipped | (ES+WC) − ES |
|---|---|---|---|---|---|---|
| tune (42) | Default | 0.808 | 0.832 | 0.840 | +0.024 (+0.012, +0.038) | +0.007 (+0.002, +0.015) |
| tune | Exploratory | 0.916 | 0.931 | 0.937 | +0.015 (+0.004, +0.027) | +0.006 (0.000, +0.012) |
| tune | Conservative | 0.709 | 0.732 | 0.778 | +0.022 (+0.011, +0.036) | +0.047 (+0.027, +0.070), **over the cap** |
| check (18) | Default | 0.878 | 0.904 | 0.904 | +0.026 (+0.005, +0.057) | 0.000 |
| check | Exploratory | 0.929 | 0.936 | 0.949 | +0.006 (0.000, +0.024) | +0.013 (0.000, +0.037) |
| check | Conservative | 0.795 | 0.808 | 0.821 | +0.013 (0.000, +0.027) | +0.013 (0.000, +0.033) |

- **ES is free.** Complete-genome recall, calls/genome and false/Mb are unchanged in every mode and on both sets. The gain sits entirely in elements cut into 3 or more fragments (Default, tune: 0.765 → 0.806), the split G5's post-hoc analysis pointed to.
- **In Default, WC recovers 4 element-seeds on the tune set (about 1.3 elements) and none on the check set.** It costs about 1 call/genome on MAG input. About a third of WC calls touch no element: 383 of 1,057 in Default, and more than half in Exploratory.
- WC never loses an element that ES finds.

**Trade-offs ES has but the selection metrics cannot see** (user asked, 2026-10-05). Covered50 and false/Mb count elements and calls, not call length. ES never adds a call, but it lengthens calls. Measured by `g5dev_compare.py` (now with bp and IoU) and `scripts/g5dev_es_distance.py`:
- **It moves most calls on MAG input.** Fragments have an N50 of 20 kb, so in Default 85–90% of MAG calls lie within 10 kb of a fragment end. About half the moved calls touch no element. The median extension is under 1 kb, but the tail is long.
- **The extra sequence is mostly off-element.** Default MAG bp called rises 12.8% (tune) / 13.0% (check), and about 88% of the added bp lies outside geNomad elements. Exploratory rises about 20%, Conservative about 11%. Calls/genome does not change.
- **Boundaries get slightly worse on MAG input.** Median IoU of detected elements: Default 0.496 → 0.486 (tune), 0.466 → 0.465 (check). Complete genomes are unaffected (IoU unchanged; 19 of 596 calls move, +96 kb in total).
- **Snapped edges are not detected boundaries.** A snapped edge reports a fragment end, not where the element stops. This needs an output flag (for example an `edge_snapped` column per side). The flag also tells the user the element may continue on another contig.
- **The endpoint cannot see the cost.** G5b scores covered50, so the longer the snap distance, the better G5b looks, and the bp cost never registers. That is a reason not to pick the distance by the endpoint alone.

**Snap distance, MAG covered50 / bp called vs shipped** (Default; Exploratory and Conservative follow the same shape):

| d | tune | check |
|---|---|---|
| 0 (shipped) | 0.808 / — | 0.878 / — |
| 1 kb | 0.819 / +2.0% | 0.891 / +1.8% |
| 2 kb | 0.823 / +3.0% | 0.891 / +2.8% |
| 5 kb | 0.827 / +6.6% | 0.904 / +6.4% |
| 10 kb | 0.832 / +12.8% | 0.904 / +13.0% |

Recall rises roughly with log d while added sequence rises roughly linearly. On the check set, 5 kb equals 10 kb in every mode. **Choosing a distance below 10 kb on bp grounds would add a criterion after the rule was fixed.** It is still on corpus data only, with G5b unseen, but the freeze amendment must say so.

### Amendment 2026-10-05: freeze of the fragmented-input fix (before prophound runs on G5b)

**Commit: `ad5b28ba26bfb6d4ab9b4f12160e9d112ea02f39`** (2026-10-05). The fix is frozen at this commit, and prophound runs on G5b from it.

- **The fix: edge snap only, on by default in all three modes.** `--edge-snap-frac 0.5` (`DEFAULT_EDGE_SNAP_FRAC`) moves a call edge within 0.5 × `--window-len` of a contig end onto that end. At the trained 10 kb window that is **5 kb**. Short-contig calls use the same distance. Off with `--edge-snap-frac 0`. No whole-contig call.
- **Flag:** a `contig_edge` column (none/start/end/both) in `regions.tsv` and the GFF3. It marks every edge on a contig end, snapped or not. `run_metadata.json` `edge_snap` records `frac`, `snap_bp`, `calls_moved` and `calls_on_contig_edge`.
- **Departures from the preregistered selection rule, all user decisions on 2026-10-05, all on corpus data with G5b unseen:**
  1. One fix for all modes, so the cap applies to every mode. The rule's per-mode pick for Default and Exploratory (ES 10 kb + WC) breaks Conservative's cap. The joint rule picks ES 10 kb.
  2. **5 kb, not 10 kb, chosen on a criterion added after the rule was fixed:** sequence called. 10 kb adds about 13% more MAG bp (Default), mostly off-element. 5 kb adds 6.5% and equals 10 kb on the check set. Default MAG covered50 on the tune set is 0.827 at 5 kb against 0.832 at 10 kb.
  3. **Distance relative to `--window-len`**, so it scales with larger windows. Snapping short-contig calls by their own 5 kb or 3 kb scale was tried and lost 1 of 156 element-seeds on the check set, so they use the coarse distance. Window lengths other than 10 kb are untested, like the thresholds.
- **Verified:** the shipped `snap_contig_edges` reproduces the measured flat 5 kb exactly, in every mode, on both sets and on the phage-free genomes (`scripts/g5dev_snap_check.py` → `analysis/g5dev/snap_check.json`).

| Default, MAG covered50 | shipped | frozen fix |
|---|---|---|
| tune (42) | 0.808 | 0.827 |
| check (18) | 0.878 | 0.904 |

  Complete-genome recall, calls/genome and false/Mb are unchanged. Unit tests: `tests/test_edge_snap.py`. End-to-end smoke run (random weights): the column, the GFF3 attribute and the metadata are written.
- **Scoring G5b with the fix:** `g5_score.py` rebuilds prophound calls offline from the tracks, so it now applies the shipped snap: `--edge_snap_frac 0.5 --lengths <table>`. The table gives exact contig lengths for every G5b condition, built on ARC like `analysis/g5dev/contig_lengths.tsv`. At 0 it reproduces G5. Cross-checked: on the 18, `mag_s0`, it matches the dev pipeline (0.846 → 0.885, 269 calls).
- **Expected G5b outcome: inconclusive or not clearing.** The gain is about a quarter of the LAMBDA deficit.
- **G5b run (2026-10-05):** `826741`. `g5_prophound.slurm` is parameterized for G5b and also runs the complete genomes, since no earlier run exists. Its preflight requires the frozen snap, and the push was verified on ARC (`predict_genome.py` md5 `62a40e2c38c0a22dbb5d8fa45e82b492`). The contig-length table `analysis/g5b/contig_lengths.tsv` was built locally from inputs that match `analysis/g5b/input_md5.tsv` (all 55 files). MAG fragment counts are 3,145–3,247, as preregistered. The comparator outputs (`824150`/`824151`) are still unread.

### G5b result (`826741`, `824150`, `824151`; scored once on 2026-10-16)

**Runs were clean.** prophound 45/45 (`.err` empty, exit 0, preflight checkpoint checksums matched). Every `run_metadata.json` has `edge_snap.frac` 0.5 / `snap_bp` 5000, Default mode, `impute`, and one checkpoint-hash set. geNomad 45/45 and Jaeger 45/45 per model, both exit 0 with empty `.err`. The comparator tables collected locally match the ARC job summaries call for call. The scorer was run once as preregistered (`g5_score.py score … --edge_snap_frac 0.5 --lengths analysis/g5b/contig_lengths.tsv`). Outputs are in `analysis/g5b/results.{json,csv}`, and every element is reachable in every condition.

**Covered50 on Philympics truth, 9 genomes / 32 elements.** MAG is pooled over three seeds. Intervals are 95% paired genome bootstrap.

| tool | complete | MAG | **retention** | MAG calls/genome | draft | draft retention |
|---|---|---|---|---|---|---|
| prophound Exploratory | 1.000 | 0.969 | 0.969 | 41.1 | 1.000 | 1.000 |
| **prophound Default** | 1.000 | **0.938** | **0.938 (0.889–0.983)** | 20.2 | 1.000 | 1.000 |
| prophound Conservative | 1.000 | 0.875 | 0.875 | 13.5 | 1.000 | 1.000 |
| **geNomad** (end-to-end) | 1.000 | **0.958** | **0.958 (0.911–1.000)** | 14.9 | 0.969 | 0.969 |
| **Jaeger 57341_1.5M** | 0.875 | **0.771** | **0.881 (0.750–1.000)** | 14.9 | 0.906 | 1.036 |
| geNomad find-proviruses only | 0.969 | 0.323 | 0.333 | 2.5 | 0.781 | 0.806 |
| Jaeger `default` model | 0.375 | 0.229 | 0.611 | 7.4 | 0.094 | 0.250 |

**Verdict: G5b does not clear.** Both preregistered conditions fail against geNomad.
- Retention, prophound Default minus each comparator: geNomad **−0.021** (−0.074 to +0.026), Jaeger **+0.057** (−0.064 to +0.165). Neither interval lies above 0.
- Absolute MAG covered50: geNomad **−0.021** (−0.074 to +0.026), so prophound is not ahead of geNomad. Jaeger **+0.167** (+0.070 to +0.256).
- Matched cost (prophound cut to the comparator's per-genome call count), MAG s0/s1/s2: 0.875 / 0.969 / 0.906 against geNomad's 0.969 / 1.000 / 0.906, and 0.844 / 1.000 / 0.844 against Jaeger's 0.781 / 0.813 / 0.719.
- This is the expected outcome stated in advance (inconclusive or not clearing). Unlike G5, it is **not a refutation**: retention against geNomad is a tie, where on LAMBDA it was −0.106 with its interval below zero. But this test was not built to show that the fix closed the gap. Both prophound and geNomad sit at the 1.000 ceiling on complete genomes, so retention reduces to absolute MAG recall. There are only 32 elements, and five of the 9 genomes are *Paracoccus*.

**What the release may say (unchanged from G5, with G5b added):**
- **Not supported:** any comparative fragmented-input claim against geNomad.
- **Supported, self-referential:** on held-out MAG-like input, prophound Default kept 0.94 of its complete-genome recall on Philympics curated prophages at 20.2 calls/genome, and kept all of it on draft-like input. Report geNomad (0.96) and Jaeger (0.88) on the same input alongside. The LAMBDA G5 figures (0.92, pre-fix) stay as reported.
- **Not claimed:** prophound's absolute lead over Jaeger's primary model (+0.167, interval above 0). Jaeger was not the axis the gate tested alone, and it misses 4 of 32 elements even on complete genomes. It is reported, not headlined.

**Post-hoc (2026-10-16, after scoring; not evidence for G5b): the fix's effect on G5b.** The same tracks were rescored with the snap off (`--edge_snap_frac 0`, giving `analysis/g5b/results_posthoc_snapoff.json`). The result is **byte-identical to the snap-on score in every mode and condition.** The snap does act: `scripts/g5b_snap_posthoc.py` shows it raises coverage on 20–26 elements per MAG seed, by at most +0.153. But none of them crosses the 50% line. Its cost matches the development estimate (+6.6% at 5 kb):

| Default | bp called, off → on | added | of which on-element |
|---|---|---|---|
| complete | 3,884 → 3,904 kb (+0.5%) | 20 kb | 2 kb (10%) |
| MAG s0 / s1 / s2 | +7.4% / +6.2% / +6.0% | 180 / 170 / 151 kb | 24 / 21 / 34 kb (12–22%) |
| draft | 3,478 → 3,556 kb (+2.3%) | 79 kb | 10 kb (13%) |

- **On fresh data the fix bought no recall at the endpoint.** Its development gain (+0.019 tune / +0.026 check) did not reproduce on 32 elements. Here an element only counts if it crosses 50% coverage, and none did.
- **So G5b's tie with geNomad is not the fix's doing.** Without the snap, prophound already scores 0.938 against geNomad's 0.958 on this set. The change from G5 (−0.106) to G5b (−0.021) comes from the genome set (Philympics: smaller, with both tools at the complete-genome ceiling), not from the fix. Never write it up as "the fix closed the gap".
- **The fix's measured cost is real:** about 6–7% more MAG sequence called, and 78–88% of the added sequence lies off-element. The `contig_edge` flag marks every snapped edge, so a user can discount them.
- The decision to ship it on by default stands as frozen. Reopening it would mean tuning after G5b, which the agreed order forbids. Whether to revisit the default post-0.1.0 is a user decision.

### Measured by the same run, where the tool is probably not better (report plainly)

| axis | existing level (LAMBDA clean) | our reference point | reading |
|---|---|---|---|
| **Boundary precision** | median IoU of detected elements: geNomad **0.921**, PIDE 0.912, VIBRANT 0.885, PHASTER 0.831. IoU ≥ 0.5 on 0.72 / 0.76 / 0.68 / 0.70 of elements | **Measured, 0.1.0 (2026-10-23):** Default median IoU 0.81 on LAMBDA (0.77 without composition edges), 0.81 on the tuning sample against geNomad; before the boundary stages about 0.69–0.74 | **Worse than the annotation tools, closer than before.** Annotation tools snap edges to genes, tRNAs and integrases. Ours come from a fine scoring pass, att direct repeats (base-exact where found, about one call in seven) and composition change-points (`docs/att-boundaries.md`). State it as a known limitation; it is a real cost of skipping gene calling. |
| **Small elements (<10 kb)** | 18 of 154 clean LAMBDA elements | 5–8 kb covered50 0.857 on geNomad truth (n=7) | Known weakness, already in the limitations. |
| **Calls/genome at the same recall** | 3.4–4.5 for most tools | Default 14.7, Conservative 7.3 (0.1.0, LAMBDA; 14.5 and 9.4 at the 2026-09-28 check) | Costlier by design. Frame it as follow-up budget, never hide it. |

### Not in the gate (would cause delay)

- **Ambiguity robustness versus other tools.** It needs spiked input run through every comparator. Our own dose-response (`821396`) is reported as a robustness fact, not a comparison.
- **Contributing prophound to the Philympics framework** (a Snakemake rule and a PR to an external repo, on their timeline). Post-v1. *Scoring* Philympics' 23 clean genomes ourselves can be added to the G1–G4 job if it costs less than a day, since the published CDS-level results are in its repo (`jupyter_notebooks/all_tptn.tsv`).
- **Anything requiring retraining** (small elements, boundary head): deferred past 0.1.0 already.

### Decision rule

- **If any of G1–G3 or G5 clears its bar**, the write-up leads with the axes that cleared: each with its bar, our number, the calls/genome and the n. G4 and G6 go alongside, never alone.
- **Every axis measured is reported**, including those where the tool is not better (boundary precision, small elements, cost per genome). Reporting the losses is what makes the wins credible.
- **If none clears**, 0.1.0 cannot claim to be new or better, and the release decision goes back to the user.

## Two layers of "known accuracy"

This is the finding that most changes how results should be written up.

| layer | what it measures | current number |
|---|---|---|
| **Agreement with geNomad** | how well the tool reproduces its label source on genomes held out from training | covered50 0.972 (0.949–0.994; ensemble, final Default, n=179). This is the **tuning sample**. **Held out, 0.1.0 pipeline (second lockbox, 2026-10-27): 0.989 (0.962–0.997 Wilson; 0.972–1.000 bootstrap; n=186, 150 genomes) at 12.4 calls/genome** (`docs/lockbox2.md` "Result"). First lockbox, pre-boundary defaults: 0.960 (0.924–0.980; n=202) at 11.0 |
| **Agreement with curated truth** | how well it finds prophages that experts annotated | **Default 0.933** of 149 curated elements on LAMBDA's clean subset (genome bootstrap 0.874–0.983; 0.920 on the 34 genomes not in the corpus), at 14.7 calls/genome, against geNomad **0.822** on the same subset ("LAMBDA result" and the 2026-10-23 amendment below) |

When this section was written (2026-09-28), geNomad was estimated to cover about 81% of expert-curated prophages, so a tool matching geNomad at 96% was expected to land near 80% on curated truth. The corrected baseline puts geNomad at 0.822. prophound came out higher, at 0.933, because it calls much more sequence: per unit of cost it is level with geNomad (release text §2). The principle stands. **Agreement with geNomad measures agreement, not accuracy. Publishing it alone would oversell the tool, so the release reports both layers.** Curated truth is the number that tells a user how much of what is really there the tool finds.

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
- **Much of the raw margin is cost.** Default calls 9.8% of each chromosome, geNomad 3.5%. Rotating each tool's calls to random positions gives a chance covered50 of 0.144 for Default and 0.050 for geNomad. **Recall above chance is 0.787 against 0.772** (Conservative 0.808, PHASTER 0.771). Per unit of cost, prophound is level with the annotation tools. It does not beat them. *Reproduced 2026-10-16 with a checked-in script (`scripts/lambda_chance.py`, 1,000 rotations, `analysis/external_benchmarks/lambda/chance.json`). It gives chance 0.141 / 0.052, geNomad's fraction called 3.7% (not 3.5%), and above chance 0.792 against 0.770. Paired on the 37 shared genomes, Default minus geNomad is +0.018 (genome bootstrap −0.037 to +0.078). The original figures had no script behind them.*
- **The gains are concentrated.** 12 of the 16 elements come from two genomes (Photorhabdus TTO1, 9; *P. aeruginosa* LESB61, 3), and only 5 genomes gain anything. Two are in the test-split EDL933.
- **The gains are not training look-alikes.** All 16 have containment ≤0.16 against training positives, 13 of them <0.01. Elements under 10 kb favour prophound, but n = 13 overall and 8 not in the corpus.
- The sample is 39 genomes, dominated by a few with many elements. Intervals overlap for every pairwise comparison except the paired one against geNomad.

**Reading for the release:**
- **Supported:** "run alongside geNomad, prophound recovered 16 additional curated prophages on LAMBDA's clean subset, at about 15 extra candidate calls per genome".
- **Supported:** it reaches the highest raw curated recall of the tools compared, *at a higher call cost*.
- **Not supported:** "more accurate than" any tool. At matched cost it ties geNomad and trails PHASTER.

**Caveat specific to prophound:** 5 of the 39 genomes are in our test split, including EDL933 (22 elements). EDL933 is also in the literature pilot, and was studied closely in the ambiguity work. The primary row is therefore the **34 genomes not in the corpus**.

### LAMBDA candidates check: preregistration (2026-10-20, before any scoring)

**Why.** Every comparative figure so far measures agreement with curated or geNomad truth. For a discovery tool, the question is whether the calls reach real elements *outside* curated truth. LAMBDA's reviewers classified the regions its language models predicted outside the curated set: LIKELY PHAGE, POSSIBLE PHAGE, WEAK PHAGE SIGNAL, GENOMIC ISLAND, GENOMIC ISLAND / ICE, UNLIKELY PHAGE. A LIKELY PHAGE candidate is a phage the reviewers believed but the curated set lacks. Requested by the user, 2026-10-20.

**Data.** The clean subset: 39 genomes (34 not in prophound's corpus) and 152 candidates (20 LIKELY, 46 POSSIBLE, 57 WEAK, 8 GI, 10 GI/ICE, 11 UNLIKELY; 127 flagged novel). The existing ensemble tracks (`analysis/lambda/ens`, the `822275` run) are rebuilt per shipped mode by `mode_recheck`/evalkit, as in the LAMBDA result. No new prophound run is needed. Script: `scripts/lambda_candidates.py` (seed 20261020; 1,000 rotations; 10,000 genome bootstraps). It was checked on synthetic data: a planted tool is far above chance and a random one is at chance.

**Primary endpoint (fixed).** prophound **Default**, the **34 not-in-corpus genomes**: the covered50 recovery of **LIKELY PHAGE** candidates **above chance**. Chance is each chromosome's calls rotated by a uniform random offset. The interval is a 95% genome bootstrap. **Positive** if the interval excludes 0.

**Secondary (reported whatever the outcome).**
- The same, for LIKELY or POSSIBLE together; for Exploratory and Conservative; for all 39 genomes; and for the novel-flagged candidates.
- **Landing:** of the calls that overlap no curated element, the fraction overlapping each class, against rotated calls (enrichment).
- **Matched cost:** prophound's calls cut per genome to geNomad's and to PHASTER's call count (by max_prob), then recovery above chance.
- **Annotation tools** (geNomad, PHASTER, PIDE, VIBRANT, PhageBoost, PhiSpy, Phigaro, ...) on the same footing, as context. Every tool is scored only on the genomes it ran on.

**Not a comparison.** The candidates were generated by LAMBDA's language models, so the LM tools are scored on regions they produced. Their rows are labelled circular and are never compared.

**Caveats fixed in advance.**
- **The sample is small.** 20 LIKELY candidates in all, fewer in the 34 genomes, and one genome can carry several.
- **The candidates are selected by LM predictions.** Phages that no LM predicted are absent, so this cannot measure what everyone missed.
- **The review classes are judgements, not validation.**
- **Higher call volume buys chance recovery.** That is why everything is reported above chance and at matched cost.

**Wording by outcome.**
- **Positive:** "prophound's calls recover X of the N phages LAMBDA's reviewers judged likely but its curated set lacks, above chance (interval), at Y calls/genome". At matched cost: say whether it holds at geNomad's call count.
- **Interval includes 0:** report the numbers plainly, without a claim.
- **Under all outcomes:** this is evidence of reach beyond curated truth, not validated discovery.

### LAMBDA candidates check: result (scored once, 2026-10-20)

Run as preregistered (`scripts/lambda_candidates.py`, 12.8 min). Outputs: `analysis/external_benchmarks/lambda/candidates.{json,csv}`. Covered50 recovery minus chance (calls rotated on each chromosome), with 95% genome-bootstrap intervals.

**Primary endpoint: positive.** On the 34 not-in-corpus genomes, prophound Default covers **18 of 19 LIKELY PHAGE candidates** (0.947), against 0.111 by chance: **0.837 above chance (0.694–0.945)**, at 14.4 calls/genome. The candidates are independent of curated truth: none of the 19 overlaps a curated element, one lies within 10 kb of one, and they are spread over 13 genomes (at most 3 each).

**The context that sets the wording: the annotation tools reach these candidates too.** On the same genomes, geNomad covers 15 of 19 (0.789; above chance 0.761, 0.543–0.927), PIDE 0.789, PhageBoost 0.778, VIBRANT 0.684, PHASTER 0.579 and PhiSpy 0.368. So LAMBDA's LIKELY candidates are mostly phages the database tools do call but its curators left out. They are not elements everyone missed.

**Paired with geNomad** (the 32 genomes both ran, LIKELY candidates):

| | covered by prophound only | by geNomad only | both | neither |
|---|---|---|---|---|
| prophound Default (14.4 calls/genome) | **3** | 0 | 15 | 1 |
| Default cut to geNomad's call count (4.4) | 3 | 2 | 13 | 1 |

The three candidates only prophound covers are **all ones LAMBDA flags as novel**.

**Secondary endpoints (not-in-corpus genomes unless stated).**

| candidates | prophound Default | Default at geNomad's cost | geNomad | PIDE | PHASTER |
|---|---|---|---|---|---|
| LIKELY or POSSIBLE (61; 57 for tools run on fewer genomes) | 0.754 (above chance 0.66, 0.51–0.78) | 0.579 (0.53, 0.39–0.66) | 0.474 (0.44, 0.28–0.58) | 0.544 | 0.263 |
| LIKELY, novel-flagged (5) | 5/5 | 4/5 *(post hoc)* | 2/5 | 3/5 | 0/5 |
| LIKELY or POSSIBLE, novel-flagged (43; 39) | 0.674 (0.57, 0.42–0.72) | 0.436 (0.38, 0.22–0.55) *(post hoc)* | 0.282 (0.24, 0.09–0.40) | 0.410 | 0.026 |
| LIKELY, all 39 genomes (20) | 0.95 (0.84, 0.70–0.94) | — | 0.75 | 0.75 | 0.55 |
| WEAK PHAGE SIGNAL (52) | 0.500 | 0.137 | 0.314 | 0.115 | 0.058 |
| UNLIKELY PHAGE (11) | 0.636 | 0.182 | 0.364 | 0.273 | 0.364 |

- Exploratory and Conservative behave the same way. On LIKELY (not in corpus): Exploratory 19/19 (above chance 0.84) and Conservative 17/19 (0.83).
- **Landing:** of Default's 383 calls outside curated truth, 13.1% overlap a LIKELY or POSSIBLE candidate, against 2.6% by chance (5.0×). geNomad's 56 off-truth calls: 57%, 17.9×. geNomad's few extra calls are concentrated on reviewed candidates. prophound's many are enriched but mostly land elsewhere. That is not evidence of error: the candidates are only the regions LAMBDA's language models flagged.
- **It also covers more of what reviewers doubted** (UNLIKELY 0.64, against geNomad's 0.36). That is the price of its call volume: at geNomad's cost it falls to 0.18.

**Reading.**
- **Supported, by the preregistered wording:** "prophound's calls recover 18 of the 19 phages LAMBDA's reviewers judged likely but its curated set lacks, above chance (0.84, 0.69–0.94), at 14 calls/genome. Cut to geNomad's call count, 16 of 19, against geNomad's 15."
- **Supported, descriptively:** the three LIKELY phages it reaches that geNomad does not are all novel-flagged by LAMBDA. On novel-flagged LIKELY or POSSIBLE candidates it covers 0.67, against geNomad's 0.28. At geNomad's call count that is 0.44 against 0.28 *(post hoc; intervals overlap, so this is a lead, not a claim)*.
- **Not supported:** that it reaches reviewed phages the annotation tools cannot. Most LIKELY candidates are called by geNomad, PIDE or PhageBoost. Nor anything about validated discovery: the classes are review judgements, and the candidates are only what LAMBDA's language models proposed.

### Amendment 2026-10-23: LAMBDA and its candidates re-scored under the 0.1.0 pipeline

The results above were scored before the boundary stages (refinement, att, composition with the 8 kb ceiling; `docs/att-boundaries.md`) and stand as scored. Re-scored from the local refinement-on run (`arc_data/att_runs/lambda`) with every 0.1.0 stage (`--genomes_dir analysis/g5/fasta` on `lambda_score.py`, `lambda_chance.py` and `lambda_candidates.py`; the last two gained the flag). Outputs `analysis/external_benchmarks/lambda/{scores,chance,candidates}_0.1.0.*`; the earlier `scores_0.1.0.*` (all stages, before the 8 kb ceiling) were renamed `scores_0.1.0_noceiling.*`.

- **Curated recall unchanged:** Default 0.933 (14.7 calls/genome), Exploratory 0.946, Conservative 0.899; G1 is the same 16 elements (+14.9 calls/genome).
- **Chance adjustment improves:** Default calls 7.9% of each chromosome (was 9.8%), chance 0.111, above chance 0.822 against geNomad 0.770. Paired: geNomad +0.048 (−0.012 to +0.109), PHASTER +0.048 (−0.007 to +0.110), PIDE +0.044 (−0.020 to +0.113). Three of eight pairings now exclude zero (Default and Conservative against VIBRANT, Conservative against PHASTER), against one before; not claimed.
- **G2 (matched call count) worsens:** Exploratory cut to EVO2 4k's count 0.792 (bar >0.785, still cleared by one element; was 0.805), to geNomad's 0.801 (was a tie at 0.822), to PHASTER's 0.770 (was 0.797). Refinement splits calls, so a count cut keeps fewer distinct regions.
- **Candidates check:** LIKELY PHAGE (not in corpus) 16 of 19 (was 18), above chance 0.753 (0.538–0.950; was 0.837); all 19 still have a Default call on them. Cut to geNomad's count 13 of 19 (was 16; geNomad 15). Paired with geNomad, prophound-only 1 (was 3), geNomad-only 0. The two that fall under 50% are two of the three novel-flagged prophound-only candidates: an att pair trims `GCF_000196155.1` 4,729,192–4,748,359 (19.2 kb) to a 9.0 kb call (0.47 covered), and composition trims `GCA_017183655.1` 936,999–965,000 (0.61 → 0.41; refinement and att alone 0.55). Candidate coordinates come from language-model windows (most ends on a 1 kb step), so which region is right is unknown. Exploratory, without composition, covers 18 of 19 (was 19). UNLIKELY PHAGE 0.55 (was 0.64).
- The release text reports both scorings (§2, §4, §9 of `docs/release-text-0.1.0.md`).
- **User decision (2026-10-23): accept and report.** No compromise rule is fitted on LAMBDA (already used); the second lockbox is the held-out test of the boundary stages.

### Amendment 2026-10-27: Exploratory re-scored with composition edges

The 2026-10-23 files were scored before Exploratory gained composition edges (gate 0.005). Re-run with the same three commands on `382efba`'s code, `chance_0.1.0.json` and `candidates_0.1.0.{json,csv}` were replaced; the earlier versions are in git. **Only Exploratory rows changed.** Default, Conservative, every other tool and `scores_0.1.0.*` are identical (byte for byte, for the scores).
- **Curated recall and calls unchanged:** Exploratory 0.946. It now calls 12.8% of each chromosome (was 15.0%), so chance falls to 0.145 (was 0.177) and recall above chance rises to 0.801 (was 0.769).
- **Candidates check, Exploratory:** LIKELY PHAGE (not in corpus) 17 of 19 (0.895; was 18 of 19), above chance 0.785 (0.573–0.950; was 0.820). All 19 still have an Exploratory call on them. UNLIKELY PHAGE 0.727 (was 0.818). Composition trims calls inward, so some candidates whose coordinates come from 1 kb language-model windows fall under 50% coverage. As in the 2026-10-23 amendment, which region is right is unknown.
- The preregistered endpoint is Default's and is unaffected. The release text quotes no Exploratory candidates figure.

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
| Design concept: recall-first screen, candidates, modes priced in calls/genome, runs alongside annotation tools | the design; survey in `docs/related-tools.md`, re-checked 2026-10-26 | "designed as". Not "first". No part is unique on its own: discovery framing (PhageBoost; geNomad's `--relaxed`, "to maximize discovery"), sensitivity settings (Jaeger `--sensitivity`, geNomad presets) and recall-leaning tools (VirSorter, VirSorter2, DBSCAN-SWA, per Philympics) all exist. What was not found in any tool of this style is the combination: an explicit recall target measured on a preregistered held-out set, and each mode priced in follow-up cost (`docs/related-tools.md` "Re-check, 2026-10-26"). |
| No gene calling and no reference database at run time | code | plain statement. Boundary ingredients (att direct repeats, composition change-points) have prior art and are cited |
| Agreement with geNomad on held-out genomes, stratified by size | **second lockbox (scored 2026-10-27, the 0.1.0 pipeline): Default 0.989 (0.962–0.997)**, so ≥95% is established by the preregistered rule; the first lockbox read 0.960 (0.924–0.980) under the pre-boundary defaults | "agrees with geNomad on 0.989 of held-out elements (covered50, 0.962–0.997), at least 95%"; **never** "99% accurate"; pooled, Default, geNomad truth only |
| Robustness facts measured so far | RC: window scores exact, and the 0.1.0 pipeline symmetric on clean sequence up to ~100 bp of composition jitter, but not with ambiguity codes (`docs/scientific-validity.md` RC-3); grid phase 1.7%; short-contig scan recovers covered 0.939 on shredded input; ambiguity `impute` | plain statements, with the numbers, including the RC-3 limitation |

### Rung 2: improvement claims that need one new run

**Outcomes (all measured; 0.1.0 figures from the 2026-10-23 re-score):**
- **B1, curated recall:** 0.933 against 0.822 for geNomad, so the highest of the tools compared, at 1.7–4.1× their calls per genome. Per unit of cost, the paired intervals span zero.
- **B1′, matched cost:** 0.792 against EVO2 4k's 0.785, cleared by one element. There is no paired separation.
- **B2, complement:** +16 / −0 over geNomad, against a bar of +5, at +14.9 calls/genome.
- **B3, fragmented input:** does not clear (G5 refuted; G5b ties geNomad).
- **B4, ambiguity:** self-referential. At most 2 of 179 elements are lost up to 3,000 ppm.
- **B5, footprint:** 16.5 MB against 1.47 GB.

Release-text §2–§6 and "Footprint" carry the wording. The table below is the plan as written on 2026-09-28.

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
| Off-target candidates are enriched for mobile-element features: 3.09× over length-matched controls (2.45–4.12), mostly conjugation and integrase/transposase rather than phage structural genes | `821371` (item 1a); release text §8 (†) |
| Where our off-truth calls land among LAMBDA's classified candidates: a call on all 19 LIKELY PHAGE candidates, covering at least half of 16 under 0.1.0 (preregistered scoring 18) | "LAMBDA candidates check" below; release text §9 |
| Recovery of geNomad-negative elements: confident in-situ inoviruses, Default 7 of 10, Exploratory 10 of 10 | `821410` in-situ arm, `821391`; release text §8 |
| Off-element att calls are mostly real integration sites: 0.495 of Default att calls off every known element have a core on a tRNA/tmRNA gene, against 0.045 for decoys and 0.004 at random | `docs/att-boundaries.md` "Are the att sites real?" (2026-10-27) |
| Retrospective discovery | item 2d, not started |

### Not claimable in 0.1.0

Validated novel elements; mechanism (1b); protein-level novelty; "first"; "guarantee"; any accuracy figure that is really agreement with geNomad.

## What makes exploratory outputs trustworthy

The user's second point: the tool itself, not only its claims. Each item below is either verified today or points to a measurement.

### Usability, checked against the ten recommendations for bioinformatics CLIs (Seemann 2013)

| recommendation | status (re-checked 2026-10-20; first checked 2026-09-28) |
|---|---|
| print usage with no arguments; `-h` | yes |
| `--version` | yes (`prophound 0.1.0`), and `run_metadata.json` records it |
| messages to stderr, not stdout | yes: a run writes nothing to stdout. Progress, warnings and the summary go to stderr; `-q` keeps warnings, errors and the summary |
| non-zero exit code on failure | yes: a fixable failure prints one `error:` line and exits 1; Ctrl-C exits 130; argparse errors exit 2 (`tests/test_cli_errors.py`) |
| validate parameters and input | yes: missing, empty, unreadable, protein, reads, alignments, duplicate contig IDs and out-of-range options fail with a reason (`src/inputs.py`, `src/io_checks.py`, `tests/test_inputs.py`) |
| no hard-coded paths; bundled data found relative to the install | yes: the model ships in the package and is found next to the code; its SHA-256s are checked at every run |
| check dependencies | n/a beyond the Python environment; the wheel declares its seven runtime dependencies |

Also, both done in the 2026-10 UI pass:
- **Flags grouped.** `--help` shows basic, advanced and experimental sections.
- **Standard output formats.** BED and GFF3 always, and a FASTA of called regions with `--write-region-fasta`.

### Per-call information, so a user can triage candidates

Under this principle, per-call evidence is what makes exploratory calls usable. The minimum set, with its status as of 2026-10-27 (`docs/output-files.md`):
- **Rank** within the genome. **Done:** `rank` in `regions.tsv`, ordered by `max_prob` then `mean_prob`.
- **Seed support** (3/3, 2/3, 1/3). Single-seed calls are 99% off-target (`docs/error-analysis.md`), so this is the most informative single column. **Done:** `support`, counted on the call the vote saw, so split pieces report their parent's support (fixed 2026-10-26).
- **Boundary evidence.** **Done:** `edges` (`att site` or `estimated`), a range on every estimated end, and each end's source in `regions.detailed.tsv`.
- **Flags**: found by a short-contig scan, contains ambiguity codes, touches a contig edge, merged from pieces or possibly conflated. **Done in part:** `note` flags calls that run off a contig end, carry 10% N or more, sit on a short contig, or have an unrefined edge. `refine_status` in `regions.detailed.tsv` marks split pieces. Ambiguity codes below 10% and possible conflation are not flagged per call.
- **Size-class caveat**: element size classes with lower measured recall (sub-10 kb) are marked, with the recall figure in the docs. **Not done** per call; the docs carry it.
- **Optional, calibrated**: an empirical "P(overlaps a known element)" per score/support bin, from data held out from tuning (the lockbox). This must be labelled as agreement with known elements, not as the probability the call is real. **Not done.** The second lockbox could supply the bins after it is scored.

### An optional known-element overlay (proposal, user decision)

A `--known <regions.bed|genomad_provirus.tsv>` flag that labels each call **known** (overlaps the supplied calls) or **candidate**. The tool stays database-free, because the user supplies another tool's output. It turns "runs alongside annotation tools" from a sentence in the README into a feature, and hands the user exactly the candidate list the concept promises.

### Run report

A short `summary.tsv` (or a block in `run_metadata.json`) with mode, number of calls, the candidate count if `--known` is used, contigs and bp skipped or short-scanned, ambiguity bases imputed, and the version and checkpoint checksums. **Done, except the `--known` count** (2026-10-20): the terminal `[summary]` lines and `run_metadata.json` carry all of it, and a batch writes `batch_summary.tsv`. The `--known` overlay above is a proposal not taken up for 0.1.0.

### Determinism

A golden `regions.tsv` on a fixed genome (the release-day checks), plus a CPU / MPS / CUDA agreement check, so the published accuracy numbers can be reproduced on the user's hardware. **Status 2026-10-19 (`docs/benchmarks.md`):**
- The check found MPS giving wrong scores (S10). It is fixed.
- CPU and MPS now give identical calls on the 30 Phase 0 reference cases. Scores agree to 1.7e-5, against CUDA to 8e-6.
- Outputs are not byte-identical across devices (float order). The release wording is "the same calls" on the verified devices: Apple M4 CPU and GPU, and NVIDIA V100S.
- MPS vs recorded CUDA on the 45 G5b runs: one internal decision flip, no released call changed.
- x86 CPU and other GPUs are untested; ARC jobs are proposed.
- The golden output is the Phase 0 harness (`scripts/reference_outputs.py`).

## Sources

- LAMBDA: [paper](https://academic.oup.com/nargab/article/8/3/lqag103/8787816), [data and visualizer files](https://github.com/leannmlindsey/lambda-benchmark) (per-genome JSONs in `public/data/`), [code](https://github.com/leannmlindsey/LAMBDA)
- Philympics: [paper](https://f1000research.com/articles/10-758), [framework and gold standard](https://github.com/linsalrob/ProphagePredictionComparisons) (`genbank/prophageCoords.tsv`; new tools are added as Snakemake rules, MIT licence)
- geNomad presets: [post-classification filtering](https://portal.nersc.gov/genomad/post_classification_filtering.html)
- Metagenomic phage-tool benchmark (contig-length effects; prophage tools not evaluated on fragments): [Frontiers 2023](https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2023.1078760/full)
- CLI usability: [Seemann 2013, *GigaScience*](https://pmc.ncbi.nlm.nih.gov/articles/PMC4076505/); [Bionitio](https://academic.oup.com/gigascience/article/8/9/giz109/5572530)
