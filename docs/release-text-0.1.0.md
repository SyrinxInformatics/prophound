# prophound 0.1.0 — release text (DRAFT)

> **Draft, 2026-10-16; lockbox figures filled in 2026-10-19** from the one scoring of `841208` (`analysis/lockbox/results.json`, `docs/lockbox.md` "Result"). **2026-10-20:** §9 (LAMBDA candidates, scored once) added, the footprint set to the bundled model, and the Phase 3 GPU timings in. Every lockbox figure below comes from that scoring. Nothing else changes when it lands, except wording the decision rule ties to it (marked ⟂). **2026-10-23: rewritten for the 0.1.0 boundary stages** (refinement on, att sites, composition edges with the 8 kb ceiling, edge ranges; `docs/att-boundaries.md`). Every figure measured on a sample with local tracks was re-scored with every 0.1.0 stage. The rest carry † (measured before the boundary stages; see "How to read"). The lockbox figures are the first lockbox's, which describe the pre-boundary defaults. They are replaced by the second lockbox's when it is scored (`docs/lockbox2.md`; marked ⟂2). **2026-10-27:** Exploratory composition edges and the novelty re-check (`docs/related-tools.md` "Re-check, 2026-10-26") folded in, and the defaults frozen for the second lockbox at `24f6c44`. **2026-10-27 (later): the second lockbox is scored** (`docs/lockbox2.md` "Result"), and its figures replace the first lockbox's wherever this text gave a held-out figure for the shipped pipeline (the ⟂2 markers are resolved). The first lockbox's figures stay only as the record of the pre-boundary defaults, marked †. Rules for this text: every claim carries its number, *n*, interval, sample role and source. Wording follows `docs/claims-and-evidence.md`'s ladder. No "first", no "guarantee", "may not call" rather than "do not call", and calls are "candidates for follow-up".

---

## What prophound is

*prophound is a recall-first prophage screen. It localizes candidate integrated elements directly from nucleotide sequence, with no gene calling or reference database at run time. It is tuned to recover already-known elements and returns everything else as candidates, at a follow-up cost the user chooses. It is designed to run alongside annotation-based tools and to surface what they may not call.*

**Method in one line:** an exploratory, database-free prophage caller that learns directly from sequence, trained on geNomad-labelled genomes. It is a reverse-complement-symmetric CNN, scored over 10 kb windows by an ensemble of three training replicates. Regions are built from the score track with no gene calling. Their ends are then placed without genes too: a finer 5 kb scoring pass, then a direct-repeat (att site) search that gives base-exact ends where it finds a pair, and a host-to-element composition change-point for the rest (every mode; a looser gate in Exploratory). Every end carries its source and, where estimated, a range.

**What is distinct is the contract, not the method.** Database-free and nucleotide-level detection already exist (PhageBoost, PIDE, Jaeger, Seeker, LAMBDA's genomic language models; `docs/related-tools.md`). So do discovery framing (PhageBoost; geNomad's `--relaxed` preset, "to maximize discovery") and sensitivity settings (Jaeger's `--sensitivity`). What none of these tools was found to offer is the combination below. prophound's contract has four parts:
1. a recall target on known elements, measured and reported with intervals;
2. three modes, each priced in calls per genome (the follow-up cost), not in precision;
3. output as candidates, not annotations;
4. use alongside annotation tools.

It rests on three design facts, each checked against the nearest tools:
- **Host-context training.** It is trained on prophages in their host genomes. Jaeger and LAMBDA's models, the nearest nucleotide-level ones, train on free phage genomes, with Jaeger masking prophages out. Gene-based callers (geNomad's provirus model, PhageBoost) also learn from genes in host genomes, so the distinction is among sequence-level models.
- **Gene-free boundaries.** Boundaries come from the score track, att direct repeats and sequence composition, never from genes. Jaeger snaps them with gene prediction. Neither ingredient is new: direct-repeat att searches go back to Phage_Finder (2006) and are in PhiSpy and DEPhT, and composition change-points have placed genomic-island edges since Alien_Hunter (2006). What is distinct is a packaged prophage caller that combines them with a model's score track and calls no genes anywhere (`docs/related-tools.md` "Prior art for the boundary stages").
- **Fragments scanned by default.** It scans contigs down to 3 kb. Jaeger's prophage mode scans contigs over 500 kb by default.

## How to read every number below

- **Recall metric: covered50.** An element counts as found when calls cover at least 50% of its bp. ≥1 bp would credit grazes. Across the ≥1 / ≥25 / ≥50 / ≥75% cut-offs, every one orders the modes identically, and ≥75% costs only 1–3 points, so covered50 is not a flattering choice (`docs/research-integrity.md` C3).
- **Intervals.** 95%, by genome bootstrap (10,000 resamples), because elements within a genome are not independent. A Wilson interval is added where a bootstrap degenerates (an estimate of 1.000) or *n* is small.
- **Cost: calls per genome.** It is the follow-up (wet-lab) budget. Precision against geNomad is not reported as accuracy. In a discovery tool, a call geNomad does not make is the intended output, not the error term.
- **Sample roles** (C1). Each figure is labelled by the role of the sample it was measured on:

| sample | role | wording used |
|---|---|---|
| first lockbox: 150 genus-disjoint genomes, 202 elements | held out from training and tuning, preregistered, scored once, under the pre-boundary defaults: the record of those defaults † | "first lockbox" |
| second lockbox: 150 further genus-disjoint genomes, 186 elements (`docs/lockbox2.md`) | **held out from training, tuning and the boundary stages**; preregistered; frozen at `24f6c44`, run once (`851599`) and scored once on 2026-10-27 | "held out" |
| 42 genomes / 179 elements | tuning: thresholds, merging, vote and modes were chosen on it | "tuning sample" |
| 18 genomes / 52 elements | confirmation, but it informed earlier decisions | "a second, previously used sample" |
| 20 genomes with no geNomad provirus | confirmation; in view for the vote decision | "geNomad-negative genomes" |
| LAMBDA clean subset: 39 genomes / 149 curated elements | external curated truth; bars fixed before the run | "LAMBDA" |
| literature pilot: 7 genomes / 57 elements | external curated truth; most contained in training (median 0.723) | "literature pilot" |
| Philympics G5b: 9 genomes / 32 elements | external curated truth; preregistered | "G5b" |

- **Pipeline version.** 0.1.0 adds three boundary stages after the region builder: refinement, att sites, and composition edges. They move ends and can split a call, so they change covered50, IoU and calls per genome. Unmarked figures were scored with every 0.1.0 stage. **†** marks a figure measured before them, on runs with no fine tracks to rebuild from (ARC runs, spike-ins, shreds). On every sample re-scored, the stages left Default covered50 unchanged and raised calls per genome by 1–5%, so † figures are expected to move little, but they were not re-measured.
- **Two layers of accuracy.** Agreement with geNomad is not accuracy, because geNomad covers 0.822 of LAMBDA's curated prophages. Both layers are reported, each labelled as what it is.
- **Novelty to training** is 21-mer containment of an element in training positives. It sees nucleotide similarity only, not protein homology.

---

## Claims

### 1. Agreement with geNomad on held-out genomes

**Default covered50 0.989 (184 of 186; Wilson 0.962–0.997; genome bootstrap 0.972–1.000), n = 186 elements in 73 of 150 genomes, at 12.4 calls per genome.** This is the second lockbox: genus-disjoint from training, from both earlier lockboxes and from every external set, preregistered, run once on the frozen 0.1.0 pipeline and scored once (`docs/lockbox2.md` "Result"). It is the primary endpoint.
- **Other modes:** Exploratory 0.984 (Wilson 0.954–0.994; bootstrap 0.956–1.000) at 37.4 calls per genome, and Conservative 0.903 (0.852–0.938; 0.858–0.949) at 5.0.
- **Covered (≥1 bp):** 1.000 / 0.989 / 0.925.
- **Misses:** Default misses 2 elements in 2 genomes, a 6.1 kb and a 34.9 kb element, both with no Default call on them. Exploratory reaches both (0.85 and 0.32 covered).
- **The first lockbox** (150 other genomes, 202 elements, scored 2026-10-19 under the pre-boundary defaults †) read 0.960 (0.924–0.980) at 11.0 calls per genome. It stays the record of those defaults (`docs/lockbox.md`). The two lockboxes are different genomes, so the difference is not an effect of the boundary stages.

Novelty to training was not measured on the second lockbox; on the first, median containment was 0.003, by the same genus-disjoint construction. The second lockbox cannot speak below 10 kb (4 elements): those rows are counts, not claims.

For comparison, on the reused samples (Default, ensemble). Their covered50 is the same before and after the boundary stages (`docs/att-boundaries.md`). Calls per genome are the 0.1.0 pipeline's, slightly higher because refinement splits some calls:

| sample | covered50 | 95% interval | calls/genome, 0.1.0 (before) |
|---|---|---|---|
| tuning (42 / 179) | 0.972 | 0.949–0.994 (bootstrap †); 0.936–0.988 (Wilson) | 14.4 (14.2) |
| second sample (18 / 52) | 1.000 | 0.931–1.000 (Wilson) | 13.6 (13.3) |
| **second lockbox (150 / 186), held out** | **0.989** | 0.972–1.000 (bootstrap); 0.962–0.997 (Wilson) | 12.4 |
| reused sample minus second lockbox: the 42 / the 18 | −0.017 / +0.011 | | |

- **≥95% is established held out, for Default.** By the rule fixed before scoring (≥95% only if the lockbox's Wilson lower bound is at least 0.95), the release may say it: the lower bound is 0.962. The claim is pooled over 186 elements, Default mode, against geNomad truth; it is agreement with geNomad, not accuracy (§2). Exploratory's lower bound, 0.954, also clears 0.95, but the rule was fixed for Default, the primary endpoint. Conservative's (0.852) does not.
- **The reused samples were not optimistic for Default.** The lockbox reads 0.017 above the tuning sample and 0.011 below the 18, both well inside either interval. For Conservative the reused samples ran higher (+0.041 and +0.059), as on the first lockbox.
- **Not carried by training look-alikes** †. On the tuning sample, the 137 elements with containment <0.1 read 0.964, and all five Default misses are in that bin (C5). Covered recall on clade-novel elements (MinHash, k = 21) is flat at 0.966–0.975 from ~3.5% out to ~20% nucleotide divergence, with the same four misses at every threshold. That was measured under the 2026-09-22 builder, at ≥1 bp (`docs/scientific-validity.md` 2a).
- **By size, second lockbox, Default** (Wilson intervals; `analysis/lockbox2/results.json`). Below 10 kb n = 2 per stratum, reported as counts, not claims. Per stratum, ≥95% is an engineering target, and no stratum is large enough to establish it on its own:

| stratum | n | covered50 | Wilson 95% |
|---|---|---|---|
| 5–8 kb | 2 | 1 of 2 | — |
| 8–10 kb | 2 | 2 of 2 | — |
| 10–15 kb | 26 | 1.000 | 0.871–1.000 |
| 15–30 kb | 39 | 1.000 | 0.910–1.000 |
| 30–60 kb | 94 | 0.989 | 0.942–0.998 |
| 60 kb+ | 23 | 1.000 | 0.857–1.000 |

- **By size, tuning sample, Default** † (Wilson intervals; `analysis/integrity/recall42.json`; the pooled figure is unchanged under 0.1.0, the strata were not re-tabulated):

| stratum | n | covered50 | Wilson 95% |
|---|---|---|---|
| 5–8 kb | 7 | 1.000 | 0.646–1.000 |
| 8–10 kb | 11 | 1.000 | 0.741–1.000 |
| 10–15 kb | 30 | 0.933 | 0.787–0.982 |
| 15–30 kb | 41 | 0.951 | 0.839–0.987 |
| 30–60 kb | 69 | 0.986 | 0.922–0.997 |
| 60 kb+ | 21 | 1.000 | 0.845–1.000 |

On the first lockbox, 15–30 kb was the weakest stratum (0.918, 4 of its 8 misses). On the second it is complete (39 of 39). Conservative is weakest there (0.821).

Wording: *"on 150 held-out genomes, Default agrees with geNomad on 0.989 of known elements (covered50, 95% interval 0.962–0.997), so at least 95% held out"*. Never "X% accurate", and never "guarantees".

### 2. Recall on expert-curated truth (LAMBDA)

On LAMBDA's clean subset (no genome in train or val), Default covers **0.933** of 149 curated elements (genome bootstrap 0.874–0.983), at 14.7 calls/genome. The primary row, 34 genomes not in the corpus, reads **0.920** (113 elements; 0.859–0.983). For the published tools on the same subset, rescored with the same rule: geNomad 0.822, PHASTER 0.824, PIDE 0.824, best genomic LM (EVO2 4k) 0.785. The bars were fixed before the run (`docs/claims-and-evidence.md` "LAMBDA result"). Re-scored with every 0.1.0 stage from a refinement-on run of the same genomes (`analysis/external_benchmarks/lambda/scores_0.1.0.json`). Covered50 is unchanged in every mode (Exploratory 0.946, Conservative 0.899), and Default's median IoU is 0.81 (0.77 without composition edges). The 8 kb composition ceiling was chosen with LAMBDA in view (`docs/att-boundaries.md`), so for that setting LAMBDA is confirmation after use.

**It is costlier, and the cost explains much of the margin.** Default calls 7.9% of each chromosome (9.8% before the boundary stages), geNomad 3.7%. Rotating each tool's calls by a random offset per chromosome keeps their lengths, count and spacing. It gives chance covered50 of 0.111 for Default and 0.052 for geNomad (1,000 draws). Recall above chance is therefore **0.822 for Default against 0.770 for geNomad**, 0.772 for PHASTER and 0.776 for PIDE (0.792 for Default before the boundary stages: tighter ends cover less sequence at the same recall).

**Paired on the genomes both tools ran, Default minus geNomad is +0.048 (genome bootstrap −0.012 to +0.109).** Against PHASTER it is +0.048 (−0.007 to +0.110), and against PIDE +0.044 (−0.020 to +0.113). Each interval spans zero, so per unit of cost prophound is not shown to be ahead of these three. Of 8 pairings (Default and Conservative against four tools), three intervals exclude zero, all small and two of them against VIBRANT: Default minus VIBRANT +0.079 (+0.009 to +0.145), Conservative minus VIBRANT +0.085 (+0.016 to +0.146) and Conservative minus PHASTER +0.055 (+0.004 to +0.109). Before the boundary stages one did (Conservative minus VIBRANT). With eight comparisons and no correction this is not claimed. Script: `scripts/lambda_chance.py --genomes_dir` (`chance_0.1.0.json`; the pre-boundary `chance.json` is kept).

- **Supported:** the highest raw curated recall of the tools compared, *at 1.7–4.1× their calls per genome* (geNomad 4.6, PHASTER 3.9, PIDE 8.4, VIBRANT 3.6), and recall above chance higher than theirs as a point estimate (0.822 against 0.770–0.776), with paired intervals that span zero.
- **Not supported:** "more accurate than" any tool.

### 3. Run alongside geNomad, it adds curated prophages (G1)

geNomad plus prophound Default recovers **16 more curated elements than geNomad alone, and loses none** (+16 / −0). On the 37 genomes both were run on, that is +0.110 covered50 (paired genome bootstrap 0.018–0.210). The cost is +14.9 calls/genome. The preregistered bar was more than +5, the best existing pairing with geNomad. Under the 0.1.0 pipeline the gain is the same 16 elements as before the boundary stages (+14.6 calls/genome then).

- **The gain is concentrated.** 12 of the 16 come from two genomes (Photorhabdus TTO1, 9; *P. aeruginosa* LESB61, 3), and only 5 genomes gain anything. The interval excludes zero, but the effect is driven by few genomes.
- **The gains are not training look-alikes.** All 16 have containment ≤0.16, and 13 of them are below 0.01.

Wording: *"run alongside geNomad, prophound recovered 16 additional curated prophages on LAMBDA's clean subset (+16 / −0), at about 15 extra candidate calls per genome"*.

### 4. Against sequence-only models at matched cost (G2)

Exploratory's calls, ranked by score and cut to EVO2 4k's per-genome call count, reach **0.792 against EVO2 4k's 0.785**, so the preregistered bar (>0.785) is still cleared, by one element. It is reported with its limits:
- at geNomad's call count it trails geNomad (0.801 vs 0.822);
- at PHASTER's lower count it trails further (0.770 vs 0.824);
- no paired interval separates it from EVO2 4k.

**The boundary stages cost this comparison.** Before them the figures were 0.805 (EVO2 4k), a tie with geNomad (0.822) and 0.797 (PHASTER). Refinement splits some calls in two, and a cut by call count then keeps fewer distinct regions. Recall and bp-based cost (§2) improved; only the count-matched view got worse.

Wording: *"at matched follow-up cost, comparable to the best published genomic language model on LAMBDA, and below the annotation tools"*. Not "exceeds".

### 5. Literature-curated pilot, split by novelty and class

There are 57 elements in 7 genomes, curated from the literature, with geNomad scored identically (`docs/final-benchmark.md`, C4 final, C5). **The pooled figure is flattered**, because 34 of the 57 elements are ≥50% contained in training. It is never quoted without the split:

| subset | n | Default | geNomad | note |
|---|---|---|---|---|
| all | 57 | 0.912 | 0.807 | paired +6 / −0, +0.105 (0.029–0.151) |
| class `prophage` (headline) | 28 | 28/28 (Wilson 0.879–1.000) | 27/28 | both at ceiling |
| novel to training (<0.1) | 13 | 8/13 (0.355–0.823) | 6/13 (0.232–0.709) | intervals overlap |
| remnants | 12 | 7/12 | 3/12 | where the difference lives |

- The six gains are Sp16, CP-933L, SpLE1, SpLE3, SpLE5 and *B. subtilis* prophage 6, and five of them are remnants or cryptic prophages.
- Every figure in the table is the same under the 0.1.0 pipeline as before it (`analysis/integrity/lit_sensitivity_0.1.0.json`; Exploratory and Conservative each lose one element, §"Operating modes").
- A permutation test against randomly placed calls is significant on both splits (Default under 0.1.0: p ≤ 1.4 × 10⁻⁵ on covered50 in each; `analysis/litpilot/results_0.1.0.json`).
- **Boundaries are closer to geNomad's but still worse:** Default median IoU 0.73 / 0.62 (test / val), against geNomad 0.89 / 0.71; before the boundary stages, 0.50 / 0.48. Where an att pair was found, ends sat a median 22 bp from the authors' coordinates (79% within 100 bp; most differences are exactly one core copy, a convention).

Wording: *"on literature-named prophages both tools are near ceiling; prophound additionally covers more remnant and cryptic prophage regions"*. This is an example, not a rate (7 genomes, two of them *E. coli* O157:H7).

### 6. Fragmented input: a design fact and a self-referential figure, no comparison

Every figure in this section is † (shred runs on ARC, with no fine tracks). The boundary stages skip contig-edge and short-contig calls, so they touch fewer calls on fragmented input than on complete genomes.

- **Design fact:** contigs of 3–10 kb are scanned by default with calibrated thresholds. Contigs under 3 kb are never scanned.
- **G5 (LAMBDA shreds, MAG-like N50 20 kb):** Default kept **0.916** (0.844–0.968) of its complete-genome recall (absolute 0.855, at 17.6 calls/genome). geNomad kept 1.022 and Jaeger 0.987. prophound's retention is lower, with the paired intervals below zero. **The comparative claim is refuted on LAMBDA.**
- **G5b (Philympics, held out, preregistered):** Default kept **0.938** (0.889–0.983), against geNomad 0.958 (difference −0.021, −0.074 to +0.026) and Jaeger 0.881 (+0.057, −0.064 to +0.165). It does not clear. On draft-like input (N50 100 kb), Default kept 1.000.
- **Fully shredded input** (tuning sample, every contig under 10 kb): Default covered50 0.626 at 20.1 calls/genome, Exploratory 0.827.
- **The shipped edge snap** extends calls within 5 kb of a contig end onto the end, and flags them `contig_edge`. On G5b it changed no element's covered50 status and added 6–7% more called sequence on MAG input, 78–88% of it off-element. It ships as frozen. Edges it moves are flagged, so users can discount them.

Wording: *"on MAG-like fragmented input, prophound Default kept 0.92–0.94 of its complete-genome recall (LAMBDA; Philympics), against geNomad 0.96–1.02 and Jaeger 0.88–0.99 on the same input"*. No comparative claim.

### 7. Robustness (self-referential)

Every row is † except the composition baseline, which is about the model, not the calls, and the reverse-complement row, which was re-run under 0.1.0.

| property | result | n / source |
|---|---|---|
| reverse complement, 0.1.0 pipeline | on sequence without ambiguity codes: window scores agree to 2e-7, coarse calls and all 15 att calls identical, composition ends differ by a median 81 bp (90% within 158 bp), and 1 of 245 refinement splits differs. **With ambiguity codes, imputation is not strand-symmetric:** on EDL933 (1,182 ppm) 4 of 28 Default calls differ, one by 8 kb | 6 genomes × 2 modes; `scripts/rc_boundary_check.py`, `docs/scientific-validity.md` RC-3 |
| ambiguity codes (IUPAC, incl. N), `impute` default | ≤2 of 179 elements lost up to 3,000 ppm scattered codes; the old encoding lost 172 | `821396`, tuning sample |
| long internal N gaps | a 20 kb gap costs 6 of 179 elements at Default. Calls ≥50% N are dropped, removing exactly the spurious gap calls with no recall change at ≤5 kb gaps | `822054`, `docs/release-scope-0.1.0.md` A.3 |
| grid phase | 1.7% of elements phase-unstable; ~5% call-count spread | `816447` |
| composition baseline | k ≤ 4 composition reaches AUPRC ≤0.16 against the CNN's 0.60–0.61 on the same split (random 0.014) | `816086`, item 1c |

The composition baseline shows the model uses more than nucleotide composition. It does not show *what* the model uses.

Earlier versions collapsed to zero on IUPAC codes (*E. coli* EDL933 at 1,182 ppm). That was fixed on 2026-09-24, and results differ from earlier runs only on sequence containing ambiguity codes.

### 8. What the off-target candidates are (exploratory evidence, not discovery)

Every figure in this section is † (ARC runs). The boundary stages trim calls (8–12% less bp called), which could change which hallmark genes a call contains.

**Candidates are enriched for mobile-element features, mostly not phage structural genes.**
- On the tuning sample, Default off-target calls made by all three checkpoints carry an MGE hallmark at **52.1% (99/190), against 16.8% (96/570) for length-matched controls** drawn clear of geNomad elements. That is **3.09× (genome bootstrap 2.45–4.12; 38 genomes)**, and Fisher one-sided p = 2 × 10⁻²⁰. Fisher treats the controls as independent, so the bootstrap is the interval to quote. These calls were made under the interim pre-vote Default, and the unanimous group is the one the shipped vote keeps.
- Calls only one or two checkpoints make are enriched 1.92× (1.46–2.78).
- **The enrichment is conjugation (16%) and integrase/transposase (13%) more than phage structural genes (4%).** Many candidates look like ICEs and other integrative elements, not missed prophages.

**Recovery beyond the labels: real inoviruses that geNomad labels negative.**
- Default covers **7 of 10** confident in-situ loci (Wilson 0.397–0.892); Exploratory covers 10 of 10 (0.722–1.000).
- These loci are not near-duplicates of training: median containment 0.001.
- One phage family, n = 10.

Wording: *"off-target candidates are enriched for mobile-element features (3.1× over background)"*. Not "phage", and never "novel elements".

### 8b. What the extra calls are made of (second lockbox, preregistered)

Every Default call on the second lockbox that touches no geNomad element (1,679 calls in 139 genomes) was sorted, by rules fixed before the data was pulled, into one bin. The same rules were applied to random host stretches of the same length (`docs/extra-calls.md`).

| bin | extra Default calls | random host stretches |
|---|---|---|
| phage-like (≥2 phage structural genes; homology v2.1) | 0.036 (0.022–0.054) | 0.005 |
| other mobile element (conjugation, plasmid or virus hallmark, integrase family, att site) | 0.243 (0.216–0.268) | 0.075 |
| island of unknown function (atypical composition only) | 0.122 (0.102–0.142) | 0.036 |
| **no detectable signal** | **0.599** (0.565–0.633) | 0.884 |

- **Most extra calls have no detectable signal.** This is the cost of the tool's recall, and users should expect it. The 0.60 is an upper bound on uninformative calls: annotation-based rules cannot tell a truly novel element from ordinary host sequence.
- **The rest are not chance.** Mobile-element and phage signals are 3–7× as frequent as in host sequence of the same length, and the phage signal survives composition-matched controls (+0.024, 0.014–0.036; `docs/candidate-phage-homology.md`).
- **Length is the best triage cue.** Calls under 10 kb show no signal 71% of the time, and calls of 60 kb and more 37% of the time (with 44% mobile-element and 13% phage-like signal).
- **Exploratory's additional calls are close to background** (0.82 no signal, against 0.88).
- **One recurring false-positive class:** ribosomal-protein operons, whose high expression gives them atypical composition. They are about 1% of the no-signal calls.

Wording: *"on held-out genomes, about 4 in 10 Default calls outside known elements carry a mobile-element or phage signal (against about 1 in 10 for random host sequence); the rest show none detectable"*. Not "the extra calls are novel elements".

### 9. LAMBDA's reviewed candidates outside its curated set (preregistered)

LAMBDA's authors list regions their language models proposed, which reviewers graded but left out of the curated truth (`docs/claims-and-evidence.md` "LAMBDA candidates check", scored once on 2026-10-20, before the boundary stages). The preregistered result stands as scored. The 0.1.0 re-score is reported beside it (`scripts/lambda_candidates.py --genomes_dir`, `analysis/external_benchmarks/lambda/candidates_0.1.0.json`).

| LIKELY PHAGE candidates, 34 genomes not in the corpus (n = 19) | preregistered scoring † | 0.1.0 pipeline |
|---|---|---|
| Default covered50 | **18 of 19** | **16 of 19** |
| chance (calls rotated on each chromosome) | 0.111 | 0.089 |
| above chance (genome bootstrap 95%) | **0.837 (0.694–0.945)** | **0.753 (0.538–0.950)** |
| calls per genome | 14.4 | 14.7 |
| candidates with any Default call on them | 19 of 19 | 19 of 19 |
| cut to geNomad's call count | 16 of 19 | 13 of 19 |
| paired with geNomad (32 genomes): prophound only / geNomad only | 3 / 0 | 1 / 0 |

- **Primary endpoint positive under both.** None of the 19 overlaps a curated element.
- **What the boundary stages cost here.** No candidate loses its call. Two fall under 50% coverage because an end moved inward onto a smaller region: an att pair trims a 19.2 kb candidate's call to 9.0 kb (`GCF_000196155.1`, 4.73 Mb), and composition trims another (`GCA_017183655.1`, 0.61 → 0.41 covered). Both are among the three candidates only prophound covered, and both are flagged novel by LAMBDA. The candidates' coordinates come from language-model windows (most ends sit on a 1 kb step), so which region is right is not known. The count-matched figure falls for the same reasons and because refinement splits calls.
- **The annotation tools reach most of them too.** geNomad covers 15 of 19, PIDE 0.789, PhageBoost 0.778. These are mostly phages the database tools call but the curators left out, not elements every tool missed.
- It also covers more of what reviewers doubted (UNLIKELY PHAGE 0.55 under 0.1.0, 0.64 before, against geNomad's 0.36), which is the price of its call volume.

Wording: *"prophound's calls overlap all 19 phages LAMBDA's reviewers judged likely but its curated set lacks, and cover at least half of 16 of them, above chance (0.75, 0.54–0.95), at 15 calls/genome. geNomad covers 15."* Not "phages the annotation tools cannot find", and not discovery: the grades are review judgements.

### 10. Boundaries: where each end comes from, and how far to trust it

Each end of a call has one source, reported per call (`start_source`, `end_source`). The stages run in this order, and each keeps the previous one's end where it does not apply (`docs/att-boundaries.md`):

1. **Refinement** (all modes): a 5 kb scoring pass at 250 bp steps moves the coarse 1 kb edge.
2. **Att site** (all modes): an exact direct repeat (core ≤300 bp, at most 2 copies, E ≤ 1e-6) with one copy within 10 kb inside / 5 kb outside each end. Where found, the region runs from the first base of the left copy to the last of the right, and that end is reported as exact.
3. **Composition** (every mode; Exploratory accepts profiles that differ less, divergence ≥ 0.005 against 0.01): the host-to-element tetranucleotide change-point inside the end's range, moved inward by at most min(15% of the call, 8 kb).

The att rule, the composition gates and the range fractions were chosen on the tuning sample. They were checked on the second sample and the literature pilot (both used before for other decisions) and on LAMBDA (in view for the 8 kb ceiling). None of these is held out for the boundary stages. The second lockbox is, and its preregistered boundary endpoints are in the last column.

| Default | tuning (42) | second sample (18) | literature (7) | LAMBDA (39) | **second lockbox (150), held out** |
|---|---|---|---|---|---|
| covered50, before → 0.1.0 | 0.972 → 0.972 | 1.000 → 1.000 | 0.912 → 0.912 | 0.933 → 0.933 | 0.989 (0.1.0) |
| median IoU, before → 0.1.0 | 0.672 → 0.812 | 0.666 → 0.821 | 0.502 → 0.724 | 0.810 (0.1.0) | 0.742 (0.1.0) |
| median end error, 0.1.0 | 2.1 kb | 2.1 kb | 1.8 kb | 1.9 kb | not an endpoint |
| ends within 1 kb, 0.1.0 | 33% | 32% | 36% | 40% | not an endpoint |
| range coverage, Default and Conservative ends (fitted as 90%) | 0.901 | 0.897 | 0.900 | 0.870 | **0.780** (Default 0.782, Wilson 0.732–0.825) |
| calls with an att end | | | | | 0.084 |

Truth is geNomad's coordinates on the 42 and the 18, authors' coordinates on the literature pilot, and LAMBDA's curation on LAMBDA. Errors against geNomad include geNomad's own gene-edge imprecision. On the tuning sample, Default calls 17% fewer Mb (20.7 → 17.3) and 22% fewer Mb off every element (14.9 → 11.6).

- **Att ends are base-exact where found, but rarely found.** Against authors' coordinates (literature pilot), att ends sat a median 22 bp away, 79% within 100 bp; most differences are exactly one core copy, a convention the core length in the output lets a user convert. Only 14–16% of calls on known elements get an att end, on all three samples. On the tuning sample, 85% of single-element calls have no exact repeat ≥16 bp near both of geNomad's edges, so this is a property of the genomes, not of the rule. From a decoy search on host sequence, an estimated 12–17% of Default att ends are chance repeats (Conservative 4–6%, Exploratory 27–32%). Against geNomad, att ends sit a median 0.5–0.7 kb away, which is geNomad's gene-edge imprecision.
- **Att cores carry the integration signature, which the search never sees.** In Default, 0.577 of att calls (94 of 163, Wilson 0.500–0.650; four samples) have a core copy on a tRNA or tmRNA gene. Random placement would give 0.004, and host repeats passing the same rule give 0.045 (3 of 66). In 91 of the 94, only one copy lies on the gene, as integration leaves it, and the genes are the known hotspots (tRNA-Arg, -Met, -Leu, tmRNA). Att calls off every known element carry it at 0.495, so most are real integration sites geNomad does not call, of prophages or other integrative elements (`docs/att-boundaries.md` "Are the att sites real?"; analysis only, no genes in the pipeline). This is evidence about the set. A single att call without a tRNA can still be chance.
- **Estimated ends are about 2 kb off (median), with a long tail.** Composition halved the median error (3.6 → 2.1 kb on the tuning sample). In Exploratory it takes the median from 3.7 to 2.0 kb (tuning) and cuts off-element bp by 12–17%, at the cost of one tuning-sample element (a partial detection just over 50%; 0.994 → 0.989), unchanged on the other three samples.
- **Ranges.** An estimated end's range runs 0.143 × the call's length outward and 0.463 × inward in every mode, fitted as the 5th–95th percentiles on the tuning sample (Default and Conservative). Coverage on the other samples is 0.88–0.92 in every mode; on LAMBDA it is lower (Default 0.845, Exploratory 0.804). Median width is 21–26 kb, so the ranges are honest but wide. Ends nothing could move stay on the 1 kb grid; they are flagged in `note` and their ranges held less often. **Held out they hold less: 0.782 of Default model edges on the second lockbox (237 of 303; Wilson 0.732–0.825), 0.765 in Exploratory and 0.777 in Conservative.** By the rule fixed before scoring they are therefore not called 90% ranges. The release states what they did: they held the true end in about 78% of cases on held-out genomes (88–92% on the development samples). The fractions were fitted on the tuning sample and are too narrow for new genera; widening them is a post-0.1.0 change.
- **Costs.** Calls per genome rise 1–5% (refinement splits some calls). That lowers count-matched comparisons (§4) and part of the LAMBDA candidates figure (§9). Run time rises 10–23% on the M4 GPU (Default; CPU not yet measured; "Run time" below). On the literature pilot, Exploratory and Conservative each lose one element (a partial detection trimmed under 50%).

Wording: *"boundaries are placed without genes: base-exact at an att site where one is found (about one call in seven on known elements), otherwise estimated to a median of about 2 kb, with a range on every estimated end that held the true end in about 78% of cases on held-out genomes"*. Not "base-exact boundaries" unqualified.

### 11. What each run reports

- **`regions.tsv`**, one row per candidate, ranked, with 1-based inclusive coordinates as genome browsers show them. Columns: `rank`, `call_id`, `contig`, `start`, `end`, `length_kb`, `edges` (`att site` or `estimated`), `start_range` and `end_range`, `support` (how many of the three checkpoints make the call) and `note`. The terminal prints the same list at the end of a run.
- **`regions.detailed.tsv`** holds every column, 0-based half-open, for scripts. It includes each end's source, the model's own edges after refinement but before att and composition, the att core, scores and flags. BED and GFF3 versions are written alongside, and `regions.filtered.tsv` lists the candidates the seed vote or the N filter set aside, with the reason.
- **Rank** is the call's highest window score, then its mean score. On the development samples this ordering put known elements first more often than the mean score alone, the earlier documented choice, or ranking att calls first. A truly novel element can therefore rank low.
- **`docs/interpreting-results.md`** is the reading guide for researchers: what each column means, how far to trust each kind of end, and how candidates feed into lab checks. Every file and column is defined in `docs/output-files.md`.

---

## Operating modes

The ensemble with the shipped seed vote. Cost is calls per genome, the follow-up budget. Rows marked † were measured before the boundary stages.

| | Exploratory | **Default** (recommended) | Conservative |
|---|---|---|---|
| rule | score 0.90, every call kept, flagged by `seed_support` | 0.95, all 3 checkpoints agree | 0.995, all 3 agree |
| short contigs (3–10 kb) | 0.95 | 0.95 | 0.99 |
| **second lockbox covered50 / calls per genome (held out)** | 0.984 / 37.4 | **0.989 / 12.4** | 0.903 / 5.0 |
| second lockbox, genomes with no geNomad element (77), calls per genome | 31.6 | 9.2 | 2.7 |
| first lockbox covered50 / calls per genome (pre-boundary defaults) † | 0.965 / 32.4 | 0.960 / 11.0 | 0.901 / 4.9 |
| tuning sample covered50 / calls per genome | 0.989 / 36.1 | 0.972 / 14.4 | 0.939 / 8.6 |
| geNomad-negative genomes, calls per genome † | 27.1 (19.1–35.5) | 7.4 (4.7–10.2) | 2.45 (1.35–3.7) |
| LAMBDA covered50 / calls per genome | 0.946 / 40.3 | 0.933 / 14.7 | 0.899 / 7.3 |
| literature pilot, all 57 | 0.947 | 0.912 | 0.877 |
| small inovirus implants, clean sites (covered50) † | 0.617 | 0.362 | — |
| median IoU, tuning sample | 0.812 | 0.812 | 0.838 |
| ends from composition | yes (looser gate) | yes | yes |

- **Exploratory** is for discovery runs: the widest coverage, at about 2.5–3× Default's calls. On the second lockbox it touched every element (covered 1.000), but its covered50 was one element below Default's (0.984 against 0.989). Its looser composition gate trimmed two partial detections under half (`docs/lockbox2.md` "Result"). Calls that only one checkpoint makes are 99% off-target relative to geNomad (`docs/error-analysis.md`), so `seed_support` is the first column to triage by.
- **Default** is the recommended mode.
- **Conservative** is for a small follow-up budget. It costs recall first on 10–15 kb elements (0.900 †) and small elements. On the tuning sample and the literature pilot, composition edges cost it one element each (0.944 → 0.939, 0.895 → 0.877): partial detections covered barely over 50%.

---

## What 0.1.0 does not claim

- **No validated novel elements.** No call has been confirmed in the lab, and enrichment for MGE features is not discovery.
- **No claim that it is more accurate than** geNomad, PHASTER or any tool. Per unit of called sequence it is not shown to be ahead of geNomad, PHASTER or PIDE, and at matched call counts it trails them (§2, §4).
- **No comparative claim on fragmented input** (§6).
- **No claim about the mechanism.** Composition alone does not explain the model, but which biology it uses is unmeasured (item 1b). Containment is nucleotide-level and cannot see protein homology.
- **No ≥95% claim per size stratum, for Exploratory or Conservative, or against curated truth.** The held-out ≥95% claim is pooled, Default mode, and against geNomad (§1). No measured claim below 10 kb from geNomad truth.
- **No claim that the edge ranges are 90% ranges.** Held out they held the true end in 78% of cases (§10).
- **No claim of base-exact boundaries in general.** Ends are base-exact only at an att site, found on about one call in seven on known elements, and an att site can be a chance repeat (§10). Estimated ends carry ranges.
- **No "first" and no "guarantee"**.

## Limitations

- **Labels are geNomad-derived.** Training and every agreement figure use geNomad 1.9.4 with DB 1.9. The model inherits geNomad's misses, and agreement figures measure fidelity to it. geNomad's ~3.5 kb floor means geNomad truth holds almost no elements under 10 kb (3 of 202 in the lockbox). Of 89 confident inovirus loci in the corpus, geNomad labels 43 negative, median 5.6 kb (`docs/data-reliability.md`).
- **Small elements** †. Real inoviruses implanted at model-independent sites are detected at **0.362** at Default (Wilson 0.317–0.409, n = 420; implants share 28 hosts, so the true interval is wider): 0.17 at 3–5 kb, 0.32 at 5–8 kb, 0.73 at 8–10 kb (C7). Evidence below 10 kb is partly synthetic and covers one phage family.
- **Boundaries are still coarser than gene-based tools'** where no att site is found, which is most calls. Estimated ends are a median ~2 kb off, with a tail of partial detections that the wide ranges reflect. In the literature pilot, Default median IoU is 0.73 / 0.62 (test / val), against 0.89 / 0.71 for geNomad. On LAMBDA, Default reaches 0.81, and the annotation tools 0.83–0.92 (§10).
- **Conflation.** Several elements may fall under one call. On the tuning sample, Default calls absorb 11 of 179 elements (6.1%) under 0.1.0, 13 before (refinement splits some). On the second lockbox, 5 calls absorb 10 of 186 elements (5.4%), the same in every mode (first lockbox †: 19 of 202, 9.4%). Covered recall cannot see conflation, which is why conflated and absorbed counts are reported.
- **Cost.** Default makes 14.4 calls per genome on the tuning sample and 14.7 on LAMBDA (12.4 on the second lockbox), against 3.6–8.4 for annotation tools on LAMBDA. On genomes with no geNomad provirus it makes 7.4 per genome on the 20 geNomad-negative genomes (4.7–10.2 †) and 9.2 on the second lockbox's 77, which is the floor cost.
- **Short contigs** † under 3 kb are never scanned. Recall on 3–10 kb contigs falls with element size: on fully shredded input, 60 kb+ elements read 0.286 at covered50, against 0.905 on complete genomes (calibration run `818235`, before the final modes; `docs/adaptive-multiscale-detection.md` Part 4).
- **Large internal N gaps** † cost recall: 6 of 179 elements at 20 kb.
- **Input outside bacteria.** Human sequence is called heavily: 57% of a 10 Mb chr21 slice, at 1.9× the bacterial call rate. Yeast is near the bacterial rate. Three eukaryotic viruses get no calls. On archaea, two hyperthermophiles have more than a third of the genome called. Input should be bacterial, and archaeal calls are untested (`docs/robustness-checks.md` check 4).
- **Plasmids are often called whole.** On the second lockbox, 18% of plasmid sequence is called against 7% of chromosomal sequence, and a single call covers most of the replicon on 32% of plasmids. Such calls describe the replicon, not an element inside it. All 7 geNomad elements on plasmids are covered (`docs/robustness-checks.md` check 3).
- **Calibration.** Mode thresholds are a user-facing trade-off between recall and cost, tuned on geNomad agreement. They are not probabilities that a call is real.
- **Sample size.** Every external benchmark is small and dominated by a few element-rich genomes: LAMBDA 39 genomes, the literature pilot 7, G5b 9 (5 *Paracoccus*). Every interval is wide, and results are read through paired comparisons.
- **Shredding** uses uniform breakpoints. Real assemblies break at repeats, so the MAG-like conditions approximate MAGs rather than reproduce them.

## Assumptions, and the evidence for each

| assumption | why it matters | evidence | status |
|---|---|---|---|
| geNomad truth is a usable proxy for "known element" | every agreement figure | geNomad covers 0.822 of LAMBDA's curated elements, so it is a proxy, not truth; hence two-layer reporting | **partly holds**; stated |
| elements within a genome are not independent | interval width | bootstrap and Wilson differ, sometimes in either direction (C2) | handled: genome bootstrap is primary |
| covered50 is not a flattering cut-off | the recall definition | the mode order is the same at ≥1/25/50/75%, and ≥75% costs 1–3 points (C3) | **holds** |
| near-duplicates of training do not carry recall | generalization | tuning-sample <0.1 containment bin reads 0.964; clade-novel recall is flat to ~20% divergence; lockbox median containment 0.003 | **holds at the nucleotide level**; protein level untested |
| the model uses more than composition | that it is not a GC detector | AUPRC ≤0.16 against 0.60 (1c) | **holds** at this budget (k ≤ 4, no ORF density) |
| calls per genome is the right cost unit | how modes are priced | a product decision (user, 2026-09-21): wet-lab follow-up scales with calls | design choice, stated |
| a random-rotation baseline is fair chance | §2's cost adjustment | rotation keeps each tool's call lengths, count and spacing; reproduces the original estimate (0.141 vs 0.144) | stated method; one choice among several |
| the tuning sample's optimism is small | how far reused numbers transfer | reused minus held out, second lockbox (the 42 / the 18): Default −0.017 / +0.011, Exploratory +0.010 / +0.016, Conservative +0.041 / +0.059 (first lockbox, Default: +0.012 / +0.040) | **holds for Default and Exploratory; larger for Conservative.** Measured by two lockboxes |
| shreds approximate MAG input | §6 | uniform breakpoints; N50 not calibrated to a catalogue | **untested**; stated |
| an exact direct repeat at both ends of a call marks the integration site | att ends replace the model's | decoy search: an estimated 12–17% of Default att ends are chance repeats; on author-curated elements, att ends match to a median 22 bp; independently, 0.58 of Default att cores sit on a tRNA/tmRNA gene against 0.004 at random and 0.045 for decoys (§10) | **holds for the set**; the chance share per call is stated, and the core is reported so a user can check it |
| a composition change-point near an estimated end is the element's edge | estimated ends, every mode | median end error halves (3.6 → 2.1 kb, tuning); it fails when the call's middle is a different element, which the 8 kb inward ceiling limits (`docs/att-boundaries.md`) | **holds on the samples measured**; held out, estimated ends were bracketed less often (range coverage 0.78 against 0.88–0.92) |
| one inovirus family stands in for small prophages | the small-element limitation | none: the corpus is mostly tailed phages | **untested**; stated as a limit |

## Footprint (G4, reported only alongside recall)

| run-time data | size | source |
|---|---|---|
| prophound: the bundled model (3 weights-only checkpoints, config, manifest) | **16.5 MB** (16,500,331 bytes) | the package's `model/` directory, 2026-10-20 |
| geNomad database v1.9 | **1.47 GB** unpacked (1,471,630,421 bytes); 842 MB compressed | `du -sb` on ARC, 2026-10-16; Zenodo 14886553 |

prophound needs about 89× less run-time data than geNomad. That matters only because it comes with the recall in §1–§3. A small tool that found less would not be better. The training checkpoints were 49.5 MB because they also held optimizer and RNG state. The bundled files hold the same weights, tensor for tensor, and give byte-identical outputs (`docs/lockbox.md`). Installed-environment size is not compared yet: prophound's runtime environment is dominated by PyTorch (~381 MB on osx-arm64).

## Run time and hardware

Measured on a MacBook Pro (Apple M4, 16 GB). This is the kind of machine most users will run it on, so it is the headline. Default mode, 3-checkpoint ensemble. The Apple GPU rows are the 2026-10-20 re-run after the Phase 3 speed work; the batch and CPU rows are from 2026-10-19 (`docs/benchmarks.md`). **Every row is †:** measured before refinement became a default. Refinement added 10–23% on the M4 GPU in Default (14.2 → 15.6 s at 2.2 Mb, 48.6 → 59.6 s at 7.2 Mb; `docs/att-boundaries.md`); CPU is expected to rise about 20% but has not been measured. The table is re-run with the 0.1.0 defaults before release.

| | wall time | memory |
|---|---|---|
| typical 5 Mb complete genome, Apple GPU (`--device auto` picks it) | **~35 s** (33.5 s median, 4.97 Mb) | ~0.5 GB |
| the same genome as a MAG-like assembly (374 contigs) | ~19 s | ~0.5 GB |
| 10 genomes, 49 Mb | 7.4 min | ~0.6 GB |
| typical 5 Mb genome, CPU only | **~45 min** (4.97 Mb, 2 runs: 45.0 min) | ~1.4 GB |
| on an NVIDIA V100S (cluster) | ~27 s for a median 3.6 Mb genome (4.5 s + 6.1 s/Mb) | |

On the same 10 genomes and the same laptop, Jaeger 1.26.5 took a median of 207 s with its G5 primary model and 56 s with `default`. It ran on CPU, because its TensorFlow Metal plugin does not load. geNomad took 342 s on a cluster CPU node (4 threads, different hardware). These are context, not a G6 claim.

**Devices.** Calls were the same on Apple M4 CPU and GPU, and scores agree with an NVIDIA V100S to about 1e-5. Outputs are not byte-identical across devices, because floating-point order differs. A call can differ only where a score sits within about 1e-5 of a threshold. Across 45 test assemblies (134,070 windows), one internal decision differed between the Apple GPU and the V100S, and no released call did. Other platforms (x86 CPUs, other GPUs, Windows) are expected to agree but are untested.

## Not yet measured (blocks wording, not the release)

- **CPU runtime (G6) against geNomad on one machine.** prophound's own hardware baseline is measured (above, `docs/benchmarks.md`). The G6 comparison still needs geNomad timed on the same machine; G5's geNomad timings came from ARC.
- **The batch row after Phase 3.** The single-genome GPU rows were re-run on 2026-10-20; the 10-genome batch was not.
- **Installed-environment size** for both tools, once the bioconda package exists.
- **Run time with the 0.1.0 defaults** on the M4 GPU and CPU (above).

## Reproduction

| figure | script | output |
|---|---|---|
| lockbox | `scripts/lockbox_score.py` | `analysis/lockbox/results.{json,elements.csv}` |
| tuning / second sample, intervals | `scripts/integrity_checks.py recall` | `analysis/integrity/recall{42,18}.json` |
| containment | `scripts/eval_containment.py` | `analysis/integrity/containment*.tsv` |
| LAMBDA | `scripts/lambda_score.py --containment` | `analysis/external_benchmarks/lambda/scores.*` |
| literature pilot | `literature_pilot.slurm`, `integrity_checks.py lit` | `analysis/litpilot/results_final.*`, `analysis/integrity/lit_sensitivity_final.*` |
| G5 / G5b | `scripts/g5_score.py` | `analysis/g5/results.*`, `analysis/g5b/results.*` |
| LAMBDA chance adjustment | `scripts/lambda_chance.py` | `analysis/external_benchmarks/lambda/chance.json` |
| off-target enrichment | `offtarget_annotation.py --control_busy elements_only`, then `scripts/offtarget_enrichment_ci.py` | `analysis/offtarget/annotation_ctrl_elements_only.*`, `enrichment_ci.json` |
| small inoviruses | `small_element_spikein.py` | `analysis/spikein_siterule/results_final.*`, `analysis/spikein_small/results_final.*` |
| LAMBDA candidates | `scripts/lambda_candidates.py` | `analysis/external_benchmarks/lambda/candidates.{json,csv}` |
| LAMBDA, chance and candidates under 0.1.0 | the three LAMBDA scripts with `--prophound_run arc_data/att_runs/lambda --genomes_dir analysis/g5/fasta` | `analysis/external_benchmarks/lambda/{scores,chance,candidates}_0.1.0.*` |
| literature pilot under 0.1.0 | `scripts/literature_pilot_eval.py --pipeline_run arc_data/att_runs/lit7`, then `integrity_checks.py lit` on its `ens_0.1.0` runs | `analysis/litpilot/results_0.1.0.*`, `analysis/integrity/lit_sensitivity_0.1.0.json` |
| boundary stages (§10) | `scripts/boundary_eval.py`, `scripts/edge_composition.py` and the cap sweep | `analysis/att/*_pipeline.csv`, `analysis/att/composition_cap_sweep.csv` |
| att cores on tRNA genes (§10) | `scripts/att_trna_check.py` (ARAGORN 1.2.41 via `pixi exec`) | `analysis/att/trna/att_trna.{summary,pairs}.tsv` |

The training checkpoints' SHA-256s are in `docs/research-integrity.md` C8, and the defaults the first lockbox scored are frozen at `09c7756`. The 0.1.0 defaults scored by the second lockbox are frozen at `24f6c44` (`docs/lockbox2.md`). The bundled model's SHA-256s are in its `manifest.json`, mapped to the C8 hashes in `docs/lockbox.md`. Every figure above whose inputs are local is re-run by `scripts/reproduce_figures.sh` and checked against these outputs by `scripts/compare_figures.py`; all reproduce exactly (2026-10-20). The 0.1.0 re-scores (the last three rows) are not in that harness yet; their inputs are the local refinement-on runs in gitignored `arc_data/att_runs/`.
