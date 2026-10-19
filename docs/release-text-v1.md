# prophound v1 — release text (DRAFT)

> **Draft, 2026-10-16; lockbox figures filled in 2026-10-19** from the one scoring of `841208` (`analysis/lockbox/results.json`, `docs/lockbox.md` "Result"). Every lockbox figure below comes from that scoring. Nothing else changes when it lands, except wording the decision rule ties to it (marked ⟂). Rules for this text: every claim carries its number, *n*, interval, sample role and source. Wording follows `docs/claims-and-evidence.md`'s ladder. No "first", no "guarantee", "may not call" rather than "do not call", and calls are "candidates for follow-up".

---

## What prophound is

*prophound is a recall-first prophage screen. It localizes candidate integrated elements directly from nucleotide sequence, with no gene calling or reference database at run time. It is tuned to recover already-known elements and returns everything else as candidates, at a follow-up cost the user chooses. It is designed to run alongside annotation-based tools and to surface what they may not call.*

**Method in one line:** an exploratory, database-free prophage caller that learns directly from sequence, trained on geNomad-labelled genomes. It is a reverse-complement-symmetric CNN, scored over 10 kb windows by an ensemble of three training replicates. Regions are built from the score track with no gene calling.

**What is distinct is the contract, not the method.** Database-free and nucleotide-level detection already exist (PhageBoost, PIDE, Jaeger, Seeker, LAMBDA's genomic language models; `docs/related-tools.md`). prophound's contract has four parts:
1. a recall target on known elements, measured and reported with intervals;
2. three modes, each priced in calls per genome (the follow-up cost), not in precision;
3. output as candidates, not annotations;
4. use alongside annotation tools.

It rests on three design facts, each checked against the nearest tools:
- **Host-context training.** It is trained on prophages in their host genomes. Jaeger and LAMBDA's models train on free phage genomes, with Jaeger masking prophages out.
- **Gene-free boundaries.** Boundaries come from the score track. Jaeger snaps them with gene prediction.
- **Fragments scanned by default.** It scans contigs down to 3 kb. Jaeger's prophage mode scans contigs over 500 kb by default.

## How to read every number below

- **Recall metric: covered50.** An element counts as found when calls cover at least 50% of its bp. ≥1 bp would credit grazes. Across the ≥1 / ≥25 / ≥50 / ≥75% cut-offs, every one orders the modes identically, and ≥75% costs only 1–3 points, so covered50 is not a flattering choice (`docs/research-integrity.md` C3).
- **Intervals.** 95%, by genome bootstrap (10,000 resamples), because elements within a genome are not independent. A Wilson interval is added where a bootstrap degenerates (an estimate of 1.000) or *n* is small.
- **Cost: calls per genome.** It is the follow-up (wet-lab) budget. Precision against geNomad is not reported as accuracy. In a discovery tool, a call geNomad does not make is the intended output, not the error term.
- **Sample roles** (C1). Each figure is labelled by the role of the sample it was measured on:

| sample | role | wording used |
|---|---|---|
| lockbox: 150 genus-disjoint genomes, 202 elements | **held out from training and tuning**, preregistered, scored once | "held out" |
| 42 genomes / 179 elements | tuning: thresholds, merging, vote and modes were chosen on it | "tuning sample" |
| 18 genomes / 52 elements | confirmation, but it informed earlier decisions | "a second, previously used sample" |
| 20 genomes with no geNomad provirus | confirmation; in view for the vote decision | "geNomad-negative genomes" |
| LAMBDA clean subset: 39 genomes / 149 curated elements | external curated truth; bars fixed before the run | "LAMBDA" |
| literature pilot: 7 genomes / 57 elements | external curated truth; most contained in training (median 0.723) | "literature pilot" |
| Philympics G5b: 9 genomes / 32 elements | external curated truth; preregistered | "G5b" |

- **Two layers of accuracy.** Agreement with geNomad is not accuracy, because geNomad covers 0.822 of LAMBDA's curated prophages. Both layers are reported, each labelled as what it is.
- **Novelty to training** is 21-mer containment of an element in training positives. It sees nucleotide similarity only, not protein homology.

---

## Claims

### 1. Agreement with geNomad on held-out genomes

⟂ **Default covered50 0.960 (194 of 202; Wilson 0.924–0.980; genome bootstrap 0.929–0.984), n = 202 elements in 80 genomes, at 11.0 calls per genome.** Exploratory 0.965 (0.930–0.983 Wilson; 0.937–0.988 bootstrap) at 32.4, and Conservative 0.901 (0.852–0.935; 0.860–0.941) at 4.9. Lockbox, held out, preregistered primary endpoint. Covered (≥1 bp): 0.975 / 0.980 / 0.916. Default's 8 misses fall in 8 different genomes; 5 have no overlapping call and 3 are covered 18–39%.

The lockbox is nearly novel to training by construction: median containment 0.003, with 6 of 202 elements at ≥0.5 (C5). It cannot speak below 10 kb (3 elements). Those rows are reported as n = 1, not as claims.

For comparison, on the reused samples (Default, ensemble):

| sample | covered50 | 95% interval | calls/genome |
|---|---|---|---|
| tuning (42 / 179) | 0.972 | 0.949–0.994 (bootstrap); 0.936–0.988 (Wilson) | 14.1 (12.1–16.3) |
| second sample (18 / 52) | 1.000 | 0.931–1.000 (Wilson) | 13.3 (9.7–17.3) |
| lockbox (150 / 202) | 0.960 | 0.929–0.984 (bootstrap); 0.924–0.980 (Wilson) | 11.0 |
| lockbox gap to the 42 / the 18 | +0.012 / +0.040 | | |

- **≥95% is a point estimate on the reused samples, not an established floor.** The tuning sample's lower bound is 0.949. ⟂ **The lockbox does not establish it:** its lower bound is 0.924 (Wilson) / 0.929 (bootstrap). By the rule fixed beforehand (≥95% only if the lockbox's lower bound is ≥0.95), the release gives the number and its interval, 0.960 (0.924–0.980), and does not say ≥95%. The optimism of the reused samples is small for Default: the lockbox reads 0.012 below the tuning sample, well inside either interval.
- **Not carried by training look-alikes.** On the tuning sample, the 137 elements with containment <0.1 read 0.964, and all five Default misses are in that bin (C5). Covered recall on clade-novel elements (MinHash, k = 21) is flat at 0.966–0.975 from ~3.5% out to ~20% nucleotide divergence, with the same four misses at every threshold. That was measured under the 2026-09-22 builder, at ≥1 bp (`docs/scientific-validity.md` 2a).
- **By size, tuning sample, Default** (Wilson intervals; `analysis/integrity/recall42.json`). The ≥95%-per-stratum target is an engineering target; no stratum below 30 kb establishes it:

| stratum | n | covered50 | Wilson 95% |
|---|---|---|---|
| 5–8 kb | 7 | 1.000 | 0.646–1.000 |
| 8–10 kb | 11 | 1.000 | 0.741–1.000 |
| 10–15 kb | 30 | 0.933 | 0.787–0.982 |
| 15–30 kb | 41 | 0.951 | 0.839–0.987 |
| 30–60 kb | 69 | 0.986 | 0.922–0.997 |
| 60 kb+ | 21 | 1.000 | 0.845–1.000 |

- **By size, lockbox, Default** (Wilson intervals; `analysis/lockbox/results.json`). Below 10 kb n = 1 per stratum, reported as counts, not claims:

| stratum | n | covered50 | Wilson 95% |
|---|---|---|---|
| <5 kb | 1 | 1 of 1 | — |
| 5–8 kb | 1 | 1 of 1 | — |
| 8–10 kb | 1 | 1 of 1 | — |
| 10–15 kb | 18 | 1.000 | 0.824–1.000 |
| 15–30 kb | 49 | 0.918 | 0.808–0.968 |
| 30–60 kb | 99 | 0.970 | 0.915–0.990 |
| 60 kb+ | 33 | 0.970 | 0.847–0.995 |

15–30 kb is the weakest lockbox stratum (4 of the 8 misses); on the tuning sample it was the second weakest.

Wording: *"agrees with geNomad on X of held-out elements (covered50, 95% interval Y–Z)"*. Never "X% accurate".

### 2. Recall on expert-curated truth (LAMBDA)

On LAMBDA's clean subset (no genome in train or val), Default covers **0.933** of 149 curated elements (genome bootstrap 0.874–0.983), at 14.4 calls/genome. The primary row, 34 genomes not in the corpus, reads **0.920** (113 elements; 0.859–0.983). For the published tools on the same subset, rescored with the same rule: geNomad 0.822, PHASTER 0.824, PIDE 0.824, best genomic LM (EVO2 4k) 0.785. The bars were fixed before the run (`docs/claims-and-evidence.md` "LAMBDA result").

**It is costlier, and the cost explains much of the margin.** Default calls 9.8% of each chromosome, geNomad 3.7%. Rotating each tool's calls by a random offset per chromosome keeps their lengths, count and spacing. It gives chance covered50 of 0.141 for Default and 0.052 for geNomad (1,000 draws). Recall above chance is therefore **0.792 for Default against 0.770 for geNomad**, 0.772 for PHASTER and 0.776 for PIDE.

**Paired on the genomes both tools ran, Default minus geNomad is +0.018 (genome bootstrap −0.037 to +0.078).** Against PHASTER it is +0.019 (−0.033 to +0.080), and against PIDE +0.015 (−0.045 to +0.081). Every interval spans zero, so per unit of cost prophound is level with the annotation tools, not ahead. Of 8 pairings (Default and Conservative against four tools), one interval excludes zero, Conservative minus VIBRANT at +0.067 (+0.001 to +0.128). That is what multiple comparisons would produce, and it is not claimed. Script: `scripts/lambda_chance.py`.

- **Supported:** the highest raw curated recall of the tools compared, *at 1.5–3.7× their calls per genome*.
- **Not supported:** "more accurate than" any tool.

### 3. Run alongside geNomad, it adds curated prophages (G1)

geNomad plus prophound Default recovers **16 more curated elements than geNomad alone, and loses none** (+16 / −0). On the 37 genomes both were run on, that is +0.110 covered50 (paired genome bootstrap 0.018–0.210). The cost is +14.6 calls/genome. The preregistered bar was more than +5, the best existing pairing with geNomad.

- **The gain is concentrated.** 12 of the 16 come from two genomes (Photorhabdus TTO1, 9; *P. aeruginosa* LESB61, 3), and only 5 genomes gain anything. The interval excludes zero, but the effect is driven by few genomes.
- **The gains are not training look-alikes.** All 16 have containment ≤0.16, and 13 of them are below 0.01.

Wording: *"run alongside geNomad, prophound recovered 16 additional curated prophages on LAMBDA's clean subset (+16 / −0), at about 15 extra candidate calls per genome"*.

### 4. Against sequence-only models at matched cost (G2)

Exploratory's calls, ranked by score and cut to EVO2 4k's per-genome call count, reach **0.805 against EVO2 4k's 0.785**, so the preregistered bar (>0.785) is cleared narrowly. It is reported with its limits:
- at geNomad's call count it ties geNomad (0.822 vs 0.822);
- at PHASTER's lower count it trails (0.797 vs 0.824);
- no paired interval separates it from EVO2 4k.

Wording: *"at matched follow-up cost, comparable to or slightly above the best published genomic language model on LAMBDA"*. Not "exceeds".

### 5. Literature-curated pilot, split by novelty and class

There are 57 elements in 7 genomes, curated from the literature, with geNomad scored identically (`docs/final-benchmark.md`, C4 final, C5). **The pooled figure is flattered**, because 34 of the 57 elements are ≥50% contained in training. It is never quoted without the split:

| subset | n | Default | geNomad | note |
|---|---|---|---|---|
| all | 57 | 0.912 | 0.807 | paired +6 / −0, +0.105 (0.029–0.151) |
| class `prophage` (headline) | 28 | 28/28 (Wilson 0.879–1.000) | 27/28 | both at ceiling |
| novel to training (<0.1) | 13 | 8/13 (0.355–0.823) | 6/13 (0.232–0.709) | intervals overlap |
| remnants | 12 | 7/12 | 3/12 | where the difference lives |

- The six gains are Sp16, CP-933L, SpLE1, SpLE3, SpLE5 and *B. subtilis* prophage 6, and five of them are remnants or cryptic prophages.
- A permutation test against randomly placed calls is significant on both splits (Fisher-combined p ≤ 7 × 10⁻⁵ on covered50; `821391`, interim modes).
- Boundaries are clearly worse than geNomad's: median IoU 0.59–0.62 against 0.71–0.94.

Wording: *"on literature-named prophages both tools are near ceiling; prophound additionally covers more remnant and cryptic prophage regions"*. This is an example, not a rate (7 genomes, two of them *E. coli* O157:H7).

### 6. Fragmented input: a design fact and a self-referential figure, no comparison

- **Design fact:** contigs of 3–10 kb are scanned by default with calibrated thresholds. Contigs under 3 kb are never scanned.
- **G5 (LAMBDA shreds, MAG-like N50 20 kb):** Default kept **0.916** (0.844–0.968) of its complete-genome recall (absolute 0.855, at 17.6 calls/genome). geNomad kept 1.022 and Jaeger 0.987. prophound's retention is lower, with the paired intervals below zero. **The comparative claim is refuted on LAMBDA.**
- **G5b (Philympics, held out, preregistered):** Default kept **0.938** (0.889–0.983), against geNomad 0.958 (difference −0.021, −0.074 to +0.026) and Jaeger 0.881 (+0.057, −0.064 to +0.165). It does not clear. On draft-like input (N50 100 kb), Default kept 1.000.
- **Fully shredded input** (tuning sample, every contig under 10 kb): Default covered50 0.626 at 20.1 calls/genome, Exploratory 0.827.
- **The shipped edge snap** extends calls within 5 kb of a contig end onto the end, and flags them `contig_edge`. On G5b it changed no element's covered50 status and added 6–7% more called sequence on MAG input, 78–88% of it off-element. It ships as frozen. Edges it moves are flagged, so users can discount them.

Wording: *"on MAG-like fragmented input, prophound Default kept 0.92–0.94 of its complete-genome recall (LAMBDA; Philympics), against geNomad 0.96–1.02 and Jaeger 0.88–0.99 on the same input"*. No comparative claim.

### 7. Robustness (self-referential)

| property | result | n / source |
|---|---|---|
| reverse complement | window scores exact; region calls identical on whole contigs | `docs/scientific-validity.md` standing check |
| ambiguity codes (IUPAC, incl. N), `impute` default | ≤2 of 179 elements lost up to 3,000 ppm scattered codes; the old encoding lost 172 | `821396`, tuning sample |
| long internal N gaps | a 20 kb gap costs 6 of 179 elements at Default. Calls ≥50% N are dropped, removing exactly the spurious gap calls with no recall change at ≤5 kb gaps | `822054`, `docs/v1-scope.md` A.3 |
| grid phase | 1.7% of elements phase-unstable; ~5% call-count spread | `816447` |
| composition baseline | k ≤ 4 composition reaches AUPRC ≤0.16 against the CNN's 0.60–0.61 on the same split (random 0.014) | `816086`, item 1c |

The composition baseline shows the model uses more than nucleotide composition. It does not show *what* the model uses.

Earlier versions collapsed to zero on IUPAC codes (*E. coli* EDL933 at 1,182 ppm). That was fixed on 2026-09-24, and results differ from earlier runs only on sequence containing ambiguity codes.

### 8. What the off-target candidates are (exploratory evidence, not discovery)

**Candidates are enriched for mobile-element features, mostly not phage structural genes.**
- On the tuning sample, Default off-target calls made by all three checkpoints carry an MGE hallmark at **52.1% (99/190), against 16.8% (96/570) for length-matched controls** drawn clear of geNomad elements. That is **3.09× (genome bootstrap 2.45–4.12; 38 genomes)**, and Fisher one-sided p = 2 × 10⁻²⁰. Fisher treats the controls as independent, so the bootstrap is the interval to quote. These calls were made under the interim pre-vote Default, and the unanimous group is the one the shipped vote keeps.
- Calls only one or two checkpoints make are enriched 1.92× (1.46–2.78).
- **The enrichment is conjugation (16%) and integrase/transposase (13%) more than phage structural genes (4%).** Many candidates look like ICEs and other integrative elements, not missed prophages.

**Recovery beyond the labels: real inoviruses that geNomad labels negative.**
- Default covers **7 of 10** confident in-situ loci (Wilson 0.397–0.892); Exploratory covers 10 of 10 (0.722–1.000).
- These loci are not near-duplicates of training: median containment 0.001.
- One phage family, n = 10.

Wording: *"off-target candidates are enriched for mobile-element features (3.1× over background)"*. Not "phage", and never "novel elements".

---

## Operating modes

The ensemble with the shipped seed vote. Cost is calls per genome, the follow-up budget.

| | Exploratory | **Default** (recommended) | Conservative |
|---|---|---|---|
| rule | score 0.90, every call kept, flagged by `seed_support` | 0.95, all 3 checkpoints agree | 0.995, all 3 agree |
| short contigs (3–10 kb) | 0.95 | 0.95 | 0.99 |
| lockbox covered50 / calls per genome | 0.965 / 32.4 | 0.960 / 11.0 | 0.901 / 4.9 |
| lockbox, genomes with no geNomad element (70), calls per genome | 24.4 | 6.8 | 2.5 |
| tuning sample covered50 / calls per genome | 0.994 / 35.6 | 0.972 / 14.1 | 0.944 / 8.5 |
| geNomad-negative genomes, calls per genome | 27.1 (19.1–35.5) | 7.4 (4.7–10.2) | 2.45 (1.35–3.7) |
| LAMBDA covered50 / calls per genome | 0.946 / 39.9 | 0.933 / 14.4 | 0.899 / 7.1 |
| literature pilot, all 57 | 0.965 | 0.912 | 0.895 |
| small inovirus implants, clean sites (covered50) | 0.617 | 0.362 | — |

- **Exploratory** is for discovery runs: the most recall at about 2.5× Default's calls. Calls that only one checkpoint makes are 99% off-target relative to geNomad (`docs/error-analysis.md`), so `seed_support` is the first column to triage by.
- **Default** is the recommended mode.
- **Conservative** is for a small follow-up budget. It costs recall first on 10–15 kb elements (0.900) and small elements.

---

## What v1 does not claim

- **No validated novel elements.** No call has been confirmed in the lab, and enrichment for MGE features is not discovery.
- **No claim that it is more accurate than** geNomad, PHASTER or any tool. Per unit of cost it is level with them (§2).
- **No comparative claim on fragmented input** (§6).
- **No claim about the mechanism.** Composition alone does not explain the model, but which biology it uses is unmeasured (item 1b). Containment is nucleotide-level and cannot see protein homology.
- **No ≥95% claim, pooled or per stratum.** The lockbox's lower bound is 0.924, so ≥95% is not established (§1). No measured claim below 10 kb from geNomad truth.
- **No "first" and no "guarantee"**.

## Limitations

- **Labels are geNomad-derived.** Training and every agreement figure use geNomad 1.9.4 with DB 1.9. The model inherits geNomad's misses, and agreement figures measure fidelity to it. geNomad's ~3.5 kb floor means geNomad truth holds almost no elements under 10 kb (3 of 202 in the lockbox). Of 89 confident inovirus loci in the corpus, geNomad labels 43 negative, median 5.6 kb (`docs/data-reliability.md`).
- **Small elements.** Real inoviruses implanted at model-independent sites are detected at **0.362** at Default (Wilson 0.317–0.409, n = 420; implants share 28 hosts, so the true interval is wider): 0.17 at 3–5 kb, 0.32 at 5–8 kb, 0.73 at 8–10 kb (C7). Evidence below 10 kb is partly synthetic and covers one phage family.
- **Boundaries are coarser than gene-based tools'.** They come from 1 kb bins of the score track. In the literature pilot, median IoU is 0.59–0.62, against 0.71–0.94 for geNomad. On LAMBDA, the annotation tools reach a median IoU of 0.83–0.92.
- **Conflation.** Several elements may fall under one call. On the tuning sample, 6 calls absorb 13 of 179 elements (7.3%; ensemble, 2026-09-24 defaults). Lockbox: 9 calls absorb 19 of 202 elements (9.4%), the same in every mode. Covered recall cannot see conflation, which is why conflated and absorbed counts are reported.
- **Cost.** Default makes ~14 calls per genome on the tuning sample and 11.0 on the lockbox, against 3.4–4.6 for annotation tools. On genomes with no geNomad provirus, it makes 7.4 per genome (4.7–10.2; lockbox 6.8), which is the floor cost.
- **Short contigs** under 3 kb are never scanned. Recall on 3–10 kb contigs falls with element size: on fully shredded input, 60 kb+ elements read 0.286 at covered50, against 0.905 on complete genomes (calibration run `818235`, before the final modes; `docs/adaptive-multiscale-detection.md` Part 4).
- **Large internal N gaps** cost recall: 6 of 179 elements at 20 kb.
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
| the tuning sample's optimism is small | how far reused numbers transfer | lockbox gap to the 42 and 18 (Default): +0.012 and +0.040; Exploratory +0.029 / +0.035; Conservative +0.043 / +0.061 | **holds for Default and Exploratory; larger for Conservative.** Measured once, by the lockbox |
| shreds approximate MAG input | §6 | uniform breakpoints; N50 not calibrated to a catalogue | **untested**; stated |
| one inovirus family stands in for small prophages | the small-element limitation | none: the corpus is mostly tailed phages | **untested**; stated as a limit |

## Footprint (G4, reported only alongside recall)

| run-time data | size | source |
|---|---|---|
| prophound: 3 checkpoints | **49.5 MB** (49,538,986 bytes) | `du -sb` on ARC, 2026-10-16 |
| geNomad database v1.9 | **1.47 GB** unpacked (1,471,630,421 bytes); 842 MB compressed | same; Zenodo 14886553 |

prophound needs about 30× less run-time data than geNomad. That matters only because it comes with the recall in §1–§3. A small tool that found less would not be better. The 49.5 MB is an upper bound if the checkpoints also hold optimizer state (to be checked at packaging). Installed-environment size is not compared yet: prophound's runtime environment is dominated by PyTorch (~381 MB on osx-arm64), and both tools' environments are measured from their install recipes at packaging.

## Not yet measured (blocks wording, not the release)

- **CPU runtime (G6)** comes from the hardware baseline after the freeze (`docs/release-checklist.md` §5). G5's timings came from different hardware and are not a comparison.
- **Installed-environment size** for both tools, at packaging.

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

Checkpoint SHA-256s are in `docs/research-integrity.md` C8, and the frozen code is `ad5b28b`.
