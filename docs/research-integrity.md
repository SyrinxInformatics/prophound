# Research integrity audit

Written 2026-09-29 (user: "ensure we are working with integrity and good scientific research"). It checks the project's methods against published guidance on benchmarking and on ML validation in biology. Each finding says what the guidance asks, what we did (with the doc that records it), how serious the gap is, and a proposed correction.

**Status: reviewed by the user 2026-09-29.** Decisions:
- **Run now:** C1, C2, C3 and C8.
- **Literature truth (C4):** the agent lists every judgment-call row and runs the sensitivity analysis. The user reviews `scripts/literature_truth.tsv` before any pilot number enters the release text.
- **Lockbox (C9):** scored once, at freeze.
- **ARC corrections (C5–C7):** all three, once the cheap fixes are in (clarified by the user the same day).

## Sources

- Weber et al. 2019, *Essential guidelines for computational method benchmarking*, Genome Biology 20:125 — [link](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1738-8). Neutrality and the self-assessment trap, fair tuning, gold standards, multiple metrics, reporting.
- Whalen, Schreiber, Noble & Pollard 2022, *Navigating the pitfalls of applying machine learning in genomics*, Nat Rev Genet — [link](https://www.nature.com/articles/s41576-021-00434-9). The five pitfalls are distributional differences, dependent examples, confounding, leaky preprocessing and unbalanced classes. They recommend group (blocked) splits, since randomized cross-validation "overestimates performance because examples in the test set can be correlated with training examples".
- Kapoor & Narayanan 2023, *Leakage and the reproducibility crisis in machine-learning-based science*, Patterns — [arXiv](https://arxiv.org/abs/2207.07048). Eight leakage types: L1.1 no test set, L1.2 pre-processing on train+test, L1.3 feature selection on train+test, L1.4 duplicates, L2 illegitimate features, L3.1 temporal, L3.2 non-independence between train and test, and L3.3 sampling bias in the test distribution. They also propose a "model info sheet".
- Walsh et al. 2021, *DOME: recommendations for supervised machine learning validation in biology*, Nat Methods — [link](https://www.nature.com/articles/s41592-021-01205-4). Reporting of data, optimization, model and evaluation.
- Cawley & Talbot 2010, *On over-fitting in model selection and subsequent selection bias in performance evaluation*, JMLR — [link](https://www.jmlr.org/papers/v11/cawley10a.html). Choosing settings on the data you report is itself overfitting; the remedy is to keep selection and evaluation data separate (nested designs).
- Dwork et al. 2015, *The reusable holdout*, Science — [link](https://www.science.org/doi/10.1126/science.aaa9375). Repeatedly consulting a holdout during adaptive analysis erodes its validity.
- Bouthillier et al. 2021, *Accounting for variance in machine learning benchmarks*, MLSys — [link](https://proceedings.mlsys.org/paper_files/paper/2021/file/0184b0cd3cfb185989f858a1d9f5c1eb-Paper.pdf). Data sampling, initialization and hyperparameter choice all add variance, so estimates should vary more than just the seed.
- Joeres, Blumenthal & Kalinina 2025, *Data splitting to avoid information leakage with DataSAIL*, Nat Commun — [link](https://www.nature.com/articles/s41467-025-58606-8). Similarity-aware splitting for biological data.
- Domain benchmarks: Roach et al. 2022, *Philympics* (F1000Research, [link](https://f1000research.com/articles/10-758)), on the lack of curated gold standards for prophage callers; Schackart et al. 2023 (Front Microbiol, [link](https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2023.1078760/full)) and Ho et al. 2023 (Microbiome, [link](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-023-01533-x)), on phage-tool benchmarking design.
- Genomic interval statistics: regioneR (Gel et al. 2016, Bioinformatics, [link](https://academic.oup.com/bioinformatics/article/32/2/289/1744157)). Circular randomization preserves region clumping; uniform nulls "could result in misleading conclusions".

## What already meets the guidance

- **Group splits:** by genome, one genome per genus (Whalen pitfall 2). The clade analysis found the test split 91% clade-novel at J≥0.2, with recall flat out to ~20% nucleotide divergence (`docs/scientific-validity.md`).
- **A confounding check:** the composition-only baseline reaches AUPRC ≤0.16 against the CNN's ~0.60, so the model uses more than composition (Whalen pitfall 3).
- **Imbalance:** evaluation keeps the natural imbalance. Calls/genome is reported as cost rather than a balanced-accuracy figure (Whalen pitfall 5).
- **Seed replication** (three seeds per change) and an explicit statement that single-run ablations are inside the noise (`docs/model-reengineering.md`).
- **Two-layer accuracy:** agreement with geNomad is reported as agreement, not accuracy. No novelty claim is made (`docs/claims-and-evidence.md`).
- **Input-robustness checks:** reverse complement, N regions, grid phase, short contigs.
- **The permutation null:** it rotates calls within each contig, which preserves clumping, as regioneR's circular randomization does.
- **A dated decision trail** (decisions dated in each doc, and an internal session log). It makes the forking paths below auditable at all.
- **Job hygiene:** launcher preflights and exit-code propagation.

## Findings

Ordered by how much each could mislead a reader of the release numbers.

### F1. No evaluation sample is untouched (serious)

*Guidance:* Kapoor L1.1 (in its adaptive form), Cawley & Talbot, Dwork. A sample used to choose settings gives an optimistically biased estimate of performance under those settings.

*What happened:*

| sample | decisions it informed |
|---|---|
| 18 genomes / 52 elements | Chose the coverage builder family (2026-09-17; `docs/cli-and-packaging.md`). The EDL933 ambiguity bug was found on it and fixed (2026-09-22). It was the confirmation sample for `frac`/`gap`/`bridge` (`docs/gap-filling.md`), which is legitimate. But on 2026-09-29 its results were in view when refinement, the seed vote and the Conservative option were decided. |
| 42 genomes / 179 elements | Chose `frac` 0.1, `gap` 50 kb and `bridge` 0.1; the mode thresholds; the seed-vote frontier; the short-contig reference load; the ambiguity dose conditions. |
| Literature pilot, val genomes | Selected the training checkpoints (early stopping on the val split). |
| Test split in general | Both eval samples were drawn from it, and it has been consulted for about two weeks. It now behaves as a validation set. |

*Consequence:* every headline number is a tuning-sample number or a partly used confirmation number. The 0.961 (42) and 0.981 (18) figures should not be presented as held-out performance.

*Correction:*
- (a) The lockbox (`docs/lockbox.md`; scored once 2026-10-19, see "C9 result").
- (b) Label every reported number by the role its sample played: tuning, confirmation-after-use, or lockbox.
- (c) Report the lockbox-minus-42 gap as the measured optimism.

### F2. Selection among many variants, best reported (moderate)

*Guidance:* Cawley & Talbot; Weber ("bias due to focusing attention on particular methods should be avoided — for example, when tuning parameters").

*What happened:* dozens of builder, merge, threshold, refinement, vote and ramp variants were compared, and the best was adopted each time. The docs usually show the alternatives, which is good. But no analysis plan fixed the candidate set in advance, and the final pipeline is a composite of per-step winners.

*Correction:* covered by F1's lockbox. Going forward, write the candidate set and the decision rule down before a sweep is run.

### F3. Truth is the training labeler (serious, known)

*Guidance:* Kapoor L3 (the test distribution differs from the claim's distribution); Weber on gold standards; Philympics on the lack of curated truth for prophage callers.

*What happened:* geNomad labels both train and evaluation. Agreement is reported as agreement (good), but it is the only large-sample measure. The non-geNomad evidence is the literature pilot (57 elements, 7 genomes) and the spike-in (synthetic, one family). The lockbox fixes reuse, not circularity.

*Correction:* the LAMBDA clean-subset run (checklist §3) is the planned fix. Until it lands, no claim should read as accuracy.

### F4. Uncertainty ignores dependence and most variance sources (moderate)

*Guidance:* Whalen pitfall 2 (dependent examples); Bouthillier (variance sources beyond seed).

*What happened:* intervals, where given, are binomial over elements. Elements cluster within genomes: two adjacent misses in `GCF_001610875` and one genome driving the 18-genome deficit (`docs/model-quality-baseline.md`). Seed variance is shown as a three-seed range. Sample-to-sample variance is visible (42 vs 18 differ by 2 pts) but not quantified.

*Correction:* genome-level cluster bootstrap intervals for every headline recall and calls/genome, offline from local tracks. The lockbox endpoints already specify this.

### F5. Homology between training and evaluation elements is checked only at the nucleotide level (moderate)

*Guidance:* Whalen pitfall 2, Kapoor L3.2, DataSAIL.

*What happened:* the genome split is genus-level (one genome per genus). Phages cross host genera, though, and the clade analysis uses k=21 MinHash, which cannot see protein-level homology (`docs/scientific-validity.md`). The spike-in enforces 21-mer disjointness; its in-situ arm (31 loci) is not filtered.

*Correction:*
- Containment of every evaluation element (42, 18, lockbox, in-situ) against training positives, reported with recall stratified by containment.
- Amino-acid clustering remains post-0.1.0, as decided. The claim must say "nucleotide-level".

### F6. Self-assessment trap and curator exposure in the literature pilot (moderate)

*Guidance:* Weber (neutral benchmarking; "comparing one method with default parameters against another with extensive tuning makes for an unfair comparison").

*What happened:*
- **Tuning is asymmetric.** prophound's settings were tuned on geNomad-labeled data; geNomad ran at its defaults. That favors geNomad-agreement for us, which matters most where geNomad is the truth.
- **The agent built the literature truth.** The 11 `author_genes` rows are not hand-placed. `build_literature_truth.py` collects the genes whose author notes name the element (by regex), groups them into clusters split at gaps over 20 kb, and takes the largest cluster. *Corrected 2026-09-29: an earlier version of this finding called them hand-set after the EDL933 tracks were seen.* The risk is a rule that matches the wrong genes, not boundaries placed with the model's output in view.
- **The gains sit in the least-defined class.** Five of the six recovered elements are remnants or prophage-like elements.

*Correction:*
- List every judgment-call truth row.
- A sensitivity analysis with those rows excluded, and with remnants excluded.
- The user (or a second curator) reviews `scripts/literature_truth.tsv` before its numbers enter the release text.

### F7. The spike-in's implant sites are chosen by the model (minor to moderate)

*What happened:* implant sites require every window within 20 kb to score below 0.1 on the ensemble track (`docs/release-scope-0.1.0.md`). That conditions the experiment on the model: implants go only where the model is confident there is nothing. It is a reasonable way to attribute calls. But it selects atypical backgrounds, and the direction of the bias on detection is unknown.

*Correction:* a sensitivity arm with a model-independent site rule (element-free by geNomad plus distance only). It needs GPU re-scoring.

### F8. The off-target annotation's categories and controls are our own (minor)

*What happened:* the hallmark categories and their first-match order were designed for this analysis. geNomad's markers partly restate the labels (noted). Controls are drawn "clear of every call", which conditions them on the model not calling there, and so can make background look emptier than random sequence.

*Correction:* a second control drawn regardless of calls (CPU, cheap), and a fixed category definition written down before any rerun.

### F9. The recall definition was chosen after seeing numbers (minor)

*What happened:* covered50 was adopted on 2026-09-24 with results in view. It is the stricter option (0.961 against 0.978), so it does not inflate, but the 50% cut is arbitrary.

*Correction:* report covered at 25 / 50 / 75% alongside, offline and cheap.

### F10. The test distribution differs from the target inputs (moderate, known)

*Guidance:* Kapoor L3.3.

*What happened:* the corpus, both eval samples and the lockbox are all complete RefSeq genomes. The stated target inputs are drafts and MAGs, tested only by shredding complete genomes.

*Correction:* the fragmented-input comparison (G5) and a statement in the limitations section. A real-MAG sample would be stronger, but its truth would again be geNomad's.

### F11. Some reported numbers cannot be reproduced from the repo (moderate)

*Guidance:* Weber on reporting; DOME; Kapoor's computational reproducibility.

*What happened:*
- **Scripts:** `docs/error-analysis.md` says its scripts "are in the session scratchpad", and the downstream stage-1 measurements were "run from the session scratchpad" (`docs/downstream-postprocessing.md`). Those sessions' scratchpads are temporary.
- **Versions:** the geNomad version and database version behind every label are not recorded anywhere (`notes/DATA_PROVENANCE.md` lists this as unresolved).
- **Weights:** no checkpoint checksums are recorded.

*Correction:*
- Re-derive the error-analysis tables with checked-in scripts. `scripts/seed_vote_confirm.py` already reproduces the vote table exactly.
- Record `genomad --version`, the DB version and the checkpoint SHA-256s (one SSH command).
- From now on, analysis code that produces a reported number lives in `scripts/`.

### F12. Analysis designed and run by the same agent that proposes the method (process note)

The analyses, their interpretation and the method changes are all proposed by the same coding agent. The user makes every decision, which helps, but the framing of the options is the agent's. The lockbox rule, that the agent sees no model output on it before the freeze, addresses this for the final number only. Recommendation: for the release claims, the user reads the lockbox result and the literature truth directly, not only through the agent's summaries.

### F13. Early local runs did not record their device, and one device was wrong (moderate; found 2026-10-19)

*What happened:*
- MPS gave wrong scores for adaptive-max models (`docs/tech-debt-audit.md` S10).
- Local evaluations from 2026-09-17 to 09-20 ran on this Mac, at least one explicitly on MPS, and their tracks were not kept.
- Every run kept in the repo records its device, and all 1,769 are CUDA.

*Consequence:* design choices drawn from those early tables (the `pool` architecture, the first coverage-builder candidates) rest on possibly wrong scores. The held-out figures do not: the lockbox, LAMBDA and the 42/18 tracks are all CUDA.

*Correction:*
- the MPS fix, verified against CPU (30/30 reference cases);
- `docs/benchmarks.md` cross-device agreement;
- going forward, any evaluation used for a decision keeps its `run_metadata.json`, which records the device.

## Proposed corrective analyses (awaiting user review; none run)

| # | analysis | addresses | cost |
|---|---|---|---|
| C1 | Label every reported number by sample role; mark the 18 as "confirmation after use" | F1 | docs only |
| C2 | Genome-level cluster bootstrap CIs for headline recall and calls/genome (42, 18, phage-free, literature) | F4 | local, offline |
| C3 | Covered at 25/50/75% for the shipped modes | F9 | local, offline |
| C4 | Literature-truth audit: list judgment-call rows; sensitivity excluding them and excluding remnants; user review of the truth file | F6 | local + user time |
| C5 | Element containment vs training positives for all eval sets; recall stratified by containment | F5 | ARC CPU |
| C6 | Off-target annotation with an unconditioned control | F8 | ARC CPU |
| C7 | Spike-in with a model-independent site rule | F7 | ARC GPU |
| C8 | Reproducibility: re-derive the scratchpad-only tables with repo scripts; record geNomad/DB versions and checkpoint checksums | F11 | local + one SSH |
| C9 | Lockbox scored once at freeze, with the lockbox-minus-42 optimism gap | F1, F2 | ARC, after freeze |

C1–C3 and C8 are cheap and would change how every release number is worded. C4 matters most for the complement claim. C9 is the only one that yields a held-out number.

## Results of the corrective analyses run 2026-09-29 (C1, C2, C3, C4, C8)

Script: `scripts/integrity_checks.py` (`recall`, `lit`, `errortables`). Outputs are in `analysis/integrity/`. Shipped modes throughout (threshold + seed vote), 10,000 genome resamples.

### C1 — sample roles

Every reported number carries one of these roles. Release text uses the role's wording.

| sample | role | wording |
|---|---|---|
| 42 genomes / 179 elements | **tuning**: settings, thresholds, vote chosen on it | "on the tuning sample" |
| 18 genomes / 52 elements | **confirmation after use**: picked the builder family, surfaced the ambiguity bug, in view for the 2026-09-29 decisions | "on a second, previously used sample" |
| 20 geNomad-negative genomes | confirmation, first used 2026-09-28; in view for the vote decision | "on genomes with no geNomad provirus" |
| literature pilot, test genomes | confirmation (EDL933 was studied closely in 2026-09-22's ambiguity work) | "literature pilot" |
| literature pilot, val genomes | checkpoint selection (early stopping) | same, with the caveat |
| lockbox | **held out**: preregistered, scored once | the only "held-out" figure |

**Applied 2026-09-29, complete:**
- The README status, the CLI mode table and the checklist claim items.
- A standard sample-roles banner at the top of the 13 docs that report sample numbers (model-quality-baseline, scientific-validity, downstream-postprocessing, model-reengineering, release-scope-0.1.0, adaptive-multiscale-detection, error-analysis, multi-model-architecture, gap-filling, claims-and-evidence, final-benchmark, boundary-refinement, related-tools).
- 29 targeted rewordings wherever the text called a sample "held-out" or "independent" without qualification. The release-claim wording in `docs/claims-and-evidence.md` now says "held out from training but used for tuning" and carries the final-mode figure.
- Historical results keep their original numbers, marked "held out at the time" where they were.

### C2 + C3 — intervals and coverage cut-offs

**Final modes (recomputed 2026-09-29 after the mode change):**

| sample, mode | calls/genome | ≥1 bp | ≥25% | **≥50%** (genome bootstrap; Wilson) | ≥75% |
|---|---|---|---|---|---|
| 42, Exploratory 0.90 | 35.6 (31.3–40.1) | 1.000 | 1.000 | **0.994** (0.982–1.000; 0.969–0.999) | 0.989 |
| 42, Default 0.95 / 3 of 3 | 14.1 (12.1–16.3) | 0.978 | 0.978 | **0.972** (0.949–0.994; 0.936–0.988) | 0.967 |
| 42, Conservative 0.995 / 3 of 3 | 8.5 (7.3–9.6) | 0.961 | 0.950 | **0.944** (0.909–0.977; 0.900–0.969) | 0.916 |
| 18, Exploratory | 39.6 (29.9–50.2) | 1.000 | 1.000 | **1.000** (Wilson 0.931–1.000) | 0.981 |
| 18, Default | 13.3 (9.7–17.3) | 1.000 | 1.000 | **1.000** (Wilson 0.931–1.000) | 0.962 |
| 18, Conservative | 6.4 (4.7–8.5) | 0.981 | 0.981 | **0.962** (0.881–1.000; 0.870–0.989) | 0.923 |
| phage-free, E / D / C | 27.1 (19.1–35.5) / 7.4 (4.7–10.2) / 2.45 (1.35–3.7) | — | — | — | — |

Default's lower bound on the 42 is 0.949, just under 95%, so the reading below still holds. The table after this one is the first pass, under the interim modes.

**First pass (interim modes, Default 0.98):**

Ensemble covered recall at four coverage cut-offs. Intervals are 95% genome bootstrap (Wilson in brackets where it differs notably).

| sample, mode | calls/genome | ≥1 bp | ≥25% | **≥50%** | ≥75% |
|---|---|---|---|---|---|
| 42, Exploratory | 29.1 (25.6–32.7) | 0.994 | 0.989 | **0.983 (0.965–1.000)** | 0.978 |
| 42, Default | 11.3 (9.7–13.0) | 0.972 | 0.966 | **0.961 (0.935–0.984)** [Wilson 0.922–0.981] | 0.950 |
| 42, Conservative | 8.5 (7.3–9.6) | 0.961 | 0.950 | **0.944 (0.909–0.977)** | 0.916 |
| 18, Exploratory | 31.2 (23.3–39.7) | 1.000 | 1.000 | **1.000** [Wilson 0.931–1.000] | 0.962 |
| 18, Default | 10.5 (7.7–13.8) | 1.000 | 1.000 | **0.981 (0.922–1.000)** | 0.942 |
| 18, Conservative | 6.4 (4.7–8.5) | 0.981 | 0.981 | **0.962 (0.881–1.000)** | 0.923 |
| phage-free, Exploratory / Default / Conservative | 20.6 (14.6–26.7) / 4.9 (3.0–7.0) / 2.45 (1.35–3.7) | — | — | — | — |

- **No sample establishes ≥95%.** Default's lower bound is 0.935 on the 42 and 0.922 on the 18. The ≥95% criterion holds for the point estimate only. The release must say so, and the lockbox is the only chance of a tighter, held-out interval. **The lockbox (2026-10-19) does not establish it either:** Default 0.960, lower bound 0.924 (Wilson) / 0.929 (bootstrap).
- **The mode order is robust to the cut-off**: every cut-off orders Exploratory > Default > Conservative. The ≥75% cut costs 1–3 points relative to ≥50%, most for Conservative. Covered50 is not a flattering choice among nearby cut-offs.
- **Genome-level intervals are not always wider than Wilson.** Misses are spread thinly, and most genomes have none. The genome bootstrap is still the one to report, because it is the right unit.
- Caveat: with 18 or 20 genomes, the bootstrap intervals are rough, and a 1.000 estimate yields a degenerate [1.000, 1.000] interval. Report the Wilson interval there.

### C4 — literature-truth sensitivity (pre-vote modes; shipped modes follow from `822053`)

Coordinate classes of the 57 truth rows:
- **author_coordinates (41):** authors' own features, or phage genomes matched exactly by k-mers.
- **gene_span_judgement (12):** the extent was chosen by the agent from the author-assigned genes. This covers 5 EDL933 cryptic prophages, 6 *B. subtilis* elements and CP-933P's O-island.
- **locus_tag_span (4):** extents from literature locus-tag ranges.

| subset | n | ensemble Default | geNomad | paired difference (genome bootstrap) |
|---|---|---|---|---|
| all | 57 | 0.912 | 0.807 | +6 / −0; +0.105 (0.000 to +0.197) |
| author coordinates only | 41 | 0.927 | 0.829 | +4 / −0; +0.098 (0.000 to +0.158) |
| no gene-span judgment | 45 | 0.933 | 0.844 | +4 / −0; +0.089 (0.000 to +0.148) |
| no remnants | 43 | 0.977 | 0.953 | +1 / −0; +0.023 (0.000 to +0.049) |
| named prophages, author coordinates | 25 | 1.000 | 0.960 | +1 / −0; +0.040 (0.000 to +0.053) |

- **The ensemble never loses an element geNomad covers**, in any subset.
- **The advantage lives in remnants.** Remove them and it shrinks from 6 elements to 1 (Sakai Sp16). Two of the six gains (*B. subtilis* prophages 3 and 6) rest on gene-span judgment rows.
- **At Default, the genome-level interval includes zero in every subset.** The gains come from 2 of the 7 genomes, so "better than geNomad on curated truth" is not established by this pilot. Exploratory's all-rows difference (+8 / −0) is the only one whose interval excludes zero (+0.029 to +0.214).
- **Verified 2026-09-29 against the raw AL009126 record** (NCBI efetch):
  - The two gains on gene-span rows reproduce the authors' own annotation exactly.
  - *B. subtilis* prophage 3: the "defective prophage 3" / "prophage region 3" genes run *ydzTc*…*ydjJ*, 651,977–675,857, identical to the truth row. The region includes the *bsuM*/*bsuR* restriction-modification genes.
  - *B. subtilis* prophage 6: the "(putative) defective prophage 6" genes run *yoaM*…*yobO*, 2,036,039–2,078,626, identical. It carries an Xre/cro repressor, phage proteins and a pre-neck appendage protein.
  - The extraction is sound for these two. What remains a judgment is **class**: whether "defective prophage" regions count as prophage truth at all.
- **User review of the 16 non-author-coordinate rows, 2026-09-29.** Of the 16, 12 are covered by both tools, 2 are missed by both, and only *B. subtilis* prophages 3 and 6 separate them (verified above). Decisions:
  - **Headline = class `prophage` only.** `cryptic_prophage`, `remnant` and the new `phage_derived` class are reported as separate rows.
  - **TLC and the R2/F2 pyocin move to `phage_derived`**: phage-derived, not prophages.
  - **CP-933P is re-derived with the gene-cluster rule**, replacing its whole O-island: 56 genes, 2,283,855–2,330,172 (46.3 kb, against the island's 58.0 kb).
  - `build_literature_truth.py` encodes all three. `scripts/literature_truth.tsv` was patched to its output: the build inputs are not in the repo, and AE005174 was re-fetched from NCBI to derive CP-933P.

  Pre-vote results under the new classes (CP-933P still at its old extent until the ARC re-evaluation):

  | class | n | ensemble Default | geNomad | paired |
  |---|---|---|---|---|
  | **prophage (headline)** | 28 | **28/28** | 27/28 | +1 / −0 (Sakai Sp16) |
  | cryptic_prophage | 15 | 14/15 | 14/15 | 0 |
  | remnant | 12 | 8/12 | 3/12 | +5 / −0 |
  | phage_derived | 2 | 2/2 | 2/2 | 0 |

  **Headline reading:** on literature-named prophages, both tools are near ceiling and the ensemble adds one. The larger difference is on remnants and defective prophages, reported as such. The shipped-mode numbers follow from re-running `literature_pilot.slurm` with the new truth: **`822071`**, submitted 2026-09-29. Pull `analysis/litpilot/results.{json,elements.csv}`.

### C8 — reproducibility

- `integrity_checks.py errortables` re-derives `docs/error-analysis.md`'s scratchpad-only tables **exactly**: seed support 152/2, 97/1, 190/165 (off/on-target at 1, 2 and 3 seeds); 7 Default misses with maximum window scores 0.988, 0.983, 1.000, 0.992, 0.991, 0.999 and 0.911. `seed_vote_confirm.py` already reproduces the vote table. The graze table in `docs/downstream-postprocessing.md` "Per-element dig" is superseded by C3's cut-off table.
- **Recorded 2026-09-29** (`analysis/integrity/provenance.txt` on ARC):
  - geNomad database **version 1.9** (`genomad_db/version.txt`).
  - Checkpoint SHA-256s:
    - `prophage_v2_pool_r1` `3259673877ced18e79a92227f898b2009f69cfe82b692c59cb355d74dc27655f`
    - `prophage_v2_pool_r2` `7a50701a0be0f9cca87f000b428341e67d55242f77807321b3237ee73f31e59c`
    - `prophage_v2_pool_r3` `541cd72f1508ca35e2829b13eaf6818727db9cc53e55399c2ad467bba1199484`
  - **geNomad software version 1.9.4**, read by the user on 2026-09-29 from the `mge_explorer` env (`genomad --version`), the env both the corpus labeling and the lockbox build use.
  - **The corpus labels were made with the same versions** (confirmed by the user 2026-09-29, who ran the original labeling). Training, evaluation and lockbox truth all come from geNomad 1.9.4 with DB 1.9, so no version difference separates them.

### C6 — model-independent off-target controls (`822067`, 2026-09-29)

The same 607 Default calls (pre-vote, `calls42.csv`) and the same categories. Only the controls change. Hallmark rate = any MGE category.

| controls drawn | control_off_3seed | control_off_1_2seed | control_on_target |
|---|---|---|---|
| clear of calls and elements (821371) | 0.119 | 0.083 | 0.149 |
| **clear of geNomad elements only** | **0.168** | **0.104** | 0.177 |
| fully random | 0.186 | 0.136 | 0.278 |

- **F8 confirmed.** Excluding called regions made the background look emptier than random sequence.
- **Against the model-independent background, unanimous off-target calls (52.1% hallmark) are enriched 3.1×, not 4.3×.** One- and two-seed off-target calls (20.1%) are enriched 1.9×, not 2.4×.
- The enrichment survives, smaller. The release uses the elements-only figure, and `docs/error-analysis.md` is corrected.

### C5 — first attempt failed (`822067`)

`literature_truth.tsv` carries both `genome_id` and `accession`, and the rename made a duplicate column. The bug is fixed, the lockbox is added as a fifth set (with a per-set genome dir), and the job is to be resubmitted. C6 in the same job was unaffected.

### C5 result (`822223`, 2026-09-30)

Rerun as step 6 of `final_modes_rescore.slurm`. The job was clean: exit 0 and an empty `.err`. `analysis/integrity/containment.tsv` holds 521 evaluation elements against 2,742 training-positive intervals. Containment is 21-mer containment of the element in training positives.

| set | n | median | ≥0.1 | ≥0.5 | ≥0.9 |
|---|---|---|---|---|---|
| sizestrat42 | 179 | 0.005 | 42 | 18 | 7 |
| confirm18 | 52 | 0.008 | 19 | 16 | 4 |
| litpilot | 57 | **0.723** | 44 | 34 | 6 |
| insitu | 31 | 0.001 | 5 | 2 | 0 |
| lockbox | 202 | 0.003 | 16 | 6 | 1 |

Recall by containment: ensemble, final modes (shipped seed vote), covered50. The 42 and the 18 are rebuilt offline with `mode_recheck.shipped_calls`, and the totals reproduce the published 0.994 / 0.972 / 0.944 and 1.000 / 1.000 / 0.962. The literature pilot comes from `litpilot/results_final.elements.csv`.

| set, containment | n | Exploratory | **Default** | Conservative | geNomad |
|---|---|---|---|---|---|
| 42, <0.1 | 137 | 0.993 | **0.964** | 0.934 | — |
| 42, 0.1–0.5 | 24 | 1.000 | 1.000 | 0.958 | — |
| 42, ≥0.5 | 18 | 1.000 | 1.000 | 1.000 | — |
| 18, <0.1 | 33 | 1.000 | 1.000 | 0.939 | — |
| 18, ≥0.1 | 19 | 1.000 | 1.000 | 1.000 | — |
| literature, <0.1 | 13 | 0.846 | **0.615** | 0.615 | 0.462 |
| literature, 0.1–0.5 | 10 | 1.000 | 1.000 | 0.900 | 0.600 |
| literature, ≥0.5 | 34 | 1.000 | 1.000 | 1.000 | 1.000 |

- **The ≥95% figure is not carried by training look-alikes.** Three-quarters of the 42 are essentially novel to training (<0.1). Default's five misses on the 42 are all in that bin, which still reads 0.964.
- **The literature pilot is the most contained set** (median 0.723). Its pooled 0.912 is flattered: 34 of 57 elements are ≥50% contained, and both tools find all 34. On the 13 novel ones, Default finds 8 and geNomad 6. **Every novel miss is a remnant.** All four novel prophage and cryptic-prophage elements (SPβ, Sp18, Sp2, CP-933I) are found by both tools; Default's +2 over geNomad here are the remnants SpLE1 and *B. subtilis* prophage 6. Default's misses are SpLE2, SpLE4, SpLE6 and *B. subtilis* prophages 1 and 3. Exploratory recovers all but SpLE6 and prophage 1.
- **In-situ loci are not near-duplicates** (median 0.001; 2 of 31 at ≥0.5). This settles the open question in `docs/release-scope-0.1.0.md` A.1: the in-situ/implant gap is not a containment artifact.
- **The lockbox is nearly novel by construction** (median 0.003; 6 of 202 at ≥0.5), as the genus-disjoint build intended.
- The release reports the literature pilot split by containment, not only pooled.

### C7 — spike-in with model-independent implant sites (`822068`, 2026-09-29)

Same 241 hold-out-checked inoviruses, now all placed twice: 482 implants in 28 hosts. Sites use geNomad elements and spacing only, with no model scores. The eval ran under the final modes. `analysis/spikein_siterule/`.

| covered50 (implants) | model-chosen sites, 821410 | model-independent, all | **model-independent, no host call within 20 kb** |
|---|---|---|---|
| Default 0.95 / 3 of 3 | 0.320 (n=406) | 0.405 (n=482) | **0.362 (n=420)** |
| Default, 3–5 / 5–8 / 8–10 kb | 0.14 / 0.27 / 0.68 | 0.23 / 0.37 / 0.75 | **0.17 / 0.32 / 0.73** |
| Exploratory 0.90 | — (old: 0.95 only) | 0.674 | **0.617 (n=303)** |
| Conservative | 0.167 | 0.195 | — |

- **Model-independent sites carry the opposite bias.** An implant inside or near a region the unmodified host already calls gets credit from that call. Inside a host call, Default detects 0.833 (30 implants). Within 20 kb of one, 0.694 (62).
- The clean figure (no host call within 20 kb, host calls from `analysis/multiscale/*_3000` under the same mode) is the one to report.
- **F7's direction:** the model-conditioned rule *understated* small-element detection a little (0.320 against 0.362). It excluded sites with moderately elevated host background, which seems to help an implant get called.
- In-situ loci are unchanged by construction (same tracks). ~~The Exploratory comparison against the old sites waits for `822223`'s `results_final`.~~ **`822223`:** Exploratory 0.90 at the model-chosen sites detects 0.584 (3–5 / 5–8 / 8–10 kb: 0.18 / 0.55 / 0.92). That is against 0.674 raw and 0.617 clean at model-independent sites, the same direction as Default. The old rule understated detection in both modes.

### C4 final — literature pilot, final modes, revised classes (`822071`)

| class | n | Exploratory | Default | Conservative | geNomad |
|---|---|---|---|---|---|
| **prophage (headline)** | 28 | 28 | **28** | 28 | 27 |
| cryptic_prophage | 15 | 15 | 15 | 14 | 14 |
| remnant | 12 | 10 | 7 | 7 | 3 |
| phage_derived | 2 | 2 | 2 | 2 | 2 |
| all (covered50) | 57 | 0.965 | **0.912** | 0.895 | 0.807 |

- Paired against geNomad: Default **+6 / −0**, +0.105 (interval 0.029–0.151, now excluding zero); Exploratory +9 / −0 (0.029–0.240).
- The gains are Sp16, CP-933L (6 kb, missed by every earlier setting), SpLE1, SpLE3, SpLE5 and *B. subtilis* prophage 6. At Default, *B. subtilis* prophage 3 is no longer among them; at Exploratory it is. CP-933P at its gene-cluster extent is fully covered.
- The headline class is at ceiling for both tools. The difference is in remnants and cryptic prophages, reported as such.

### Next

C5, C6 and C7 go to ARC now that the cheap fixes are in (user, 2026-09-29). Built the same day:
- **C5:** `scripts/eval_containment.py`. Tested on a planted fixture: copy 1.0, novel 0.0, half 0.4995.
- **C6:** `offtarget_annotation.py --control_busy {elements_only,none}`.
- **C7:** `small_element_spikein.py build --site_rule genomad`.
- Launchers: `scripts/integrity_cpu.slurm` (C5 + C6), submitted as **`822067`**, and `scripts/spikein_siterule.slurm` (C7, GPU), submitted as **`822068`**.

Both new options default to the old behavior, so earlier results still reproduce. C8's version and checksum record goes in the same SSH call as the submission.

### C9 result — lockbox (`841208`, scored once 2026-10-19)

Default covered50 **0.960** (194/202; Wilson 0.924–0.980, genome bootstrap 0.929–0.984) at 11.0 calls/genome. Exploratory 0.965 / 32.4, and Conservative 0.901 / 4.9. Full table and pre-scoring checks: `docs/lockbox.md` "Result".

- **F1 (no untouched sample) is now closed for the headline:** the release carries one held-out figure from a preregistered, scored-once sample.
- **Measured optimism (F1(c)).** Reused figure minus lockbox: Default +0.012 against the 42 and +0.040 against the 18; Exploratory +0.029 / +0.035; Conservative +0.043 / +0.061. For Default the tuning sample's optimism is small, inside either interval. It grows with stringency, so Conservative's tuning-sample figure should not be quoted without the lockbox's.
- **≥95% is not established held-out.** By the rule fixed beforehand, the release gives 0.960 (0.924–0.980) and not "≥95%".
- **What it does not fix.** Truth is still geNomad (F3, circularity), and the lockbox holds 3 elements under 10 kb.
