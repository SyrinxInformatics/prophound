# Research integrity

How prophound's evaluation guards against the usual ways a machine-learning method in biology overstates itself: tuning on the data it reports, leakage between training and test, circular truth and unreported variance. The project's methods were audited against published guidance. This doc gives each finding, what was done about it and what it showed, and what remains open.

## The guidance

- Weber et al. 2019, *Essential guidelines for computational method benchmarking*, Genome Biology 20:125 ([link](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1738-8)): neutrality and the self-assessment trap, fair tuning, gold standards, several metrics, reporting.
- Whalen, Schreiber, Noble and Pollard 2022, *Navigating the pitfalls of applying machine learning in genomics*, Nat Rev Genet ([link](https://www.nature.com/articles/s41576-021-00434-9)): distributional differences, dependent examples, confounding, leaky preprocessing and unbalanced classes; group splits over random ones.
- Kapoor and Narayanan 2023, *Leakage and the reproducibility crisis in machine-learning-based science*, Patterns ([arXiv](https://arxiv.org/abs/2207.07048)): eight types of leakage, from having no real test set to sampling bias in the test distribution.
- Walsh et al. 2021, *DOME: recommendations for supervised machine learning validation in biology*, Nat Methods ([link](https://www.nature.com/articles/s41592-021-01205-4)).
- Cawley and Talbot 2010, *On over-fitting in model selection and subsequent selection bias in performance evaluation*, JMLR ([link](https://www.jmlr.org/papers/v11/cawley10a.html)): choosing settings on the data you report is itself overfitting.
- Dwork et al. 2015, *The reusable holdout*, Science ([link](https://www.science.org/doi/10.1126/science.aaa9375)): consulting a holdout repeatedly erodes it.
- Bouthillier et al. 2021, *Accounting for variance in machine learning benchmarks*, MLSys ([link](https://proceedings.mlsys.org/paper_files/paper/2021/file/0184b0cd3cfb185989f858a1d9f5c1eb-Paper.pdf)).
- Joeres, Blumenthal and Kalinina 2025, *Data splitting to avoid information leakage with DataSAIL*, Nat Commun ([link](https://www.nature.com/articles/s41467-025-58606-8)).
- Phage-tool benchmarking: Roach et al. 2022, *Philympics* ([link](https://f1000research.com/articles/10-758)); Schackart et al. 2023 ([link](https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2023.1078760/full)); Ho et al. 2023 ([link](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-023-01533-x)).
- Genomic interval statistics: regioneR (Gel et al. 2016, [link](https://academic.oup.com/bioinformatics/article/32/2/289/1744157)), on nulls that preserve region clumping.

## What every figure is labeled with: the role of its sample

The single largest risk the audit found was that no evaluation sample was untouched: the samples used to report performance had also been used to choose settings. So every figure in the docs carries the role its sample played.

| sample | role | what it was used for |
|---|---|---|
| 42 genomes, 179 elements (size-stratified) | **tuning** | the call builder's settings, merging, the mode thresholds, the seed vote |
| 18 genomes, 52 elements | **confirmation after use** | chose the call builder's family, surfaced an ambiguity-code bug, and was in view when the modes were set |
| 20 genomes with no geNomad element | confirmation | the call floor; in view for the seed vote |
| literature pilot (57 elements, 7 genomes) | confirmation | curated truth from the literature; its validation-split genomes also selected the training checkpoints |
| first held-out test (150 genomes) | **held out**, then superseded | the defaults before the boundary stages; preregistered and scored once |
| **second held-out test** (150 genomes, 186 elements) | **held out** | the shipped 0.1.0 pipeline; preregistered and scored once ([the held-out test](lockbox2.md)) |

"Held out" is used only for the last two. A figure from the 42 is "on the tuning sample", and from the 18 "on a second, previously used sample".

## Findings and what was done

| finding | risk | what was done | result |
|---|---|---|---|
| **No sample was untouched** | headline figures would be tuning figures | two preregistered, genus-disjoint test sets, each scored once; every figure labeled by role | the shipped pipeline covers 0.989 (0.962–0.997) held out; the tuning sample's optimism is measured ([below](#how-optimistic-were-the-tuning-figures)) |
| **The best of many variants was kept** | selection bias | covered by the held-out tests; candidate sets and decision rules written down before later sweeps | the held-out figures include any selection effect |
| **Truth is the labeler** (geNomad labels both training and evaluation) | agreement read as accuracy | agreement reported only as agreement; non-geNomad truth from LAMBDA's expert curation and the literature | geNomad covers about 82% of LAMBDA's curated prophages; prophound 0.933 ([Comparisons](tool-comparisons.md)) |
| **Intervals ignored dependence** (elements cluster in genomes) | intervals too narrow | genome-level bootstrap intervals beside Wilson intervals for every headline figure | [below](#intervals-and-coverage-cut-offs) |
| **Homology between training and test elements** was checked only by genus | look-alikes inflate recall | 21-mer containment of every evaluation element in the training positives; recall reported by containment | recall is not carried by look-alikes ([below](#are-the-test-elements-novel-to-training)) |
| **The literature truth was assembled by the analyst**, and tuning was asymmetric | self-assessment trap | judgment-call rows listed and reviewed; headline restricted to the "prophage" class; sensitivity analysis | the advantage over geNomad is in remnants, not named prophages ([below](#the-literature-truth)) |
| **Spike-in implant sites were chosen by the model** | conditioning on the model | a second arm with sites set by geNomad elements and spacing only | the model-chosen sites slightly understated detection ([below](#small-elements-spike-in-sites)) |
| **Off-target controls were drawn clear of calls** | background looks emptier than it is | a control drawn clear of geNomad elements only | enrichment 3.1×, not 4.3× ([below](#off-target-enrichment-controls)) |
| **The recall definition (covered50) was chosen with results in view** | a flattering cut-off | coverage reported at 1 bp, 25%, 50% and 75% | the mode order holds at every cut-off ([below](#intervals-and-coverage-cut-offs)) |
| **The test distribution differs from the target inputs** (complete genomes, not MAGs) | performance on drafts unknown | shredded-genome comparisons (G5, G5b) and a stated limitation | retention 0.92–0.94, no comparative claim ([Comparisons](tool-comparisons.md)) |
| **Some figures could not be reproduced from the repository** | unverifiable numbers | analysis code moved into `scripts/`; geNomad and database versions and checkpoint checksums recorded | [below](#provenance) |
| **The same analyst designed, ran and interpreted the analyses** | framing bias | preregistration with blinding for the held-out tests; the project owner reads the held-out results and the literature truth directly | the first held-out test fell short of its rule (lower bound 0.924) and was reported that way; two fragmented-input comparisons with bars set in advance (G5, G5b) did not clear and are reported that way |
| **Early local runs did not record their device, and one device was wrong** (Apple GPU pooling gave wrong scores) | decisions on wrong scores | the bug was fixed and verified against CPU; every run now records its device | every held-out, LAMBDA and tuning figure ran on CUDA and is unaffected |

### What already met the guidance

- **Group splits:** one genome per genus. On the test split, 91% of elements are clade-novel at Jaccard 0.2, and recall stays flat out to about 20% nucleotide divergence.
- **A confounding check:** a composition-only model reaches AUPRC 0.16 or less against the network's 0.60 on the same split, so the model uses more than composition.
- **Imbalance:** evaluation keeps the natural imbalance and reports candidates per genome as the cost.
- **Seed replication:** three seeds per training change. Single runs fall inside the spread between seeds.
- **The permutation null** rotates calls within each contig, which preserves their clumping, as regioneR's circular randomization does.
- **Input robustness:** reverse complement, runs of N, window-grid phase and short contigs are checked as a standing rule.

## Results of the corrective analyses

Scripts: `scripts/integrity_checks.py`, `scripts/eval_containment.py`, `scripts/offtarget_annotation.py`, `scripts/small_element_spikein.py`. Outputs: [`analysis/integrity/`](../prophound/models/prophage/analysis/integrity/), [`analysis/litpilot/`](../prophound/models/prophage/analysis/litpilot/), [`analysis/spikein_siterule/`](../prophound/models/prophage/analysis/spikein_siterule/). The figures in this section are for the modes as shipped (threshold plus seed vote), before the boundary stages were added; the boundary stages' effect is measured in the [held-out test](lockbox2.md).

### Intervals and coverage cut-offs

Covered recall at four cut-offs, with 95% genome-bootstrap intervals (10,000 resamples) and Wilson intervals:

| sample, mode | candidates per genome | ≥1 bp | ≥25% | **≥50%** (bootstrap; Wilson) | ≥75% |
|---|---|---|---|---|---|
| 42, Exploratory | 35.6 (31.3–40.1) | 1.000 | 1.000 | **0.994** (0.982–1.000; 0.969–0.999) | 0.989 |
| 42, Default | 14.1 (12.1–16.3) | 0.978 | 0.978 | **0.972** (0.949–0.994; 0.936–0.988) | 0.967 |
| 42, Conservative | 8.5 (7.3–9.6) | 0.961 | 0.950 | **0.944** (0.909–0.977; 0.900–0.969) | 0.916 |
| 18, Exploratory | 39.6 (29.9–50.2) | 1.000 | 1.000 | **1.000** (Wilson 0.931–1.000) | 0.981 |
| 18, Default | 13.3 (9.7–17.3) | 1.000 | 1.000 | **1.000** (Wilson 0.931–1.000) | 0.962 |
| 18, Conservative | 6.4 (4.7–8.5) | 0.981 | 0.981 | **0.962** (0.881–1.000; 0.870–0.989) | 0.923 |
| no-element genomes, E / D / C | 27.1 / 7.4 / 2.45 | | | | |

- **No development sample established 95%:** Default's lower bound is 0.949 on the 42. Only the held-out test could, and did ([the held-out test](lockbox2.md)).
- **The mode order is the same at every cut-off,** and 75% costs 1–3 points against 50%. Covered50 is not a flattering choice among nearby cut-offs.
- **Genome-level intervals are not always wider than Wilson's,** because misses are spread thinly. With 18 or 20 genomes, a 1.000 estimate gives a degenerate bootstrap interval, so the Wilson interval is reported there.

### How optimistic were the tuning figures?

Each sample's figure minus the held-out figure, covered50:

| | Exploratory | Default | Conservative |
|---|---|---|---|
| first held-out test: the 42 minus held out | +0.029 | +0.012 | +0.043 |
| first held-out test: the 18 minus held out | +0.035 | +0.040 | +0.061 |
| second held-out test: the 42 minus held out | +0.010 | −0.017 | +0.041 |
| second held-out test: the 18 minus held out | +0.016 | +0.011 | +0.059 |

For Default the tuning sample's optimism is small, inside either interval. It grows with stringency, so Conservative's tuning figures should not be quoted without the held-out one.

### Are the test elements novel to training?

21-mer containment of each evaluation element in the training positives (521 elements against 2,742 training intervals; [`analysis/integrity/containment.tsv`](../prophound/models/prophage/analysis/integrity/containment.tsv)):

| set | elements | median containment | ≥0.1 | ≥0.5 | ≥0.9 |
|---|---|---|---|---|---|
| tuning (42) | 179 | 0.005 | 42 | 18 | 7 |
| confirmation (18) | 52 | 0.008 | 19 | 16 | 4 |
| literature pilot | 57 | **0.723** | 44 | 34 | 6 |
| in-situ inovirus loci | 31 | 0.001 | 5 | 2 | 0 |
| first held-out test | 202 | 0.003 | 16 | 6 | 1 |

Recall by containment, covered50:

| set, containment | elements | Exploratory | **Default** | Conservative | geNomad |
|---|---|---|---|---|---|
| 42, under 0.1 | 137 | 0.993 | **0.964** | 0.934 | |
| 42, 0.1–0.5 | 24 | 1.000 | 1.000 | 0.958 | |
| 42, 0.5 and over | 18 | 1.000 | 1.000 | 1.000 | |
| 18, under 0.1 | 33 | 1.000 | 1.000 | 0.939 | |
| 18, 0.1 and over | 19 | 1.000 | 1.000 | 1.000 | |
| literature, under 0.1 | 13 | 0.846 | **0.615** | 0.615 | 0.462 |
| literature, 0.1–0.5 | 10 | 1.000 | 1.000 | 0.900 | 0.600 |
| literature, 0.5 and over | 34 | 1.000 | 1.000 | 1.000 | 1.000 |

- **The 95% target is not carried by training look-alikes.** Three-quarters of the tuning elements are essentially novel to training, and Default's misses there are all in that group, which still reads 0.964.
- **The literature pilot is the most contained set,** so its pooled figure is flattered and is reported split by containment. On its 13 novel elements Default finds 8 and geNomad 6; every novel miss is a remnant.
- **The first held-out set is nearly novel to training** (median 0.003), as its genus-disjoint construction intended. The second was built the same way; its containment was not measured.
- This is nucleotide-level containment. Protein-level homology, which k-mers cannot see, is untested, so claims say "nucleotide-level".

### The literature truth

The 57 truth rows: 41 use the authors' own coordinates, 12 take an element's extent from the genes the authors assign to it, and 4 come from locus-tag ranges in the literature. Every one of the 16 non-author rows was reviewed. The two that separate prophound from geNomad (*B. subtilis* prophages 3 and 6) were checked against the raw GenBank record and reproduce the authors' annotation exactly. Elements are reported by class, with "prophage" as the headline:

| class | elements | Exploratory | Default | Conservative | geNomad |
|---|---|---|---|---|---|
| **prophage** (headline) | 28 | 28 | **28** | 28 | 27 |
| cryptic prophage | 15 | 15 | 15 | 14 | 14 |
| remnant | 12 | 10 | 7 | 7 | 3 |
| phage-derived (tailocins, pyocins) | 2 | 2 | 2 | 2 | 2 |
| all (covered50) | 57 | 0.965 | **0.912** | 0.895 | 0.807 |

- Paired against geNomad, Default is +6 / −0 (+0.105, interval 0.029–0.151).
- **Both tools are at ceiling on named prophages;** the difference is in remnants and cryptic prophages, and is reported as such. With 7 genomes, this is an example, not a general claim.

### Small elements: spike-in sites

Inoviruses (3–10 kb, checked against the training data) implanted into host genomes, with sites chosen two ways:

| covered50 of implants | sites chosen by the model (background below 0.1) | sites from geNomad elements and spacing only, no host call within 20 kb |
|---|---|---|
| Default | 0.320 (406 implants) | **0.362** (420) |
| Default, 3–5 / 5–8 / 8–10 kb | 0.14 / 0.27 / 0.68 | 0.17 / 0.32 / 0.73 |
| Exploratory | 0.584 | **0.617** (303) |

- The model-chosen sites slightly understated detection, in both modes.
- Sites inside or near a host call get credit from that call (0.833 inside one), so the clean figure excludes them.
- These are synthetic insertions of one family, so small-element recall is reported with that caveat.

### Off-target enrichment: controls

The share of Default calls outside geNomad elements that carry a mobile-element feature, against random host stretches of the same length, drawn three ways:

| controls drawn | controls with a feature | unanimous off-target calls (52.1%) enriched |
|---|---|---|
| clear of calls and elements | 0.119 | 4.3× |
| **clear of geNomad elements only** | **0.168** | **3.1×** |
| fully random | 0.186 | 2.8× |

Excluding called regions made the background look emptier than random sequence. The release uses the elements-only control, so the enrichment is 3.1×.

### Provenance

- **Labels:** geNomad 1.9.4 with database 1.9 for the training labels, every evaluation sample and both held-out tests, so no version difference separates them.
- **Training checkpoints** (SHA-256):
  - `prophage_v2_pool_r1` `3259673877ced18e79a92227f898b2009f69cfe82b692c59cb355d74dc27655f`
  - `prophage_v2_pool_r2` `7a50701a0be0f9cca87f000b428341e67d55242f77807321b3237ee73f31e59c`
  - `prophage_v2_pool_r3` `541cd72f1508ca35e2829b13eaf6818727db9cc53e55399c2ad467bba1199484`

  The bundled weights-only files hold the same tensors and give byte-identical outputs; their own hashes are in `prophound/models/prophage/model/manifest.json`.
- **Devices:** every run behind a reported figure records its device in `run_metadata.json`.

## What remains open

None of these is in scope for 0.1.0. The release makes no claim that depends on them, and the [release text](release-text-0.1.0.md) states each one as a limitation or a non-claim. Each needs new data or new labels, so each is planned for after 0.1.0.

- **Circular truth.** The large-sample figures are agreement with geNomad, which also labeled the training data. The curated benchmarks (LAMBDA, 39 genomes; the literature pilot, 7) are small. Closing it needs a larger literature-curated benchmark across taxa, scored with the permutation test, and training labels that do not come from one tool.
- **Novelty.** No test measures whether prophound finds elements that no database contains, so 0.1.0 claims no novel elements. One computational test is possible: elements that a newer database release adds but an older one missed, called by a model trained only on the older labels. Whether geNomad's database versions allow this has not been checked.
- **Mechanism.** Composition alone does not explain the model's scores, but which sequence features it uses is unmeasured. The planned test asks whether high-scoring windows are enriched for phage gene classes against matched background.
- **Protein-level homology** between training and test elements. Containment is nucleotide-level, so novelty-to-training claims say "nucleotide-level". Amino-acid-level clustering would close it.
- **Real draft assemblies and MAGs,** as opposed to shredded complete genomes. This needs real MAGs with matched complete genomes.
- **Small elements** under 10 kb: thin in every held-out set, because geNomad truth holds few of them, and measured mainly on synthetic implants of one family. This needs a second source of real small elements, which the planned retraining also uses.

## Lab validation

No prophound candidate has been confirmed in the lab, and lab work is not part of the project's plans. Every figure in these docs is computational. A candidate is a region worth checking, and checking it is up to the researcher using the tool ([how to read the results](interpreting-results.md#checking-a-candidate-in-the-lab)).
