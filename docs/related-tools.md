# Related tools and positioning

Living doc, updated in place. Built 2026-09-28 from a web survey, with Jaeger's prophage code read directly from source (v1.26.5, released 2026-09-04). It is the source for the related-work and positioning parts of the planned GitHub write-up (`docs/release-checklist.md` §4). No journal publication is planned, so this repo write-up is the only place the positioning is argued.

## Positioning (user direction, 2026-09-28)

The tool is presented as **a new tool in its own class: novel in concept, even though its methods resemble existing tools.** Wording rule (user, 2026-09-28): **never "guarantee"**. Recall on known elements is a *target that is measured and reported*, not a promise. The complementarity claim is "may not call", never "do not call".

**Positioning line** (also in the checklist):

> *mge-detect is a recall-first prophage screen. It localizes candidate integrated elements directly from nucleotide sequence, with no gene calling or reference database at run time. It is tuned to recover already-known elements and returns everything else as candidates, at a follow-up cost the user chooses. It is built to run alongside annotation-based tools and to surface what they may not call.*

### What is not new

- **The method.** Neural networks on raw nucleotides (Seeker, Jaeger, the LAMBDA gLMs), database-free ML (PhageBoost, PIDE) and sliding-window genome scans (LAMBDA, Jaeger, DRAD, DBSCAN-SWA) all exist.
- **Being database-free, and the word "discovery".** PhageBoost's paper is titled *Rapid discovery of novel prophages*. Prophage Hunter offers a similarity-free option explicitly to find novel phages.
- **Being recall-leaning, and having modes.** Philympics characterizes VirSorter, VirSorter2 and DBSCAN-SWA as trading precision for recall. geNomad ships `--conservative` and `--relaxed` presets for its post-classification filters. What is distinct here is stating the recall target explicitly and pricing each mode in follow-up cost (calls/genome), not having a recall-leaning mode at all.
- **The idea that unannotated calls can be real phages.** LAMBDA reviewed 305 genome-wide calls that fell outside its curated truth and judged 22 of them "likely phage". This appears there as a *finding* about a benchmark, not as a tool's design.

### What is distinct, as far as this survey found

1. **The contract.** Every tool surveyed is tuned so that its calls are trustworthy annotations. None sets a recall target on known elements, prices its modes as follow-up cost (calls/genome), or presents itself as a candidate generator meant to run *alongside* annotation tools. This is the main basis for the "own class" claim.
2. **It learns integrated elements in their host context.** Jaeger masks prophages out of its training genomes with N, so it learns free phage versus host. LAMBDA and Seeker train on free phage genomes versus bacterial genomes. mge-detect trains on prophages *as they sit in the chromosome*, with host sequence within 5 kb of each element as hard negatives. It learns the actual task, which is telling an integrated element from its own flanking host. The honest counterpart belongs next to it: those tools' positive labels are independent of any prophage caller, while ours come from geNomad.
3. **Boundaries come from the score track alone, with no gene calling anywhere in the pipeline.** Jaeger's prophage mode always runs gene prediction (pyrodigal-gv) to snap boundaries to intergenic positions, and it uses tRNAs and integrases when an annotation is supplied. geNomad, PhageBoost and PIDE all start from gene calls. LAMBDA's pipeline is gene-free, but it is not a packaged tool (see 4).
4. **It is a packaged tool that works on fragmented input.** Jaeger's prophage mode only scans contigs longer than **500 kb** by default (`--lc 500000`), so it is aimed at near-complete genomes. mge-detect scans contigs down to 3 kb with calibrated per-scale thresholds (`--short-contigs scan`). LAMBDA is a SLURM research pipeline: its fine-tuned checkpoints are available on request, and the released ones need the 7 B-parameter Evo2 base model.

Items 2–4 are facts about design, and each should be stated as such. Item 1 is the concept. What the contract delivers in practice rests on the evidence in checklist §3 ("Exploratory-claim evidence").

## The landscape

| tool | model input | gene calling? | reference DB at run time? | how regions and boundaries are set | how it relates |
|---|---|---|---|---|---|
| **geNomad** `find-proviruses` | genes | yes | yes (marker profiles; MMseqs2 for integrases) | CRF over marker-gene specificity; edges extended to tRNAs (ARAGORN) and integrases | This tool's label source. Annotation-based. Tuned for precise calls. |
| **Jaeger** (v1.26.5, 2026) | 2 kb nucleotide windows, six-frame translated, dilated CNN | **yes, for boundaries** (pyrodigal-gv); an annotation is optional | no (homology-free) | See "Jaeger in detail" | Nearest *method* neighbour. It localizes prophages, but only on contigs >500 kb by default and with gene-snapped edges. It trained on free phage versus host, with prophages masked. |
| **LAMBDA** (2026 benchmark) | 2/4/8 kb windows, 9 genomic LMs (1 M–7 B params) | no | no | See "LAMBDA in detail" | **Closest pipeline**, but a research benchmark, not a tool. Its curated truth set is useful to us (below). |
| **PhageBoost** | genes | yes (Pyrodigal) | **no** (gene features + XGBoost) | per-gene probability >0.9, ≥10 genes, allowed gaps, Parzen smoothing | **Closest in framing**: database-free, "discovery of novel prophages". Labels came from pVOG vs COG gene stretches. Needs gene calls, and has no recall target. |
| **PIDE** (2025) | proteins | yes | no homology search; fine-tuned ESM-2 (650 M) | predicted phage ORFs clustered by gene density | Database-free, protein-level ML. Handles MAGs. |
| **PhageTAILor** (2026 preprint) | genes | yes | yes (6 homology detectors) + LightGBM | homology candidates, then classified | Adds classes of tail-like elements. Homology-based. |
| **PHI** (2025 preprint) | — | — | — | Galaxy workflow that wraps existing detectors | An integration layer, not a new detector. |
| **Seeker** | sequence | no | no (LSTM) | whole sequence, phage vs bacterium | The authors state it was not trained or tested for prophages inside bacterial genomes. |
| **DeepVirFinder / PPR-Meta** | contigs | no | no | per-contig viral vs not | Contig classifiers that don't localize. §3's comparison needs caveats for this. |
| **DEPhT** | genes | yes | partial (architecture features to find regions; targeted homology to extract) | shell-genome trimming | Discovery-oriented. Trained per genus (Mycobacterium). Degrades on assemblies with more than 100 contigs. |
| **PhiSpy, DRAD, DBSCAN-SWA** | mixed | mostly yes | mixed | composition or sliding windows | Predecessors. Composition-based detection is decades old. |
| **Prophage Hunter** | genes | yes | optional | regions plus an activity score | Its similarity-free option exists "to find novel phages". |

### Jaeger in detail (read from `postprocess/prophages.py` and `prophage_boundaries.py`)

- **Scoring.** Windows of 2 kb (`--fsize`), stride 2 kb by default, 4-class logits (phage / bacteria / eukarya / archaea). Prophage-containing training genomes had their prophages masked with N.
- **Segmentation** (`segment`), run only on contigs longer than `--lc` (default 500 kb):
  1. Seeds are windows with phage logit > `--sensitivity` (default 1.5).
  2. Each seed extends while the region's mean phage logit stays above its mean host logit.
  3. A kernel change-point detector (`ruptures.KernelCPD`, penalty chosen by knee point) proposes segments whose mean exceeds the sensitivity.
  4. The two region sets are merged. Filters then keep regions of 6–250 kb with coefficient of variation <0.6 whose mean phage logit exceeds both the host mean and the sensitivity.
- **Boundaries** (`refine_prophage_boundaries`, always run):
  1. pyrodigal-gv predicts genes on the contig, and each end is moved outward to the nearest intergenic position (at most 2 × `fsize`).
  2. If a GenBank annotation is supplied, ends snap to the nearest integrase (matched by product name) or tRNA within 20 kb, provided the phage score still exceeds the host score there.
  3. An optional `--filter-low-gene-density` drops regions sparser in genes than the host (complete genomes only).
- **Accuracy.** No published prophage-level figures were found. The repo has a Casjens 2003 validation script, but its results are not in the docs.

**Consequence:** Jaeger localizes prophages from a nucleotide-input model, so "localizes from sequence" alone is not distinct. The distinctions that hold against Jaeger are items 2–4 above, plus the contract.

### LAMBDA in detail

- Windows of 2/4/8 kb with 50% overlap, collapsed by majority vote, then per-genome z-score, bidirectional exponentially weighted smoothing, a threshold, gap merging and a length filter. Hyperparameters are tuned per model. No gene calling. Boundaries are at window resolution.
- Training data is free phage genomes (INPHARED, clustered with vclust) versus GTDB bacteria, split by cluster and genus.
- Truth is **380 expert-curated regions in 80 genomes** (Casjens gold standard, PHASTER manual annotations, Gauthier et al.). geNomad scores MCC 0.794 / recall 87.3%. The best gLM, Evo2, scores MCC 0.722 / recall 71.8%, and ProkBERT-mini-long scores MCC 0.704. The paper concludes that homology tools still win.
- Code: [github.com/leannmlindsey/LAMBDA](https://github.com/leannmlindsey/LAMBDA). Data: Zenodo 10.5281/zenodo.19236552, including per-window labels for the 80 genomes and `prophage_reference_locations.csv`.

### Philympics (2021)

A community gold standard: 78 manually curated genomes with CDS-level `is_phage` tags, and a Snakemake framework that accepts new tools as rules (MIT licence). Published CDS-level results: PhiSpy has the best F1 (0.733; recall 0.731), VirSorter2 the highest recall (0.766, precision 0.399), and Seeker, the only nucleotide deep-learning tool included, F1 0.093. geNomad is not in it. Contributing mge-detect as a rule is a credible external route to publishing known accuracy. Genome overlap with our corpus, and the LAMBDA baseline computed from its published per-tool predictions, are in `docs/claims-and-evidence.md`.

## Consequences for the write-up

- **Don't say "first"**, and don't say "guarantee". Say the tool is *designed for a different job*, and cite the tools above as related work.
- **Lead with the contract, then the three design facts** (host-context training, gene-free boundaries, fragmented input). The method comes last.
- **State the label trade-off against LAMBDA and Jaeger openly.** Their labels are caller-independent, ours come from geNomad, and ours are in host context. Each choice buys something different.
- **Claims, benchmarks and the evidence each claim needs are in `docs/claims-and-evidence.md`**, including the two-layer point: 96% agreement with geNomad is not 96% accuracy, since geNomad itself covers about 81% of LAMBDA's curated prophages.
- **Run on LAMBDA's truth set.** Its 380 curated regions are not geNomad-derived, and published geNomad and gLM numbers exist on the same set, with per-window labels already on Zenodo. It is the cheapest comparison not bounded by geNomad. Check genome overlap with our training split first (*E. coli* K-12 at least is in training) and report the non-overlapping subset. Its 22 reviewed "likely phage" off-truth regions also give a small, external check on candidate calls.
- **Recall figures are not comparable across truth sets or matching rules.** Our 0.961 covered50 is measured against geNomad on held-out genomes. Only a run on LAMBDA's own truth set, with its matching rule, supports a side-by-side comparison.

## Sources

- geNomad: [pipeline docs](https://portal.nersc.gov/genomad/pipeline.html), [provirus identification](https://portal.nersc.gov/genomad/provirus_identification.html), [Nat Biotechnol 2023](https://www.nature.com/articles/s41587-023-01953-y)
- Jaeger: [bioRxiv 2024](https://www.biorxiv.org/content/10.1101/2024.09.24.612722v1), [GitHub](https://github.com/Yasas1994/Jaeger) (source read: `src/jaeger/postprocess/prophages.py`, `prophage_boundaries.py`, `commands/predict.py`, `docs/_source/usage.md`), [docs](https://jaeger.readthedocs.io/en/latest/)
- LAMBDA: [NAR Genom Bioinform 2026](https://academic.oup.com/nargab/article/8/3/lqag103/8787816), [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC13041943/), [GitHub](https://github.com/leannmlindsey/LAMBDA), [interactive benchmark](https://leannmlindsey.github.io/lambda-benchmark/)
- PhageBoost: [NAR Genom Bioinform 2021](https://academic.oup.com/nargab/article/3/1/lqaa109/6066536), [GitHub](https://github.com/ku-cbd/PhageBoost)
- PIDE: [Genome Biology 2025](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-025-03733-0)
- PhageTAILor: [bioRxiv 2026](https://www.biorxiv.org/content/10.64898/2026.08.24.746745v1.full); PHI: [bioRxiv 2025](https://www.biorxiv.org/content/10.64898/2025.12.02.691814v1)
- Seeker: [NAR 2020](https://academic.oup.com/nar/article/48/21/e121/5921300)
- DEPhT: [NAR 2022](https://academic.oup.com/nar/article/50/13/e75/6572362)
- PhiSpy: [NAR 2012](https://academic.oup.com/nar/article/40/16/e126/1027055); DRAD: [PLOS One](https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0001193); DBSCAN-SWA: [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9061938/)
- Prophage Hunter: [NAR 2019](https://academic.oup.com/nar/article/47/W1/W74/5494712)
- Philympics: [paper](https://f1000research.com/articles/10-758), [framework](https://github.com/linsalrob/ProphagePredictionComparisons)
- geNomad presets: [post-classification filtering](https://portal.nersc.gov/genomad/post_classification_filtering.html)
