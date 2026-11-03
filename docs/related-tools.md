# Related tools

How prophound relates to other prophage and phage detectors: what it shares with them, what is distinct, and where its parts come from. Head-to-head results are in [Comparisons with established tools](tool-comparisons.md).

**How this was surveyed.** A web survey of published and released tools, repeated once after the boundary stages were added, plus a direct read of Jaeger's prophage code (v1.26.5: `postprocess/prophages.py`, `prophage_boundaries.py`). It rests on what each tool documents, not on a full read of every paper. "Not found" means not found in this survey, not that it does not exist.

## Where prophound sits

> prophound is a recall-first prophage screen. It localizes candidate integrated elements directly from nucleotide sequence, with no gene calling or reference database at run time. It is tuned to recover already-known elements and returns everything else as candidates, at a follow-up cost the user chooses. It is built to run alongside annotation-based tools and to surface what they may not call.

## What is shared with existing tools

- **The method.** Neural networks on raw nucleotides (Seeker, Jaeger, the genomic language models in LAMBDA), database-free machine learning (PhageBoost, PIDE) and sliding-window genome scans (LAMBDA, Jaeger, DRAD, DBSCAN-SWA) all exist.
- **Being database-free, and aiming at discovery.** PhageBoost's paper is titled *Rapid discovery of novel prophages*. Prophage Hunter offers a similarity-free option explicitly to find novel phages.
- **Leaning toward recall, and offering modes.** Philympics describes VirSorter, VirSorter2 and DBSCAN-SWA as trading precision for recall. geNomad ships `--conservative` and `--relaxed` presets, and its developer writes that "by default, precision is preferred over recall. If you want to maximize discovery, use the `--relaxed` parameter" (issue #110). Jaeger's prophage extraction has a `--sensitivity` setting (0–4, default 1.5).
- **The idea that unannotated calls can be real phages.** LAMBDA reviewed 305 genome-wide calls outside its curated truth and judged 22 of them "likely phage", as a finding about its benchmark.

## What is distinct

Among the database-free, sequence-level callers (PhageBoost, PIDE, Jaeger, Seeker and LAMBDA's models), each ingredient alone has a precedent:

| ingredient | precedent |
|---|---|
| discovery framing | PhageBoost; Jaeger (phages "hidden within metagenomes"); Prophage Hunter's similarity-free option; geNomad's `--relaxed` |
| adjustable sensitivity | Jaeger `--sensitivity`; geNomad presets; DBSCAN-SWA's cluster parameters |
| recall over precision | DBSCAN-SWA (described in its preprint as "highly nonconservative"; not confirmed in the published version); LAMBDA's 22 "likely phage" regions, as a benchmark finding |
| accuracy framing (the opposite) | PIDE ("highly accurate", "precise boundaries"); most tools headline precision, F1 or MCC |

**No tool of this style was found to combine:**

- a recall target on known elements, reported with an interval on a preregistered held-out set ([the held-out test](lockbox2.md));
- modes priced in candidates per genome, the follow-up budget, rather than in precision or F1;
- ranked candidates with triage fields for each call (checkpoint support, the source of each end, a range on every estimated end), meant to be followed up alongside annotation tools.

That combination, the contract, is the distinct part. Whether it helps discovery is untested: no candidate has been confirmed in the lab.

**Three design facts** set prophound apart from its nearest neighbors. Each is a fact about design, not a measured advantage.

1. **It learns integrated elements in their host context.** Jaeger masks prophages out of its training genomes with N, so it learns free phage against host. LAMBDA's models and Seeker train on free phage genomes against bacterial genomes. prophound trains on prophages as they sit in the chromosome, with host sequence within 5 kb of each element as hard negatives, which is the actual task: telling an integrated element from its own flanking host. The counterpart: their positive labels are independent of any prophage caller, while prophound's come from geNomad. Gene-based callers also learn from genes in host genomes (geNomad's provirus model; PhageBoost), so this holds among nucleotide-level models.
2. **Boundaries use no gene calling anywhere in the pipeline.** Jaeger always runs gene prediction (pyrodigal-gv) to move its ends to intergenic positions, and uses tRNAs and integrases when an annotation is given. geNomad, PhageBoost and PIDE all start from gene calls. LAMBDA's pipeline is gene-free, but it is a research pipeline rather than a packaged tool, and its edges sit on the window grid.
3. **It scans fragmented input by default.** Jaeger's prophage mode scans only contigs over 500 kb by default (`--lc 500000`). prophound scans contigs down to 3 kb, with thresholds calibrated for each window size. This is a difference in defaults, not in performance: with `--lc` lowered, Jaeger keeps as much recall on MAG-like assemblies, and so does geNomad ([Comparisons](tool-comparisons.md#fragmented-input-g5-g5b-)).

## Where the boundary stages come from

- **Att sites from direct repeats** are long established. Phage_Finder (Fouts 2006) found them with FASTA33, MUMmer or BLASTN and reports attL and attR. PhiSpy writes att coordinates and sequences. DEPhT finds repeats on both flanks with BLASTN and scores the pairs (cores of 8–39 bp, mismatches allowed). Prophage Tracer gets exact att ends from split sequencing reads. Each pairs the repeat search with gene calls, homology or reads. prophound accepts exact repeats only, so it misses mismatched cores that DEPhT would accept. That is the cost of exactness, and part of why only about one call in seven on known elements gets an att end.
- **Composition change-points at element edges** come from genomic-island tools. Alien_Hunter (Vernikos and Parkhill 2006) refines island boundaries with a change-point step over composition scores. The cumulative GC profile places domain boundaries to a single base (Zhang and Zhang), and Zisland Explorer and MTGIpick build on it. All are gene-free.
- **Ranges on estimated ends:** none of the surveyed tools reports a boundary uncertainty.

What is distinct is the combination: a packaged prophage caller that sets its ends from a model's score track, att repeats and composition, with no gene calls.

## The landscape

| tool | model input | gene calling? | database at run time? | how regions and boundaries are set | relation to prophound |
|---|---|---|---|---|---|
| **geNomad** `find-proviruses` | genes | yes | yes (marker profiles; MMseqs2 for integrases) | a conditional random field over marker-gene specificity; edges extended to tRNAs (ARAGORN) and integrases | the source of prophound's training labels; annotation-based; tuned for precise calls |
| **Jaeger** (v1.26.5) | 2 kb nucleotide windows, six-frame translated, dilated CNN | yes, for boundaries (pyrodigal-gv) | no | see [Jaeger in detail](#jaeger-in-detail) | the nearest method; localizes prophages, on contigs over 500 kb by default, with gene-snapped edges |
| **LAMBDA** (benchmark) | 2/4/8 kb windows, 9 genomic language models (1 M–7 B parameters) | no | no | see [LAMBDA in detail](#lambda-in-detail) | the closest pipeline, but a research benchmark; its curated truth set is used in prophound's evaluation |
| **PhageBoost** | genes | yes (Pyrodigal) | no (gene features + XGBoost) | per-gene probability over 0.9, at least 10 genes, gaps allowed, smoothing | the closest framing (database-free, discovery); needs gene calls; no recall target |
| **PIDE** | proteins | yes | no homology search; fine-tuned ESM-2 (650 M) | predicted phage ORFs clustered by gene density | database-free, protein-level; handles MAGs |
| **PhageTAILor** (preprint) | genes | yes | yes (6 homology detectors) + LightGBM | homology candidates, then classified | adds classes of tail-like elements; homology-based |
| **PHI** (preprint) | | | | a Galaxy workflow wrapping existing detectors | an integration layer, not a detector |
| **Seeker** | sequence | no | no (LSTM) | whole sequence, phage or bacterium | its authors state it was not trained or tested for prophages inside bacterial genomes |
| **DeepVirFinder**, **PPR-Meta** | contigs | no | no | each contig viral or not | contig classifiers; they do not localize |
| **DEPhT** | genes | yes | partly | shell-genome trimming | discovery-oriented; trained per genus (*Mycobacterium*); degrades on assemblies of over 100 contigs |
| **PhiSpy**, **DRAD**, **DBSCAN-SWA** | mixed | mostly yes | mixed | composition or sliding windows | predecessors; composition-based detection is decades old |
| **Prophage Hunter** | genes | yes | optional | regions plus an activity score | its similarity-free option exists "to find novel phages" |

### Jaeger in detail

Read from its source (v1.26.5):

- **Scoring.** Windows of 2 kb (`--fsize`), stride 2 kb by default, with four classes (phage, bacteria, eukarya, archaea). Prophages in its training genomes were masked with N.
- **Segmentation,** only on contigs longer than `--lc` (default 500 kb):
  1. Seeds are windows with a phage logit above `--sensitivity` (default 1.5).
  2. Each seed extends while the region's mean phage logit stays above its mean host logit.
  3. A kernel change-point detector proposes segments whose mean exceeds the sensitivity.
  4. The two region sets are merged and filtered to 6–250 kb with a coefficient of variation under 0.6.
- **Boundaries** (always run):
  1. pyrodigal-gv predicts genes, and each end moves outward to the nearest intergenic position (at most 4 kb).
  2. With a GenBank annotation, ends snap to the nearest integrase or tRNA within 20 kb where the phage score still exceeds the host score.
  3. Optionally, regions sparser in genes than the host are dropped (complete genomes only).
- **Accuracy:** no published prophage-level figures were found.

So localizing prophages from a nucleotide model is not distinct in itself; the distinctions that hold against Jaeger are the three design facts and the contract.

### LAMBDA in detail

- Windows of 2, 4 and 8 kb with 50% overlap, collapsed by majority vote, then a per-genome z-score, smoothing, a threshold, gap merging and a length filter. No gene calling; boundaries at window resolution.
- Training data: free phage genomes (INPHARED) against GTDB bacteria, split by cluster and genus.
- Truth: **380 expert-curated regions in 80 genomes** (the Casjens gold standard, PHASTER's manual annotations and Gauthier et al.). geNomad scores MCC 0.794 and recall 87.3%; the best genomic language model, Evo2, MCC 0.722 and recall 71.8%. The paper concludes that homology-based tools still lead.
- Its files carry every published tool's predictions on the curated genomes, which is what makes the LAMBDA comparison in [Comparisons](tool-comparisons.md#curated-recall-lambda) possible.
- Code: [github.com/leannmlindsey/LAMBDA](https://github.com/leannmlindsey/LAMBDA). Data: Zenodo 10.5281/zenodo.19236552.

### Philympics

A community gold standard: 78 manually curated genomes with gene-level phage tags, and a Snakemake framework that accepts new tools (MIT license). Published gene-level results: PhiSpy has the best F1 (0.733; recall 0.731), VirSorter2 the highest recall (0.766, precision 0.399), and Seeker, the only nucleotide deep-learning tool included, F1 0.093. geNomad is not in it. prophound's fragmented-input test (G5b) used Philympics genomes ([Comparisons](tool-comparisons.md#fragmented-input-g5-g5b-)).

## Comparing figures across tools

- **Recall figures are not comparable across truth sets or matching rules.** prophound's headline figures count agreement with geNomad on held-out genomes. Only a run on another benchmark's own truth set, scored by one rule for every tool, supports a side-by-side comparison; that is what [Comparisons](tool-comparisons.md) reports.
- **Agreement with geNomad is not accuracy.** geNomad itself covers about 82% of LAMBDA's curated prophages ([Claims and evidence](claims-and-evidence.md)).

## Sources

- geNomad: [pipeline docs](https://portal.nersc.gov/genomad/pipeline.html), [provirus identification](https://portal.nersc.gov/genomad/provirus_identification.html), [post-classification filtering](https://portal.nersc.gov/genomad/post_classification_filtering.html), [Nat Biotechnol 2023](https://www.nature.com/articles/s41587-023-01953-y), [issue #110](https://github.com/apcamargo/genomad/issues/110)
- Jaeger: [bioRxiv 2024](https://www.biorxiv.org/content/10.1101/2024.09.24.612722v1), [GitHub](https://github.com/Yasas1994/Jaeger), [releases](https://github.com/Yasas1994/Jaeger/releases), [docs](https://jaeger.readthedocs.io/en/latest/)
- LAMBDA: [NAR Genom Bioinform 2026](https://academic.oup.com/nargab/article/8/3/lqag103/8787816), [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC13041943/), [GitHub](https://github.com/leannmlindsey/LAMBDA), [interactive benchmark](https://leannmlindsey.github.io/lambda-benchmark/)
- PhageBoost: [NAR Genom Bioinform 2021](https://academic.oup.com/nargab/article/3/1/lqaa109/6066536), [GitHub](https://github.com/ku-cbd/PhageBoost)
- PIDE: [Genome Biology 2025](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-025-03733-0), [preprint](https://www.biorxiv.org/content/10.1101/2025.02.26.640259.full.pdf)
- PhageTAILor: [bioRxiv 2026](https://www.biorxiv.org/content/10.64898/2026.08.24.746745v1.full); PHI: [bioRxiv 2025](https://www.biorxiv.org/content/10.64898/2025.12.02.691814v1)
- Seeker: [NAR 2020](https://academic.oup.com/nar/article/48/21/e121/5921300)
- DEPhT: [NAR 2022](https://academic.oup.com/nar/article/50/13/e75/6572362), [att scheme](https://pmc.ncbi.nlm.nih.gov/articles/PMC9303363/figure/F5)
- PhiSpy: [NAR 2012](https://academic.oup.com/nar/article/40/16/e126/1027055), [PyPI](https://pypi.org/project/PhiSpy/4.1.7); DRAD: [PLOS One](https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0001193); DBSCAN-SWA: [published](https://pmc.ncbi.nlm.nih.gov/articles/PMC9061938/), [preprint](https://www.biorxiv.org/content/10.1101/2020.07.12.199018.full.pdf)
- Prophage Hunter: [NAR 2019](https://academic.oup.com/nar/article/47/W1/W74/5494712); Prophage Tracer: [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC8682789/); Phage_Finder: [documentation](https://phage-finder.sourceforge.net/documentation.htm)
- Alien_Hunter: [manual](https://manpages.debian.org/testing/alien-hunter/alien_hunter.1.en.html); Zisland Explorer: [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC6454522); MTGIpick: [publication](https://utsouthwestern.elsevierpure.com/en/publications/mtgipick-allows-robust-identification-of-genomic-islands-from-a-s)
- Philympics: [paper](https://f1000research.com/articles/10-758), [framework](https://github.com/linsalrob/ProphagePredictionComparisons)
