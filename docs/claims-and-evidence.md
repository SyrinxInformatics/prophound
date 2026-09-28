# Claims and evidence

Living doc, updated in place. Built 2026-09-28. It says which novelty and improvement claims the v1 write-up can make, what evidence each needs, and what the tool itself must do so that its exploratory outputs can be trusted. The positioning argument and related tools are in `docs/related-tools.md`. Release tasks are in `docs/release-checklist.md`.

## The principle (user, 2026-09-28)

The tool should make a novelty or improvement claim confidently without overselling. It should **not rest on the claim alone.** It must be technically robust and easy to use, and its accuracy on *known* elements must be measured and published, so that its exploratory outputs can be considered with confidence. The order that follows from this: known accuracy first, then robustness and usability, then the improvement claims, and the concept last. The concept is only as credible as the measured layer under it.

## Release gate: new, or measurably better (user, 2026-09-28)

> "If this isn't doing anything new, it at least needs to be better."

The concept (rung 1 below) is new only as a *contract*, and the methods are not new. So v1 has to show it is **better than existing tools on at least one measured axis**, on independent truth or identical input, not only on agreement with geNomad. **Scope rule (user, 2026-09-28): an axis is in the gate only if measuring it will not cause major delay.** The axes are sorted by what they cost to measure.

### In the gate: all measured by one run, no new tools to install

One ensemble run of mge-detect (all three modes) on LAMBDA's 39 clean genomes, scored locally against the published predictions of the other tools. Prerequisites: a 21-mer containment check against training positives (reuse `small_element_spikein.py`'s `disjoint`) and confirming LAMBDA's coordinate frame. The cost is about one GPU job the size of the 42-genome runs (1–2 h). Scoring against the other tools' predictions is plain Python and **can run on ARC** in the same job or a CPU follow-up (user, 2026-09-28), so it adds no local turnaround. The LAMBDA JSONs are about 20 MB and are pushed alongside.

| # | axis | bar to clear (LAMBDA clean subset) | outlook |
|---|---|---|---|
| G1 | **Complement to geNomad** | geNomad ∪ mge-detect gains **more than +8** curated elements over geNomad alone (0.809, 123/152). +8 is the best existing pairing (PHASTER), at +3.9 calls/genome. | **Plausible**, most likely in Exploratory. It is the direct test of the concept. Report the extra calls/genome with it. |
| G2 | **Sequence-only, at matched cost** | top-*k* calls per genome (k matched to ~3.5) reach covered50 **above 0.727** (best genomic LM, NTv2; Evo2 0.724) | **Plausible.** It is trained on integrated elements in host context, whereas the LMs trained on free phage genomes. |
| G3 | **Overall curated recall** | covered50 **above 0.818** (PHASTER; PIDE 0.812, geNomad 0.809), with calls/genome | **Unlikely at Default**, because the model inherits geNomad's misses. Possible in Exploratory, at about 7× the calls/genome, and must be stated that way. |
| G4 | **Footprint** | install size versus geNomad's database: `genomad_db_v1.9.tar.gz` is **842 MB compressed** (Zenodo 14886553; measure the unpacked size). Our side is 3 checkpoints plus the package. | **Likely.** It counts only alongside a G1–G3 recall figure, since a small tool that finds less is not better. |

### Time-boxed: in the gate only if setup goes smoothly

| # | axis | bar | cost and fallback |
|---|---|---|---|
| G5 | **Fragmented input** | on identical MAG-like shreds of the 39 clean genomes, mge-detect retains more complete-genome recall than geNomad and Jaeger `-p` (`--lc` lowered, and stated) | **The most likely win**, but the only axis that needs installing other tools. ARC lacks MMseqs2 and Prodigal (preflight), so run the comparators locally, where both are conda/pip-installable. **Time-box: about one day for installs.** If it slips, fall back to a self-referential figure ("retains X% of its complete-genome recall"), which is measured already on eval genomes (covered 0.939) and needs no comparator. That is not a "better" claim. |
| G6 | **CPU runtime** | wall-clock per genome versus geNomad on the same machine | Comes free with G5's local geNomad install. Without it, use no runtime comparison; published figures come from other hardware. |

### Measured by the same run, where the tool is probably not better (report plainly)

| axis | existing level (LAMBDA clean) | our reference point | reading |
|---|---|---|---|
| **Boundary precision** | median IoU of detected elements: geNomad **0.921**, PIDE 0.912, VIBRANT 0.885, PHASTER 0.831. IoU ≥ 0.5 on 0.72 / 0.76 / 0.68 / 0.70 of elements | against geNomad truth, our median IoU is about 0.69–0.74 (`docs/scientific-validity.md`) | **Probably worse.** Annotation tools snap edges to genes, tRNAs and integrases; ours come from 1 kb bins of the score track. State it as a known limitation. It is a real cost of skipping gene calling. |
| **Small elements (<10 kb)** | 18 of 154 clean LAMBDA elements | 5–8 kb covered50 0.857 on geNomad truth (n=7) | Known weakness, already in the limitations. |
| **Calls/genome at the same recall** | 3.4–4.5 for most tools | Default 14.5, Conservative 9.4 | Costlier by design. Frame it as follow-up budget, never hide it. |

### Not in the gate (would cause delay)

- **Ambiguity robustness versus other tools.** It needs spiked input run through every comparator. Our own dose-response (`821396`) is reported as a robustness fact, not a comparison.
- **Contributing mge-detect to the Philympics framework** (a Snakemake rule and a PR to an external repo, on their timeline). Post-v1. *Scoring* Philympics' 23 clean genomes ourselves can be added to the G1–G4 job if it costs less than a day, since the published CDS-level results are in its repo (`jupyter_notebooks/all_tptn.tsv`).
- **Anything requiring retraining** (small elements, boundary head): deferred past v1 already.

### Decision rule

- **If any of G1–G3 or G5 clears its bar**, the write-up leads with the axes that cleared: each with its bar, our number, the calls/genome and the n. G4 and G6 go alongside, never alone.
- **Every axis measured is reported**, including those where the tool is not better (boundary precision, small elements, cost per genome). Reporting the losses is what makes the wins credible.
- **If none clears**, v1 cannot claim to be new or better, and the release decision goes back to the user.

## Two layers of "known accuracy"

This is the finding that most changes how results should be written up.

| layer | what it measures | current number |
|---|---|---|
| **Agreement with geNomad** | how well the tool reproduces its label source on held-out genomes | covered50 0.961 (ensemble, Default, n=179) |
| **Agreement with curated truth** | how well it finds prophages that experts annotated | **not yet measured**. geNomad itself reaches only **0.81** covered50 on LAMBDA's clean subset (below) |

geNomad covers about 81% of expert-curated prophages. A tool that matches geNomad at 96% will therefore land somewhere near 80% on curated truth. The 0.961 figure measures agreement, not accuracy. **Publishing it alone would oversell the tool, and the release must report both layers.** Curated truth is the number that tells a user how much of what is really there the tool finds.

## External curated benchmarks available now (checked 2026-09-28)

Two community gold standards exist. Their genome lists were matched against our corpus by assembly accession. Overlap tables: `src/mge_detect/models/prophage/analysis/external_benchmarks/`.

| benchmark | genomes | curated elements | not in our corpus | in our test split | in train/val (exclude) |
|---|---|---|---|---|---|
| **LAMBDA** (2026; Casjens gold standard, PHASTER manual annotations, Gauthier et al.) | 76 | 344 regions | 34 genomes (118 regions) | 5 | 37 |
| **Philympics** (2021; 78 manually curated genomes, CDS-level `is_phage` tags) | 73 files with coordinates (69 mapped) | 275 rows (includes 1-gene fragments) | 17 genomes (77 rows) | 6 | 46 |

**LAMBDA is the stronger instrument.** Its per-genome files include the curated truth, **every published tool's predicted regions** (geNomad, PHASTER, VIBRANT, VirSorter2, Phigaro, PhiSpy, PhageBoost, PIDE, and eight genomic language models) and 292 off-truth candidate regions that its authors classified from PHROG gene evidence. A head-to-head comparison therefore needs only an mge-detect run. The other tools are already scored on the same genomes.

### LAMBDA baseline on the clean subset (39 genomes, 154 curated regions, 18 under 10 kb)

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

These are the bars. **mge-detect's complement claim holds only if it gains more than +8 elements over geNomad** (the best existing pairing), and the gain is reported with its extra calls/genome.

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
| Agreement with geNomad on held-out genomes, stratified by size | `docs/model-quality-baseline.md` | "agrees with geNomad on 96% of held-out elements (covered50, CI)", **never** "96% accurate" |
| Robustness facts measured so far | RC exact; grid phase 1.7%; short-contig scan recovers covered 0.939 on shredded input; ambiguity `impute` | plain statements, with the numbers |

### Rung 2: improvement claims that need one new run

| claim | run needed | wording by outcome |
|---|---|---|
| **B1. Recall on independent curated truth** | mge-detect (ensemble, all three modes) on LAMBDA's clean subset, and on Philympics' clean subset | Above 0.818 covered50: "highest recall of the N tools compared, at X calls/genome versus Y". Within the CI of the leaders: "comparable to annotation-based tools without gene calling or a database". Below: report it plainly, since it is still known accuracy. |
| **B1'. Against sequence-only models at matched cost** | same run. Rank calls by score and keep the top *k* per genome, to match the genomic LMs' ~3.5 calls/genome | if higher: "exceeds published sequence-only (genomic language model) approaches at matched follow-up cost". **Matched cost is required.** Comparing Exploratory at 29 calls/genome with models at 3.5 would be an unfair comparison. |
| **B2. Complement to geNomad** | same run, geNomad ∪ mge-detect | above +8 elements: "adds more curated prophages to geNomad's calls than any other single tool tested, at +Z calls/genome". Otherwise report the gain against the table above. **This is the direct test of the concept.** |
| **B3. Fragmented input** | shred the clean LAMBDA genomes to a MAG-like contig-length distribution. Run mge-detect, Jaeger `-p` (with `--lc` lowered, and stated) and geNomad on the same shreds | "retains X% of its complete-genome recall on fragmented assemblies". A comparative claim is allowed only with the others run on identical input. It is a good candidate: DEPhT degrades on assemblies with more than 100 contigs, and Jaeger scans only contigs over 500 kb by default. |
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
| **`--version`** | **missing**: `mge-detect --version` is rejected. It is needed for reproducibility and should also be written to `run_metadata.json`. |
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
- **Optional, calibrated**: an empirical "P(overlaps a known element)" per score/support bin, from held-out data. This must be labelled as agreement with known elements, not as the probability the call is real.

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
