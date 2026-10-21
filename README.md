# prophound

`prophound` is a recall-first prophage screen for prokaryotic genomes. It localizes candidate integrated elements directly from nucleotide sequence, with no gene calling or reference database at run time. It is tuned to recover already-known elements and returns everything else as candidates, at a follow-up cost the user chooses. It is built to run alongside annotation-based tools and to surface what they may not call.

The name is a prophage hound, and a homophone of "profound". The project was called `mge-detect` until 2026-09-29. Other mobile-element types may follow the prophage detector, which is why the code keeps a `models/prophage/` layout.

**Status: version 0.1.0, not yet released.** This is research software. Read [What it does not claim](#what-it-does-not-claim) before using its calls.

## How it works

A reverse-complement-symmetric CNN scores 10 kb windows along each contig, sliding 1 kb at a time. Three training replicates run as an ensemble, and they ship with the package (16.5 MB in all). Calls are built from the score track: stretches where enough windows score high, merged across dips inside large elements. The model was trained on genomes labelled by geNomad 1.9.4, but it never sees geNomad's inputs (gene calls, marker HMMs), and nothing is looked up at run time.

## Install

From source, with [pixi](https://pixi.sh):

```bash
git clone https://github.com/SyrinxInformatics/mge-detect.git
cd mge-detect
pixi install
pixi run prophound --version
```

With pip, into an environment with Python 3.11 or later (installs PyTorch from PyPI):

```bash
pip install git+https://github.com/SyrinxInformatics/mge-detect.git
```

A bioconda package is planned; the draft recipe is in `recipes/prophound/`. Under pixi, prefix the commands below with `pixi run`.

## Quickstart

```bash
prophound predict --input examples/hin_rd_kw20_region.fna --outdir example_out
```

The example is 220 kb of *Haemophilus influenzae* Rd KW20 holding its two geNomad-called prophage regions. It takes a few seconds on a GPU and a few minutes on CPU. The terminal ends with the candidates, best first:

```
2 candidate prophage regions (default mode), best first:
  rank  contig              start-end (1-based)  length   edges      support
  1     hin_rd_kw20_region  121,179-152,504      31.3 kb  estimated  3/3 models
  2     hin_rd_kw20_region  56,749-76,926        20.2 kb  estimated  3/3 models
```

`regions.tsv` holds the same table with a range for each edge. **[How to read your results](docs/interpreting-results.md)** explains every column, what the edges mean, and how to take a candidate to the lab.

**Boundaries.** Where an element has an att site (a direct repeat left at both ends by integration), `edges` is `att site` and the ends are exact to the base, with no database or gene calling. About 15% of known elements have one. Elsewhere `edges` is `estimated`: the model's estimate, moved in Default and Conservative mode to where the sequence's composition turns from host-like to element-like (median error about 2 kb on our test samples), with a range for each end that holds the true end about 90% of the time. Neither example element has an att site; the calls' ends are 1–4 kb from the elements' (58,837–73,073 and 119,962–154,275, 1-based), inside the ranges. [`examples/README.md`](examples/README.md) walks through the output, and `pixi run smoke` checks it.

Other inputs:

```bash
prophound predict --input genome.gbk.gz --outdir out            # GenBank, EMBL, GFF3+FASTA; compressed is fine
prophound predict --input genomes/ --outdir out                 # a batch: every genome in a directory
prophound predict --input a.fna --input b.fna --outdir out      # or repeat --input, or --input-list FILE
prophound predict --input mag.fna --outdir out --mode exploratory
```

`prophound predict --help` lists every option; the advanced ones change thresholds, the ensemble or the call builder.

## Modes

Each mode is priced in calls per genome, the follow-up budget. Figures are from the held-out test set: 150 complete RefSeq genomes from genera absent from training, 202 geNomad elements, scored once.

| mode | rule | geNomad elements covered | calls per genome |
|---|---|---|---|
| `exploratory` | window score 0.90; every call kept | 0.965 | 32.4 |
| **`default`** | 0.95; kept only if all three checkpoints make it | **0.960** (95% interval 0.924–0.980) | 11.0 |
| `conservative` | 0.995; all three checkpoints | 0.901 | 4.9 |

"Covered" means calls cover at least half of an element (covered50). On the 70 held-out genomes with no geNomad element, Default still makes 6.8 calls per genome. That is the floor cost of a run.

In Exploratory mode, start triage with `seed_support`: calls only one checkpoint makes are almost never known elements.

## Outputs

| file | contents |
|---|---|
| `regions.tsv` | **start here:** the candidates, best first, in 1-based coordinates: where each is, how exact its ends are (`att site` or `estimated`, with ranges), model support and notes |
| `regions.detailed.tsv` | the same candidates with every column (scores, the model's own edges, the att core), 0-based like BED; what scripts should read |
| `regions.bed`, `regions.gff3` | the same calls; BED score is `round(1000 × mean_prob)` |
| `regions.filtered.tsv` | calls a filter set aside (the seed vote, or mostly N), with the reason |
| `windows.tsv` | every window's score: the track the calls are built from |
| `fine_windows.tsv` | the fine pass's window scores, around each call |
| `run_metadata.json` | settings, model and input checksums, versions, timing |
| `regions.fna` | the called sequences (`--write-region-fasta`) |

A batch writes one subdirectory per genome, plus combined calls and `batch_summary.tsv`. Every column is described in [`docs/output-files.md`](docs/output-files.md). An outdir that already holds outputs is refused unless `--force` is given.

`rank` puts the candidates most like known prophages first, so a novel element can rank low. The scores behind it are not probabilities that a call is a prophage.

## What has been measured

Every figure, its interval and its source are in [`docs/release-text-0.1.0.md`](docs/release-text-0.1.0.md). Briefly:

- **Agreement with geNomad, held out:** Default covers 0.960 of 202 elements (0.924–0.980). The test set was preregistered and scored once ([`docs/lockbox.md`](docs/lockbox.md)). Its lower bound is under 0.95, so 95% recall is not established.
- **Expert-curated truth (LAMBDA, 39 genomes):** Default covers 0.933 of 149 curated prophages, against 0.82 for geNomad, PHASTER and PIDE. It makes 1.5–3.7× their calls per genome, though, and per unit of cost it is level with them: after adjusting for chance, its paired differences to all three span zero.
- **Alongside geNomad:** adding prophound Default to geNomad recovers 16 more curated LAMBDA prophages and loses none, at about 15 extra calls per genome. Twelve of the 16 come from two genomes.
- **LAMBDA's own candidates:** prophound covers 18 of the 19 phages LAMBDA's reviewers judged likely but left out of its curated set (0.84 above chance). geNomad covers 15. Most of them are called by annotation tools too.
- **Off-target calls** carry mobile-element features 3.1× as often as matched background, mostly conjugation and integrase/transposase genes rather than phage structural genes. Many look like other integrative elements, not missed prophages.

## What it does not claim

- **No validated novel elements.** No call has been confirmed in the lab. Calls are candidates for follow-up.
- **Not more accurate than** geNomad or any other tool. It finds more at a higher cost.
- **Every agreement figure is with geNomad,** which also supplied the training labels. The model inherits its blind spots.
- **Small elements:** recall under 10 kb is low (0.36 on implanted inoviruses at Default) and partly measured on synthetic data.
- **Fragmented input:** contigs of 3–10 kb are scanned, and contigs under 3 kb never are. On MAG-like assemblies Default keeps 0.92–0.94 of its complete-genome recall, which does not beat geNomad.
- **Several elements may fall under one call** (9% of held-out elements).
- **Platforms:** results are verified on Apple M4 (CPU and GPU) and an NVIDIA V100S. Other platforms are expected to agree but are untested.

## Run time

On a MacBook Pro (Apple M4, 16 GB), Default mode:

| input | time | memory |
|---|---|---|
| a 5 Mb genome, Apple GPU (picked automatically) | ~35 s | ~0.5 GB |
| the same genome as a 374-contig assembly | ~19 s | ~0.5 GB |
| a 5 Mb genome, CPU only | ~45 min | ~1.4 GB |

A GPU makes a large difference. On an NVIDIA V100S a median 3.6 Mb genome takes ~27 s. Details: [`docs/benchmarks.md`](docs/benchmarks.md).

## Citation

See [`CITATION.cff`](CITATION.cff). Changes are listed in [`CHANGELOG.md`](CHANGELOG.md).

## Development

Training and large-scale evaluation run on the University of Texas at San Antonio's HPC cluster (ARC). This repository holds the CLI, the packaging and the evaluation record.

- `src/prophound/` is the package. The prophage detector is in `src/prophound/models/prophage/`, including the bundled model in `model/` and the dev-only evaluation scripts in `scripts/`.
- `docs/` holds the living topic docs, and [`SESSIONS.md`](SESSIONS.md) indexes the work. [`CLAUDE.md`](CLAUDE.md) says where everything is documented.
- `pixi run test` runs the fast tests, and `pixi run smoke` runs the example end to end.

## License and acknowledgements

Apache License 2.0; see `LICENSE` and `NOTICE`.

This work received computational support from UTSA's high-performance computing cluster operated by Tech Solutions.
