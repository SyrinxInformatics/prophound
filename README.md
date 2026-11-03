# Prophound

`prophound` is a prophage screen for bacterial genomes, built for exploration. It marks regions that look like integrated prophages directly from the DNA sequence, with no gene calling and no reference database, and returns them ranked as candidates for follow-up. Its purpose is to point researchers at regions worth checking in the lab, including ones databases don't contain. It is meant to run in supplement to annotation tools.

## What it is not

- **Not an annotation tool.** It reports where a candidate is and how sure its ends are. It names no genes and does not classify phages.
- **Not a source of confirmed new prophages.** No candidate has been confirmed in the lab.
- **Not more accurate than annotation tools.** It finds more known prophages by making more calls. Per call, it is not shown to be ahead.
- **Not precise, by design.** Most candidates outside known elements show no detectable signal. That is the price of finding nearly every known element.
- **Not for eukaryotic or archaeal input.** It was trained and tested on bacteria.

## How it works

1. **Score.** A neural network scores 10 kb windows along each contig, 1 kb apart, for how prophage-like they look. Three separately trained copies score every window.
2. **Call.** Runs of high-scoring windows become candidates. In the default mode, all three copies must agree.
3. **Place the ends.** Where a candidate is flanked by an att site (the repeat a phage leaves when it integrates), its ends sit exactly on it. Otherwise each end goes where the DNA's composition shifts from host-like to element-like, with a range.
4. **Rank.** Candidates are listed best first.

No genes are predicted and nothing is looked up.

## What it is based on

- **Training data:** 2,674 complete RefSeq bacterial genomes from over 2,200 genera, with prophages labeled by geNomad 1.9.4 (database 1.9). 1,871 of them trained the model. The model never sees geNomad's gene calls or marker HMMs, but it learned geNomad's idea of a prophage, blind spots included.
- **Model:** a convolutional network that reads both DNA strands identically. Three training replicates run as an ensemble and ship with the package (16.5 MB).
- **Ends:** att direct-repeat searches (as in Phage_Finder, PhiSpy and DEPhT) and composition change-points (as in Alien_Hunter), applied to the model's score track.
- **Testing:** settings were chosen on development samples, then tested once on 150 genomes from unseen genera, with the test fixed in advance ([`docs/lockbox2.md`](docs/lockbox2.md)). It was also compared with published tools on expert-curated and literature-named prophages.

## Install

### Install pixi

prophound is developed with [pixi](https://pixi.sh), which installs Python, PyTorch and every other dependency into a project environment. Follow the [pixi installation instructions](https://pixi.sh/latest/#installation), then check it:

```bash
pixi --version
```

### Install prophound

```bash
git clone https://github.com/SyrinxInformatics/mge-detect.git
cd mge-detect
pixi install
```

The model (three checkpoints, 16.5 MB) ships with the package, so there is nothing else to download.

**NOTE:** under pixi, run every command below as `pixi run prophound ...` from the `mge-detect` directory, or open a shell in the environment first with `pixi shell`.

**With pip instead:** into an environment with Python 3.11 or later (PyTorch is installed from PyPI):

```bash
pip install git+https://github.com/SyrinxInformatics/mge-detect.git
```

### Check the install

```bash
prophound --version
pixi run smoke        # runs the example genome end to end and checks its calls
```

## Usage

### Quickstart

```bash
prophound predict --input examples/hin_rd_kw20_region.fna --outdir example_out
```

The example is 220 kb of *Haemophilus influenzae* Rd KW20 holding its two geNomad-called prophage regions. It takes a few seconds on a GPU and a few minutes on CPU. The terminal ends with the candidates, best first:

```text
2 candidate prophage regions (default mode), best first:
  rank  contig              start-end (1-based)  length   edges      support
  1     hin_rd_kw20_region  121,179-152,504      31.3 kb  estimated  3/3 models
  2     hin_rd_kw20_region  56,749-76,926        20.2 kb  estimated  3/3 models
```

[`examples/README.md`](examples/README.md) walks through this output.

### Step by step

#### 1. Prepare the input

- **One genome per file.** `--input` takes FASTA, GenBank/RefSeq, EMBL, or GFF3 with an embedded `##FASTA` section (as Prokka and Bakta write), optionally gzip, bzip2 or xz compressed. The format is read from the content, so NCBI downloads work as they are.
- **Complete genomes, drafts and MAGs all work.** Contigs of 3–10 kb are scanned with their own thresholds. Contigs under 3 kb are never scanned, and the run says how many there were.
- **Bacterial input only.** prophound was trained and tested on bacteria. Remove host or eukaryotic contigs from a metagenome or assembly first, and treat calls on archaea with extra caution. A warning appears when a whole input looks eukaryotic ([Warnings and errors](#warnings-and-errors)).
- **Assembled sequence only.** Reads (FASTQ), alignments (BAM, SAM), protein sequence and GFF3 without sequence are refused with a message saying what to give instead.

#### 2. Choose a mode

`--mode` sets how many candidates to make, that is, the follow-up budget:

| mode | rule | candidates per genome | known elements found | when to use it |
|---|---|---|---|---|
| `exploratory` | window score 0.90; every call kept | about 37 | 0.98 | discovery runs, where a missed element costs more than a wasted check |
| **`default`** | 0.95; kept only if all three checkpoints make it | about 12 | 0.99 | most uses |
| `conservative` | 0.995; all three checkpoints | about 5 | 0.90 | a small follow-up budget |

The figures are from a held-out test of 150 complete genomes from genera unseen in training ([About prophound](docs/overview.md#modes-measured)). In Exploratory mode, triage candidates by `support` first.

#### 3. Run it

One genome. *E. coli* O157:H7 Sakai carries about 18 known prophages:

```bash
curl -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/865/GCF_000008865.2_ASM886v2/GCF_000008865.2_ASM886v2_genomic.fna.gz
prophound predict --input GCF_000008865.2_ASM886v2_genomic.fna.gz --outdir sakai_out
```

This takes about a minute on an Apple M4 GPU and gives 27 Default candidates, several of them with att-site ends.

Many genomes. Each genome gets its own subdirectory, plus combined calls and a `batch_summary.tsv`:

```bash
prophound predict --input genomes/ --outdir out                 # every genome in a directory
prophound predict --input a.fna --input b.fna --outdir out      # or repeat --input
prophound predict --input-list genomes.txt --outdir out         # or list paths in a file, one per line
```

- **Device.** `--device auto` (the default) uses a GPU when there is one: CUDA, then Apple MPS, then Intel XPU, then CPU. A 5 Mb genome takes under a minute on a laptop GPU and about an hour on CPU ([Run time](#run-time)). If a GPU runs out of memory, lower `--batch-size`.
- **On a cluster,** give each job a batch (`--input-list`) and `--device cuda`. The model ships with the package, so jobs need no network access.
- **Re-running** into an outdir that already holds outputs is refused; add `--force` to replace them.

#### 4. Check the terminal summary

Before the candidate list, the run prints `[summary]` lines: what was read (format, contigs, ambiguity codes), how many calls were made and how many the seed vote set aside, and the device and run time. `[warn]` lines flag anything that needs attention; the next section lists them. `-q` keeps only the warnings, errors and summary.

#### 5. Read the candidates

Open `regions.tsv`. [Interpreting the results](#interpreting-the-results) below gives the short version, and [How to read your results](docs/interpreting-results.md) the full guide, including how to take a candidate to the lab.

### Warnings and errors

| message begins | what it means | what to do |
|---|---|---|
| `this input may not be a prokaryotic genome` | eukaryotic telomere repeats at two or more contig ends, or over 17 Mb in total | remove eukaryotic contigs; for a metagenome or several genomes in one file, ignore it, or add `--no-domain-check` |
| `N contigs shorter than 10,000 bp` | short contigs were scanned with shorter windows; those under 3 kb were not scanned | nothing, unless much of the assembly is under 3 kb |
| `N calls set aside in regions.filtered.tsv: at least 50% N` | calls lying mostly on gap sequence were dropped | check `regions.filtered.tsv` if a gap may hide an element |
| `N U bases read as T (RNA sequence?)` | the input holds U | check that it is a genome |
| `--outdir ... already holds prophound outputs` | a previous run's outputs are there | choose another `--outdir`, or add `--force` |
| `... looks like FASTQ sequencing reads`, `... is a BAM alignment file` | the input is reads or alignments | assemble first and give the contigs |
| `... looks like protein, not nucleotide, sequence` | the input is protein | give the genome sequence |
| `... has duplicate contig IDs` | two records share an ID | rename them |
| `N of M genomes failed` (batch) | some genomes could not be read; the others ran | the `[warn]` line for each failed genome says why |

### Common options

| option | what it does |
|---|---|
| `--outdir DIR` | where outputs go (required). A directory that already holds outputs is refused. |
| `--force` | replace the outputs already in `--outdir` |
| `--mode MODE` | `exploratory`, `default` (the default) or `conservative` |
| `--device DEV` | `auto` (the default), or `cuda`, `rocm`, `mps`, `xpu`, `cpu` |
| `--write-region-fasta` | also write the called sequences to `regions.fna` |
| `--no-domain-check` | skip the non-prokaryote check and its warning |
| `-q` / `-v` | quieter (warnings, errors and the summary) / more verbose (stage timings, full tracebacks) |

`prophound predict --help` lists every option. The advanced ones change thresholds, the ensemble, the call builder and the boundary stages; [Command line and packaging](docs/cli-and-packaging.md) describes them.

### Run time

A 5 Mb genome takes under a minute (about 48 s) on an Apple M4 GPU and about an hour on its CPU. A median 3.6 Mb genome takes about 27 s on an NVIDIA V100S. Details: [Benchmarks](docs/benchmarks.md).

## Output

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

A batch writes one subdirectory per genome, plus combined calls and `batch_summary.tsv`. Every column is defined in [`docs/output-files.md`](docs/output-files.md).

## Interpreting the results

Start with `regions.tsv`. [How to read your results](docs/interpreting-results.md) covers every column and how to take a candidate to the lab. In short:

- **Rank.** Candidates most like known prophages come first, so a novel element can rank low. In a discovery run, read past the top of the list.
- **Ends.** `att site` ends sit on an integration repeat and are exact to the base; about 1 call in 7 on known elements gets them, and a repeat can be chance. `estimated` ends are a median of about 2 kb off and come with a range.
- **Support.** How many of the three checkpoints make the call. In Exploratory mode, calls only one checkpoint makes are almost never known elements.
- **Length.** Long candidates carry a mobile-element or phage signal far more often than short ones, and may hold more than one element.
- **Scores** are not probabilities that a call is a prophage.

## Disclaimers

- **No candidate has been confirmed in the lab.** Calls are candidates for a researcher to check.
- **Agreement is not accuracy.** "Known prophages" in the headline figures are geNomad's calls, and geNomad also labeled the training data. It covers about 82% of expert-curated prophages, and prophound inherits its blind spots.
- **A score alone is not evidence.** Shuffled genomes with no genes left still got calls, and above 55% GC most of the shuffled sequence was called. Real bacterial genomes show no GC effect, but on unusual sequence a high score proves little.
- **Input should be bacterial.** Over half of a 10 Mb stretch of human chromosome was called, and over a third of two archaeal genomes.
- **Plasmids are often called whole** (a third of plasmids). Such a call marks the replicon and locates nothing inside it.
- **Small elements and fragments.** Default finds about a third of implanted 3–10 kb inoviruses. Contigs under 3 kb are never scanned. On MAG-like assemblies Default keeps 0.92–0.94 of its complete-genome recall.
- **Ends.** Estimated ranges held the true end in about 78% of cases on held-out genomes. About 5% of held-out elements shared a call with another element.
- **Platforms.** Verified on Apple M4 (CPU and GPU), an NVIDIA V100S and Linux x86 CPU; others are untested.

Every limitation, with its numbers, is in [About prophound](docs/overview.md#what-it-does-not-claim) and the [release text](docs/release-text-0.1.0.md).

## Why it is built this way

- **No database.** Annotation tools recognize prophages by resemblance to known genes and proteins, so they can miss elements unlike anything cataloged. A model that reads sequence directly does not need a complete catalog.
- **Recall first.** A missed element cannot be recovered, while an extra candidate costs one check. So prophound aims to find nearly every known element and pays with extra candidates.
- **Cost in candidates per genome.** Precision against geNomad would count every call geNomad doesn't make as an error, and those calls are where anything new would be. So each mode is priced in candidates per genome, the follow-up budget, chosen with `--mode`.
- **Alongside annotation tools.** On LAMBDA, adding prophound's calls to geNomad's recovered 16 more curated prophages and lost none, at about 15 extra candidates per genome (12 of the 16 from two genomes).

**What the evidence shows.**

- On 150 held-out genomes, Default found 99% of geNomad's prophages (95% interval 96–100%) at about 12 candidates per genome.
- On 39 genomes of expert-curated prophages (LAMBDA), Default found 93% and geNomad, PHASTER and PIDE about 82%. It made 1.7–4.1 times as many calls, and allowing for that, it is not shown to be ahead.
- About 4 in 10 Default candidates outside known elements carry a mobile-element or phage signal, against 1 in 10 for random stretches of the same genomes.

[Comparisons with established tools](docs/tool-comparisons.md) has every head-to-head result.

**What is new.** Database-free detection, nucleotide-level models and discovery framing already exist (PhageBoost, Jaeger, LAMBDA's models, geNomad's `--relaxed`). No tool of this style was found to combine a recall target measured on a preregistered held-out set with modes priced in follow-up cost. Whether that helps discovery is untested in the lab ([Related tools](docs/related-tools.md)).

## Documentation

| doc | what it covers |
|---|---|
| [About prophound](docs/overview.md) | how it works, what has been measured, what it does not claim, run time |
| [How to read your results](docs/interpreting-results.md) | the researcher's guide to the outputs, and taking a candidate to the lab |
| [Output files](docs/output-files.md) | every file and column |
| [Command line and packaging](docs/cli-and-packaging.md) | every option, the bundled model, installation and distribution |
| [Release text](docs/release-text-0.1.0.md) | a plain-language summary, then every figure with its interval, sample and source |
| [Claims and evidence](docs/claims-and-evidence.md) | which claims the evidence supports, and how each was tested |
| [The held-out test](docs/lockbox2.md) | the preregistered test behind the headline figures, and its result |
| [Comparisons with established tools](docs/tool-comparisons.md) | every head-to-head result against geNomad, PHASTER, PHASTEST, PIDE, Jaeger and others, with its verdict |
| [Related tools](docs/related-tools.md) | neighboring tools and how prophound differs |
| [Benchmarks](docs/benchmarks.md) | run time, memory and cross-device agreement |
| [Research integrity](docs/research-integrity.md) | how the evaluation guards against tuning on its own test data |
| [Project principles](docs/project-principles.md) | the framing and working rules behind every decision |

## Development

- **Layout.** The package is `prophound/`; the prophage detector is in `prophound/models/prophage/`, with the shipped code in `src/`, the bundled model in `model/`, and the evaluation scripts and their recorded outputs in `scripts/` and `analysis/` (not part of the package). The bioconda recipe is in `bioconda/recipes/prophound/`.
- **Tests.** `pixi run test` runs the fast tests; `pixi run smoke` runs the example end to end.
- **Commit hooks.** `pixi run install-hooks` once per clone. The hooks refuse commits and pushes that carry secrets, credential files, files over 10 MB, or data and checkpoints outside the places that ship them (`scripts/check_commit.py`).
- **CI** (`.github/workflows/ci.yml`) runs errors-only pylint, the tests and the example, the commit check and a gitleaks scan of the full history, and a bioconda lint and build.

## Versioning

prophound follows [Semantic Versioning 2.0.0](https://semver.org/). While it is in initial development (0.y.z), the command line, output files and defaults may change in any minor release. Changes are listed in [`CHANGELOG.md`](CHANGELOG.md).

## Citation

See [`CITATION.cff`](CITATION.cff).

## License and acknowledgments

Apache License 2.0; see `LICENSE` and `NOTICE`.

This work received computational support from UTSA's high-performance computing cluster operated by Tech Solutions.
