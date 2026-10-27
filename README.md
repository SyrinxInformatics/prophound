# prophound

`prophound` is a recall-first prophage screen for prokaryotic genomes. It finds candidate integrated elements directly from nucleotide sequence, with no gene calling or reference database at run time, and returns them ranked, for follow-up alongside annotation-based tools.

**NOTE:** version 0.1.0 is not yet released, and this is research software. Its calls are candidates, not annotations. Read [About prophound](docs/overview.md), and especially [what it does not claim](docs/overview.md#what-it-does-not-claim), before using them.

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

A bioconda package is planned; the draft recipe is in `recipes/prophound/`.

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

```
2 candidate prophage regions (default mode), best first:
  rank  contig              start-end (1-based)  length   edges      support
  1     hin_rd_kw20_region  121,179-152,504      31.3 kb  estimated  3/3 models
  2     hin_rd_kw20_region  56,749-76,926        20.2 kb  estimated  3/3 models
```

[`examples/README.md`](examples/README.md) walks through this output.

### A complete genome

*E. coli* O157:H7 Sakai carries about 18 known prophages. prophound reads compressed input, so the NCBI file can be used as downloaded:

```bash
curl -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/865/GCF_000008865.2_ASM886v2/GCF_000008865.2_ASM886v2_genomic.fna.gz
prophound predict --input GCF_000008865.2_ASM886v2_genomic.fna.gz --outdir sakai_out
```

This takes about a minute on an Apple M4 GPU and gives 27 Default candidates, several of them with att-site ends.

### Input

`--input` takes FASTA, GenBank/RefSeq, EMBL, or GFF3 with an embedded `##FASTA` section (as Prokka and Bakta write), optionally gzip, bzip2 or xz compressed. The format is read from the content.

```bash
prophound predict --input genome.gbk.gz --outdir out            # any supported format, compressed or not
prophound predict --input genomes/ --outdir out                 # a batch: every genome in a directory
prophound predict --input a.fna --input b.fna --outdir out      # or repeat --input
prophound predict --input-list genomes.txt --outdir out         # or list paths in a file, one per line
```

Contigs of 3–10 kb are scanned with their own thresholds, so draft assemblies and MAGs work as input. Contigs under 3 kb are never scanned, and the run reports how many there were.

### Modes

`--mode` sets how many candidates to make, that is, the follow-up budget:

| mode | rule | when to use it |
|---|---|---|
| `exploratory` | window score 0.90; every call kept | discovery runs, where a missed element costs more than a wasted check; start triage with `support` |
| **`default`** | 0.95; kept only if all three checkpoints make it | the recommended mode |
| `conservative` | 0.995; all three checkpoints | a small follow-up budget |

On a held-out test of 150 complete genomes these made about 32, 11 and 5 calls per genome. Recall and cost per mode are in [About prophound](docs/overview.md#modes-measured).

### Common options

| option | what it does |
|---|---|
| `--outdir DIR` | where outputs go (required). A directory that already holds outputs is refused. |
| `--force` | replace the outputs already in `--outdir` |
| `--mode MODE` | `exploratory`, `default` (the default) or `conservative` |
| `--device DEV` | `auto` (the default: CUDA, then Apple MPS, then Intel XPU, then CPU), or `cuda`, `rocm`, `mps`, `xpu`, `cpu` |
| `--write-region-fasta` | also write the called sequences to `regions.fna` |
| `-q` / `-v` | quieter (warnings, errors and the summary) / more verbose (stage timings, full tracebacks) |

`prophound predict --help` lists every option. The advanced ones change thresholds, the ensemble, the call builder and the boundary stages. They are described in [`docs/cli-and-packaging.md`](docs/cli-and-packaging.md).

### Run time

A GPU makes a large difference. A 5 Mb genome takes about 35 s on an Apple M4 GPU and about 45 min on its CPU, and a median 3.6 Mb genome about 27 s on an NVIDIA V100S. See [`docs/benchmarks.md`](docs/benchmarks.md).

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

### Reading the results

- **Ends.** `edges` = `att site` means both ends sit on a direct repeat left by integration, exact to the base. `estimated` ends can be several kb off, and `regions.tsv` gives a range for each.
- **Rank.** `rank` puts the candidates most like known prophages first, so a novel element can rank low. The scores behind it are not probabilities that a call is a prophage.
- **Support.** `support` is how many of the three checkpoints make the call. In Exploratory mode, calls only one checkpoint makes are almost never known elements.

[How to read your results](docs/interpreting-results.md) explains every column, how far to trust each kind of end, and how to take a candidate to the lab.

## Documentation

| doc | what it covers |
|---|---|
| [About prophound](docs/overview.md) | how it works, what has been measured, what it does not claim, run time, development |
| [How to read your results](docs/interpreting-results.md) | the researcher's guide to the outputs |
| [Output files](docs/output-files.md) | every file and column |
| [Release text](docs/release-text-0.1.0.md) | every figure with its interval, sample and source |
| [Claims and evidence](docs/claims-and-evidence.md) | which claims the evidence supports, and how each was tested |
| [Related tools](docs/related-tools.md) | neighbouring tools and how prophound differs |
| [Benchmarks](docs/benchmarks.md) | run time, memory and cross-device agreement |
| [Project principles](docs/project-principles.md) | the framing and working rules behind every decision |
| [Research integrity](docs/research-integrity.md) | the methods audit and the held-out test design |

The other files in `docs/` are the living record of each topic, from boundaries to training data.

## Versioning

prophound follows [Semantic Versioning 2.0.0](https://semver.org/). While it is in initial development (0.y.z), the command line, output files and defaults may change in any minor release. Changes are listed in [`CHANGELOG.md`](CHANGELOG.md).

## Citation

See [`CITATION.cff`](CITATION.cff).

## License and acknowledgements

Apache License 2.0; see `LICENSE` and `NOTICE`.

This work received computational support from UTSA's high-performance computing cluster operated by Tech Solutions.
