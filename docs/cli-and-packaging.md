# Command line and packaging

How the `prophound` command works, what each option does, and how the package is built and distributed. For a first run, the [README](../README.md#usage) is the place to start; every output file and column is in [Output files](output-files.md).

## The command

```bash
prophound predict --input <genome> --outdir <dir> [--mode exploratory|default|conservative] [options]
prophound --version
```

`predict` is the only subcommand. The model ships with the package, so `--input` and `--outdir` are all a run needs.

- **Output streams.** Nothing goes to stdout. Progress bars, `[summary]`, `[warn]` and `[done]` lines go to stderr.
- **Exit codes.** 0 on success; 1 for an error the user can fix, reported as one `error:` line; 130 on Ctrl-C. A batch exits 1 if any genome failed.
- **`-q`** prints only warnings, errors and the summary. **`-v`** adds each stage's timing and the full traceback of an unexpected error.

## What a run does

1. **Read and check the input.** The format is detected from the content (FASTA, GenBank/RefSeq, EMBL, or GFF3 with `##FASTA`; gzip, bzip2 or xz). Reads, alignments, protein and sequence-free files are refused with a reason. Duplicate contig IDs, empty records and headers without an ID are errors. Ambiguity codes are counted per code; U is read as T, with a warning.
2. **Check that the input looks prokaryotic.** A warning is printed if two or more contig ends carry a eukaryotic telomere repeat, or if the input totals over 17 Mb. The run continues either way ([below](#the-non-prokaryote-warning)).
3. **Cut windows.** 10 kb windows every 1 kb on each contig. Contigs of 3–10 kb get 5 kb or 3 kb windows instead, whichever fits, at thresholds calibrated to the same false-call rate. Contigs under 3 kb are not scanned, and a warning reports them.
4. **Score.** Each of the three checkpoints scores every window on both strands. The ensemble score is the highest of the three (`--ensemble max`).
5. **Build calls.** A position is inside a call when at least 10% of the windows covering it reach the extend threshold (`--coverage-frac 0.1`). Calls within 50 kb of each other are merged when every window between them stays at or above 0.1 (`--max-gap-bp`, `--bridge-prob`), which repairs large elements whose score dips partway through.
6. **Vote.** A call is kept only if enough checkpoints make an overlapping call on their own (the mode's seed vote). Calls the vote drops go to `regions.filtered.tsv` with the reason.
7. **Refine.** A fine pass re-scores each call and its surroundings with 5 kb windows every 250 bp, splits calls across a clear background gap, and places the edges between the fine and coarse estimates.
8. **Finish the calls.** An edge within 5 kb of a contig end moves onto it and is flagged in `contig_edge`. Calls that are at least half N are set aside.
9. **Place the ends.** Each call is searched for an att site, a direct repeat with one copy near each edge. Where one is accepted, the call's ends are its two copies, exact to the base. Every other end moves to where the composition shifts from host-like to element-like, and gets a range.
10. **Rank and write.** Calls are ranked by their highest window score, then their mean, and every output is written atomically.

## Modes

Each mode is a seed threshold plus a seed vote.

| mode | seed threshold | checkpoints that must agree | short-contig threshold | held out: known elements covered / candidates per genome |
|---|---|---|---|---|
| `exploratory` | 0.90 | 1 (every call kept; `support` shows the count) | 0.95 | 0.984 / 37.4 |
| **`default`** | 0.95 | 3 | 0.95 | 0.989 / 12.4 |
| `conservative` | 0.995 | 3 | 0.99 | 0.903 / 5.0 |

- **The extend threshold** follows the seed threshold, capped at 0.99.
- **The figures** are from the [held-out test](lockbox2.md): 150 genomes from unseen genera, scored once.
- **The vote needs three checkpoints.** With fewer (`--checkpoint`), it votes with as many as there are, and with one it does not apply. A warning says so, because the thresholds were calibrated with the vote.
- **Short-contig thresholds** are the lowest at which the false-call rate on shredded genomes stays at or below that of full-length windows, with the mode's vote.

## Options

`prophound predict --help` lists every option in three groups. The defaults are the tested configuration: changing an advanced or experimental option gives results the [held-out test](lockbox2.md) does not describe.

### Basic

| option | default | what it does |
|---|---|---|
| `--input PATH` | | a genome file; repeat it, or give a directory, for a batch |
| `--input-list FILE` | | a file listing genomes, one path per line (relative to the file; `#` comments) |
| `--outdir DIR` | required | where outputs go; a directory already holding outputs is refused |
| `--force` | off | replace the outputs in `--outdir` (every file prophound writes there is removed first) |
| `--mode` | `default` | `exploratory`, `default` or `conservative` |
| `--device` | `auto` | `auto` (CUDA, then Apple MPS, then Intel XPU, then CPU), `cpu`, `cuda`, `rocm`, `mps`, `xpu` (`intel` is the same) |
| `--write-region-fasta` | off | also write the calls' sequences to `regions.fna` |
| `-q` / `-v` | | quieter / more verbose |

### Advanced: thresholds, model, input handling

| option | default | what it does |
|---|---|---|
| `--checkpoint FILE ...` | the bundled three | other model checkpoints; weights-only files load with no flag |
| `--config FILE` | the bundled `model.yaml` | the model config for `--checkpoint` files |
| `--trust-checkpoint` | off | allow full training checkpoints, which load by unpickling arbitrary Python objects; only for files you trust |
| `--threshold`, `--seed-threshold`, `--extend-threshold` | the mode's | override the window scores needed to start and to extend a call |
| `--min-seed-support K` | the mode's | keep a call only if at least K checkpoints make it; 1 keeps every call |
| `--write-seed-scores` | off | also write each checkpoint's window scores to `windows.tsv`, so any vote can be re-applied offline |
| `--ambiguity-policy` | `impute` | how IUPAC codes reach the model; `impute` replaces each with a fixed, reproducible base and changes nothing on clean sequence |
| `--batch-size N` | 128 | windows scored at once; lower it if the GPU runs out of memory |
| `--num-workers N` | 0 | extra processes reading windows |
| `--coverage-frac F` | 0.1 | the share of covering windows that must pass for a position to be inside a call |
| `--max-gap-bp N`, `--bridge-prob F` | 50000, 0.1 | merge nearby calls across gaps that stay element-like |
| `--edge-snap-frac F` | 0.5 | move edges within this fraction of a window of a contig end onto it |
| `--short-contigs` | `scan` | `scan` contigs of 3–10 kb with shorter windows, or `skip` them |
| `--short-contig-threshold F` | calibrated | override the short-contig threshold |
| `--no-domain-check` | off | skip the non-prokaryote check and its warning |
| `--att-search` | `on` | search each call for an att site |
| `--composition-edges` | `auto` | move estimated ends to composition change-points (`auto`: every mode, with a looser gate in Exploratory) |
| `--refine-window-len N` | 5000 | the fine pass's window; 0 turns refinement off |
| `--refine-stride`, `--refine-threshold`, `--refine-margin` | 250, 0.6, one window | the fine pass's step, threshold and the flank it re-scores |
| `--edge-alpha F`, `--edge-alpha-ramp LO HI` | 0.6, 15000 60000 | where refined edges sit between the fine core (0) and the coarse edge (1), scaled with element size |
| `--split-min-gap-bp`, `--split-max-prob` | 4000, 0.05 | the shortest and the emptiest fine-score gap that can split a call |

### Experimental: not validated for default use

| option | default | what it does |
|---|---|---|
| `--ensemble` | `max` | combine the checkpoints' scores by `max` or `mean`; `mean` costs recall |
| `--grid-offset N` | 0 | shift the window grid, to test how calls depend on window placement |
| `--window-len`, `--stride` | 10000, 1000 (as trained) | change the window; the mode thresholds were set at the trained length, and a warning says so |

## Input handling

- **Ambiguity codes.** Two ambiguous bases in a row could drop a window's score from 1 to 0 under the training-time encoding. The default `impute` replaces each with one concrete base, drawn with a seed taken from the window's own content, so a window always encodes the same way and clean sequence is unaffected. `n-channel` gives the training-time behavior.
- **Short contigs.** Scanned by default with shorter windows, which makes draft assemblies and MAGs usable. Contigs under 3 kb are never scanned, and the warning gives their count and total length.
- **Contig ends.** On fragmented input an element cut by a contig break is scored only partly near the break, so an edge near a contig end is moved onto it, and `contig_edge` marks every call touching one: the element may continue on another contig.
- **Gaps.** Calls at least half N are set aside to `regions.filtered.tsv`.
- **Both strands.** The model reads both strands identically, so a contig and its reverse complement give the same window scores. With ambiguity codes they can differ slightly, since each strand's windows draw their own replacement bases.

### The non-prokaryote warning

prophound was trained on bacteria, and on eukaryotic sequence it calls large blocks. Before scoring, it reads the first and last 500 bp of every contig of at least 10 kb, and the input's total length. It warns when:

- two or more contig ends carry a known eukaryotic telomere repeat in tandem (vertebrate and fungal TTAGGG, plant TTTAGGG, insect TTAGG, yeast and *Plasmodium* variants); or
- the input totals over 17 Mb, more than any known prokaryotic genome (the largest is about 16 Mb).

Neither fired on 166 bacteria (including linear *Streptomyces* and *Borrelia* replicons and genomes up to 16 Mb) or 5 archaea. One or the other fired on 8 of 12 eukaryotes. It missed compact eukaryotic genomes under 17 Mb without assembled telomeres (*S. pombe*, *Candida albicans*, the microsporidian *Encephalitozoon*) and a 10 Mb slice of human sequence. It looks at the whole input, so eukaryotic contigs mixed into a bacterial assembly are not detected. A metagenome or several genomes in one file can exceed 17 Mb; `--no-domain-check` skips the check. The check reads only contig ends, so it adds no measurable time ([release text](release-text-0.1.0.md), input outside bacteria).

## The bundled model

- **What ships.** `prophound/models/prophage/model/` holds three weights-only checkpoints (5.5 MB each), the model config `model.yaml`, and `manifest.json` with the SHA-256 of each file. The manifest is checked at every run, and `run_metadata.json` records the hashes used.
- **Safety.** Checkpoints load with `torch.load(weights_only=True)`, which cannot run code. A full training checkpoint, which holds pickled Python objects, is refused unless `--trust-checkpoint` is given. `scripts/strip_checkpoint.py` makes a weights-only file from one.
- **Provenance.** The three checkpoints are three training replicates of one configuration, with different seeds, trained on 1,871 genomes. The bundled files hold the same weights, tensor for tensor, as the training checkpoints, and give byte-identical outputs.

## Installation and packaging

- **pixi** (development and the recommended local setup): `pixi install` builds the environment from `pyproject.toml` and `pixi.lock`; `pixi run prophound ...`, `pixi run test` and `pixi run smoke` run the command, the fast tests and the example end to end.
- **pip**: `pip install git+https://github.com/SyrinxInformatics/prophound.git` into Python 3.11 or later. The runtime dependencies are torch, numpy, pandas, pyfaidx, tqdm, pyyaml and biopython.
- **The wheel** holds only the package code and the bundled model (about 15.5 MB). The evaluation scripts, their recorded outputs and the training configs are excluded. The sdist holds the package, tests, the example and the license files.
- **Bioconda.** The recipe is [`bioconda/recipes/prophound/meta.yaml`](../bioconda/recipes/prophound/meta.yaml): noarch Python, built from the release's GitHub archive, with the model inside the package, so nothing else is downloaded. Its tests run `prophound --version`, `prophound predict --help` and the model's checksum verification. CI lints it with `bioconda-utils` and builds it.
- **A conda environment file,** `metadata/conda/environment.yaml`, mirrors the pixi dependencies for users who prefer conda directly.
- **CI** runs on every push and pull request: errors-only pylint, the fast tests and the example on CPU, the commit check and a gitleaks scan over the full history, and the bioconda lint and build.

## Versioning

prophound follows [Semantic Versioning](https://semver.org/). In 0.y.z, the command line, outputs and defaults may change in any minor release; every change is in the [changelog](../CHANGELOG.md). `run_metadata.json` records the version, the command, the settings and the input's SHA-256, so any result can be traced to the code and data that made it.
