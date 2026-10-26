# Changelog

All notable changes to prophound are recorded here. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow [Semantic Versioning](https://semver.org/). A change to the bundled model, a mode's thresholds or any other default that changes calls is a minor version at least, and is listed under **Calls changed**.

prophound is in initial development (0.y.z). Until 1.0.0 the command line, output files and defaults may change in any minor release; 1.0.0 will mark the first interface we commit to keeping stable.

## [0.1.0] — unreleased

The first release. Earlier development versions were called `mge-detect`.

### Model and calls

- A reverse-complement-symmetric CNN scores 10 kb windows (1 kb stride). The bundled model is an ensemble of three training replicates (`prophage_v2_pool`, weights only, SHA-256s in `model/manifest.json`), trained on geNomad 1.9.4 labels (DB 1.9).
- Calls are built from the score track with no gene calling: a coverage rule, merging of fragments of large elements across score dips, and dropping calls that are at least half N.
- Three modes, priced in calls per genome on 150 held-out genomes:
  - Exploratory: seed threshold 0.90, every call kept, 32.4 calls per genome.
  - Default: 0.95, kept only when all three checkpoints make the call, 11.0 calls per genome.
  - Conservative: 0.995, all three checkpoints, 4.9 calls per genome.
- Contigs of 3–10 kb are scanned with their own thresholds. Contigs under 3 kb are never scanned, and are reported.
- A call within half a window of a contig end is extended to the end and flagged in `contig_edge`.
- IUPAC ambiguity codes are replaced by fixed, reproducible bases (`--ambiguity-policy impute`). This changes nothing on sequence without them. The training-time encoding made two adjacent ambiguous bases drop a window's score to zero.
- **Boundaries.** Where a call has an att site (a direct repeat with one copy at each edge, left by integration), its edges are those two copies, exact to the base (`boundary` = `att`). No database or gene calling is used. Elsewhere, a fine pass (5 kb windows on a 250 bp step) places the edges (`boundary` = `model`; `--refine-window-len 0` turns it off), and each edge is then moved to where the sequence's composition turns from host-like to element-like, which roughly halves the median edge error (`--composition-edges`; Exploratory accepts profiles that differ less). Every call gives a range for each edge (`start_lo` … `end_hi`): exact for att edges, a 90% range for model edges. About 15% of known elements have a findable att site. `docs/att-boundaries.md`.

### Input and output

- Input can be FASTA, GenBank, EMBL or GFF3 with embedded sequence, optionally gzip, bzip2 or xz compressed. The format is read from the content.
- Batch runs: repeat `--input`, give a directory, or use `--input-list`.
- Outputs: `regions.tsv` (the candidates for a researcher: ranked, plain columns, 1-based), `regions.detailed.tsv` (every column, 0-based), `regions.bed`, `regions.gff3`, `windows.tsv` and `run_metadata.json`. The terminal lists the ranked candidates, and `docs/interpreting-results.md` explains how to read them. Calls a filter set aside go to `regions.filtered.tsv`, so a discovery user can recover them without a re-run. `regions.fna` is optional. The full contract is in `docs/output-files.md`.
- An outdir that already holds outputs is refused unless `--force`, and outputs are written atomically.
- Checkpoints load with `weights_only=True`. A full training checkpoint needs `--trust-checkpoint`.
- Devices: CUDA, Apple MPS, Intel XPU or CPU, picked automatically.

### Fixed during development

- Apple MPS gave wrong scores, through adaptive max pooling, until 2026-10-19. CUDA and CPU were unaffected.
- IUPAC ambiguity codes silently collapsed scores to zero, fixed 2026-09-24.
- Contigs under 10 kb were skipped without a warning. They are now scanned down to 3 kb.

### Known limitations

These are in the release text (`docs/release-text-0.1.0.md`). In short:

- Every agreement figure is measured against geNomad, which also supplied the training labels.
- Boundaries are coarser than gene-based tools'.
- Recall on elements under 10 kb is low and partly measured on synthetic implants.
- No call has been validated in the lab, and the release makes no claim to novel elements.
