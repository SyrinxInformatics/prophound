# mge-detect

`mge-detect` is a database-independent CLI tool for detecting mobile genetic elements (MGEs) in prokaryotic sequences. It combines classical sequence-analysis heuristics with machine learning — the ML model does the heavy lifting — so it can flag both recognized MGE types and novel or divergent ones that reference-database tools miss.

## Status

This is research software under active development. The first (and currently only) detector targets **prophages**:

- a reverse-complement CNN scored over sliding windows of nucleotide sequence, with three user-selectable confidence modes (`Exploratory`, `Default`, `Conservative`)
- wired into a working CLI: `mge-detect predict --input <fasta> --outdir <dir> --config <config.yaml> --checkpoint <best.pt>`
- evaluated against independent geNomad ground truth (not just held-out window-level metrics) on a held-out, **element-size-stratified** sample of 42 genomes / 179 elements, across three training seeds: **96–97% region recall in `Default` mode and ~99% in `Exploratory`**, at roughly 11 and 21 calls per genome respectively. A second sample of 18 genomes / 52 regions, which played no part in tuning, sits at 94% — every miss on it falls in a single unusually prophage-dense genome (*E. coli* O157:H7 EDL933)
- region boundaries were the weak point and are much improved: a fragment-merging stage takes 60 kb+ elements from median IoU 0.38 to 0.85 (median predicted length 33 kb → 70 kb, against 70.6 kb of truth) and removes fragmentation almost entirely, while *reducing* the number of calls. See `docs/gap-filling.md`
- **precision against geNomad (~35–42%) should be read as a wet-lab cost proxy, not an error rate.** This is a discovery tool: a call a database-driven tool does not make is the kind of output it exists to produce. Distinguishing genuine discoveries from genuine mistakes among those calls is open work, not a solved problem

### What is not yet established

The tool's central claim — that it flags *novel* elements a database tool would miss — is **not yet tested**. Every number above measures agreement with geNomad, which also supplied 100% of the training labels, and the train/test split is by genome rather than by phage homology, so a test genome's prophage may be a close relative of a training one. Whether the model keys on phage biology (structural genes, integrases) or on a correlate such as nucleotide composition is likewise unmeasured. `docs/scientific-validity.md` sets out the program that addresses both; the first experiments are running. Treat this as research software whose accuracy is characterized and whose novelty claim is not.

Model training and large-scale evaluation run on the University of Texas at San Antonio's HPC cluster (ARC); this repository is where the CLI, packaging, and portable inference path are developed. See [`SESSIONS.md`](SESSIONS.md) for current progress and [`CLAUDE.md`](CLAUDE.md) for where everything is documented.

## Repository layout

- `src/mge_detect/` — the `mge-detect` package and CLI entrypoint
- `src/mge_detect/models/prophage/` — the prophage detector: model/training/inference code, dev-only validation scripts, and its own detailed research notes
- `CLAUDE.md` — orientation for AI coding agents working in this repo
- `SESSIONS.md` — running log of what's been done and what's next

## Acknowledgements

This work received computational support from UTSA's high-performance computing cluster operated by Tech Solutions.
