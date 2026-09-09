# mge-detect

`mge-detect` is a database-independent CLI tool for detecting mobile genetic elements (MGEs) in prokaryotic sequences. It combines classical sequence-analysis heuristics with machine learning — the ML model does the heavy lifting — so it can flag both recognized MGE types and novel or divergent ones that reference-database tools miss.

## Status

This is research software under active development. The first (and currently only) detector targets **prophages**:

- a reverse-complement CNN scored over sliding windows of nucleotide sequence, with three user-selectable confidence modes (`Exploratory`, `Default`, `Conservative`)
- wired into a working CLI: `mge-detect predict --input <fasta> --outdir <dir> --config <config.yaml> --checkpoint <best.pt>`
- evaluated against independent geNomad ground truth (not just held-out window-level metrics): region-level recall/precision are currently ~77%/~70%, with region boundaries typically truncated inside the true extent — both are real, open limitations, not yet a finished detector
- model-side re-engineering (an isolated architecture ablation, plus a new boundary-regression head aimed at the truncation problem) is implemented and training on the HPCC as of this writing

Model training and large-scale evaluation run on the University of Texas at San Antonio's HPC cluster (ARC); this repository is where the CLI, packaging, and portable inference path are developed. See [`SESSIONS.md`](SESSIONS.md) for current progress and [`CLAUDE.md`](CLAUDE.md) for where everything is documented.

## Repository layout

- `src/mge_detect/` — the `mge-detect` package and CLI entrypoint
- `src/mge_detect/models/prophage/` — the prophage detector: model/training/inference code, dev-only validation scripts, and its own detailed research notes
- `CLAUDE.md` — orientation for AI coding agents working in this repo
- `SESSIONS.md` — running log of what's been done and what's next

## Acknowledgements

This work received computational support from UTSA's high-performance computing cluster operated by Tech Solutions.
