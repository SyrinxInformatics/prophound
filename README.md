# mge-detect

`mge-detect` is a database-independent CLI tool for detecting mobile genetic elements (MGEs) in prokaryotic sequences. It combines classical sequence-analysis heuristics with machine learning — the ML model does the heavy lifting — so it can flag both recognized MGE types and novel or divergent ones that reference-database tools miss.

## Status

This is research software under active development. The first (and currently only) detector targets **prophages**:

- a reverse-complement CNN scored over sliding windows of nucleotide sequence
- validated on held-out test data with three user-selectable confidence modes (`Exploratory`, `Default`, `Conservative`)
- not yet wired into the `mge-detect` CLI — see [`src/mge_detect/models/prophage/notes/V1_INFERENCE_CONTRACT.md`](src/mge_detect/models/prophage/notes/V1_INFERENCE_CONTRACT.md) for the target interface and [`SESSIONS.md`](SESSIONS.md) for current progress

Model training and large-scale evaluation run on the University of Texas at San Antonio's HPC cluster (ARC); this repository is where the CLI, packaging, and portable inference path are developed.

## Repository layout

- `src/mge_detect/` — the `mge-detect` package and CLI entrypoint
- `src/mge_detect/models/prophage/` — the prophage detector: model/training/inference code, dev-only validation scripts, and its own detailed research notes
- `CLAUDE.md` — orientation for AI coding agents working in this repo
- `SESSIONS.md` — running log of what's been done and what's next

## Acknowledgements

This work received computational support from UTSA's high-performance computing cluster operated by Tech Solutions.
