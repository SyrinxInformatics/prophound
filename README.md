# mge-detect

`mge-detect` is a database-independent CLI tool for detecting mobile genetic elements (MGEs) in prokaryotic sequences. It combines classical sequence-analysis heuristics with machine learning — the ML model does the heavy lifting — so it can flag both recognized MGE types and novel or divergent ones that reference-database tools miss.

## Status

This is research software under active development. The first (and currently only) detector targets **prophages**:

- a reverse-complement CNN scored over sliding windows of nucleotide sequence, with three user-selectable confidence modes (`Exploratory`, `Default`, `Conservative`)
- wired into a working CLI: `mge-detect predict --input <fasta> --outdir <dir> --config <config.yaml> --checkpoint <best.pt>`
- evaluated against independent geNomad ground truth (not just held-out window-level metrics) on an 18-genome stratified set: the shipped `v1` checkpoint reaches ~83% region recall / ~39% precision, and the best-precision checkpoint measured so far (~77% / ~70%) came from a confounded run whose gains are being isolated — see `docs/model-reengineering.md`. Region boundaries are typically truncated inside the true extent. These are real, open limitations, not yet a finished detector
- model-side re-engineering is ongoing on the HPCC: an isolated `pos_weight`/`amp` ablation is training as of this writing; a boundary-regression head was tried and shelved. Downstream post-processing (boundary truncation, gap-filling, anti-overextension) is scoped in `docs/downstream-postprocessing.md` and developed in parallel

Model training and large-scale evaluation run on the University of Texas at San Antonio's HPC cluster (ARC); this repository is where the CLI, packaging, and portable inference path are developed. See [`SESSIONS.md`](SESSIONS.md) for current progress and [`CLAUDE.md`](CLAUDE.md) for where everything is documented.

## Repository layout

- `src/mge_detect/` — the `mge-detect` package and CLI entrypoint
- `src/mge_detect/models/prophage/` — the prophage detector: model/training/inference code, dev-only validation scripts, and its own detailed research notes
- `CLAUDE.md` — orientation for AI coding agents working in this repo
- `SESSIONS.md` — running log of what's been done and what's next

## Acknowledgements

This work received computational support from UTSA's high-performance computing cluster operated by Tech Solutions.
