# mge-detect

`mge-detect` is a database-independent CLI tool for detecting mobile genetic elements (MGEs) in prokaryotic sequences. It combines classical sequence-analysis heuristics with machine learning — the ML model does the heavy lifting — so it can flag both recognized MGE types and novel or divergent ones that reference-database tools miss.

## Status

This is research software under active development. The first (and currently only) detector targets **prophages**:

- a reverse-complement CNN scored over sliding windows of nucleotide sequence, with three user-selectable confidence modes (`Exploratory`, `Default`, `Conservative`)
- wired into a working CLI: `mge-detect predict --input <fasta> --outdir <dir> --config <config.yaml> --checkpoint <r1.pt> <r2.pt> <r3.pt>`. Passing the three training replicates runs them as an ensemble, which is the recommended and evaluated configuration.
- evaluated against geNomad ground truth at the region level on a held-out, **element-size-stratified** sample of 42 genomes / 179 elements. Recall counts an element as found only when calls cover at least half of it. With the three-checkpoint ensemble: **96.1% in `Default` (about 11 calls per genome), 98.3% in `Exploratory` (about 29) and 94.4% in `Conservative` (about 8.5)**. Default and Conservative keep a call only when all three checkpoints make it (`docs/cli-and-packaging.md` "Seed vote per mode"). Operating settings were chosen on this sample. On an independent second sample (18 genomes / 52 elements) the same settings give 98.1% / 100% / 96.2% at about 10.5 / 31 / 6.4 calls per genome.
- that second sample earlier sat at 94%, with every miss in one genome (*E. coli* O157:H7 EDL933). The cause was a bug, not the genome's prophage density: IUPAC ambiguity codes in the sequence silently collapsed the model's scores to zero. Ambiguity codes are now replaced by in-distribution bases by default (`--ambiguity-policy impute`, a no-op on clean sequence), which turns the missed 37 kb prophage into a confident call. Contigs of 3–10 kb, which were previously skipped without warning, are now scanned with calibrated thresholds. Both matter for draft assemblies and MAGs
- region boundaries were the weak point and are much improved: a fragment-merging stage takes 60 kb+ elements from median IoU 0.38 to 0.85 (median predicted length 33 kb → 70 kb, against 70.6 kb of truth) and removes fragmentation almost entirely, while *reducing* the number of calls. See `docs/gap-filling.md`
- **agreement with geNomad (about 28% of `Default` ensemble calls overlap a geNomad element) should be read as a wet-lab cost proxy, not an error rate.** This is a discovery tool: a call a database-driven tool does not make is the kind of output it exists to produce. Telling genuine discoveries from genuine mistakes among those calls is open work. One signal is already reported: each ensemble call carries `seed_support`, the number of checkpoints that make it on their own, and calls only one checkpoint makes are almost never known elements (`docs/error-analysis.md`)

### What is not yet established

The tool's central claim — that it flags *novel* elements a database tool would miss — is **not yet tested, and the first release does not make it**. Every number above measures agreement with geNomad, which also supplied 100% of the training labels. Two early checks are in (`docs/scientific-validity.md`). First, recall holds on elements with no close relative in training: 91% of test elements have none at the nucleotide level, and the ensemble overlaps 97.5% of them with at least one call. Nucleotide clustering cannot see protein-level homology, though. Second, nucleotide composition alone gets nowhere near the model (AUPRC ≤0.16 against ~0.60), so the model uses more than composition. *Which* biology it uses is still unmeasured. Treat this as research software whose agreement with geNomad is characterized and whose novelty claim is not.

Model training and large-scale evaluation run on the University of Texas at San Antonio's HPC cluster (ARC); this repository is where the CLI, packaging, and portable inference path are developed. See [`SESSIONS.md`](SESSIONS.md) for current progress and [`CLAUDE.md`](CLAUDE.md) for where everything is documented.

## Repository layout

- `src/mge_detect/` — the `mge-detect` package and CLI entrypoint
- `src/mge_detect/models/prophage/` — the prophage detector: model/training/inference code, dev-only validation scripts, and its own detailed research notes
- `CLAUDE.md` — orientation for AI coding agents working in this repo
- `SESSIONS.md` — running log of what's been done and what's next

## Acknowledgements

This work received computational support from UTSA's high-performance computing cluster operated by Tech Solutions.
