# About prophound

What prophound is, how it works, what has been measured, and what it does not claim. The manual (installing and running it) is the [README](../README.md). Every figure here, with its interval and source, is in [`release-text-0.1.0.md`](release-text-0.1.0.md).

**Status: version 0.1.0, not yet released.** This is research software. Read [What it does not claim](#what-it-does-not-claim) before using its calls.

## What it is

`prophound` is a recall-first prophage screen for prokaryotic genomes. It localizes candidate integrated elements directly from nucleotide sequence, with no gene calling or reference database at run time. It is tuned to recover already-known elements and returns everything else as candidates, at a follow-up cost the user chooses. It is built to run alongside annotation-based tools and to surface what they may not call. How this differs from neighbouring tools is in [`related-tools.md`](related-tools.md), and the framing behind it is in [`project-principles.md`](project-principles.md).

The name is a prophage hound, and a homophone of "profound". The project was called `mge-detect` until 2026-09-29. Other mobile-element types may follow the prophage detector, which is why the code keeps a `models/prophage/` layout.

## How it works

A reverse-complement-symmetric CNN scores 10 kb windows along each contig, sliding 1 kb at a time. Three training replicates run as an ensemble, and they ship with the package (16.5 MB in all). Calls are built from the score track: stretches where enough windows score high, merged across dips inside large elements. The model was trained on genomes labelled by geNomad 1.9.4, but it never sees geNomad's inputs (gene calls, marker HMMs), and nothing is looked up at run time.

**Boundaries.**
- **Att site.** Where an element has an att site (a direct repeat left at both ends by integration), `edges` is `att site` and the ends are exact to the base, with no database or gene calling. About 15% of known elements have one. More than half of these att cores sit on a tRNA or tmRNA gene, the classic integration sites, which the search never looks at (`att-boundaries.md` "Are the att sites real?").
- **Estimated.** Elsewhere `edges` is `estimated`: the model's estimate, moved to where the sequence's composition turns from host-like to element-like. The median error is about 2 kb on our test samples, and each end has a range that holds the true end about 90% of the time.

Neither element in the bundled example has an att site. The calls' ends are 1–4 kb from the elements' (58,837–73,073 and 119,962–154,275, 1-based), inside the ranges. [`examples/README.md`](../examples/README.md) walks through the output.

## Modes, measured

Each mode is priced in calls per genome, the follow-up budget. Figures are from the held-out test set: 150 complete RefSeq genomes from genera absent from training, 202 geNomad elements, scored once. They describe the defaults before the 0.1.0 boundary stages, and a second held-out test of the shipped pipeline replaces them ([`lockbox2.md`](lockbox2.md)).

| mode | rule | geNomad elements covered | calls per genome |
|---|---|---|---|
| `exploratory` | window score 0.90; every call kept | 0.965 | 32.4 |
| **`default`** | 0.95; kept only if all three checkpoints make it | **0.960** (95% interval 0.924–0.980) | 11.0 |
| `conservative` | 0.995; all three checkpoints | 0.901 | 4.9 |

"Covered" means calls cover at least half of an element (covered50). On the 70 held-out genomes with no geNomad element, Default still makes 6.8 calls per genome. That is the floor cost of a run.

In Exploratory mode, start triage with `seed_support`: calls only one checkpoint makes are almost never known elements.

## What has been measured

- **Agreement with geNomad, held out:** Default covers 0.960 of 202 elements (0.924–0.980). The test set was preregistered and scored once ([`lockbox.md`](lockbox.md)). Its lower bound is under 0.95, so 95% recall is not established.
- **Expert-curated truth (LAMBDA, 39 genomes):** Default covers 0.933 of 149 curated prophages, against 0.82 for geNomad, PHASTER and PIDE. It makes 1.7–4.1× their calls per genome, though, and per unit of cost it is level with them: after adjusting for chance, its paired differences to all three span zero.
- **Alongside geNomad:** adding prophound Default to geNomad recovers 16 more curated LAMBDA prophages and loses none, at about 15 extra calls per genome. Twelve of the 16 come from two genomes.
- **LAMBDA's own candidates:** prophound has a call on all 19 phages LAMBDA's reviewers judged likely but left out of its curated set, and covers at least half of 16 (0.75 above chance; 18 in the preregistered scoring, before the 0.1.0 boundary stages trimmed two calls). geNomad covers 15. Most of them are called by annotation tools too.
- **Off-target calls** carry mobile-element features 3.1× as often as matched background, mostly conjugation and integrase/transposase genes rather than phage structural genes. Many look like other integrative elements, not missed prophages.

## What it does not claim

- **No validated novel elements.** No call has been confirmed in the lab. Calls are candidates for follow-up.
- **Not more accurate than** geNomad or any other tool. It finds more at a higher cost.
- **Every agreement figure is with geNomad,** which also supplied the training labels. The model inherits its blind spots.
- **Small elements:** recall under 10 kb is low (0.36 on implanted inoviruses at Default) and partly measured on synthetic data.
- **Fragmented input:** contigs of 3–10 kb are scanned, and contigs under 3 kb never are. On MAG-like assemblies Default keeps 0.92–0.94 of its complete-genome recall, which does not beat geNomad.
- **Several elements may fall under one call** (9% of held-out elements).
- **Strand symmetry:** on sequence with IUPAC ambiguity codes, a contig and its reverse complement can give slightly different calls ([`scientific-validity.md`](scientific-validity.md) RC-3).
- **Platforms:** results are verified on Apple M4 (CPU and GPU) and an NVIDIA V100S. Other platforms are expected to agree but are untested.

## Run time

On a MacBook Pro (Apple M4, 16 GB), Default mode, measured before refinement became a default (it adds 10–23% on the Apple GPU; the table is re-run before release):

| input | time | memory |
|---|---|---|
| a 5 Mb genome, Apple GPU (picked automatically) | ~35 s | ~0.5 GB |
| the same genome as a 374-contig assembly | ~19 s | ~0.5 GB |
| a 5 Mb genome, CPU only | ~45 min | ~1.4 GB |

A GPU makes a large difference. On an NVIDIA V100S a median 3.6 Mb genome takes ~27 s. Details: [`benchmarks.md`](benchmarks.md).

## Development

Training and large-scale evaluation run on the University of Texas at San Antonio's HPC cluster (ARC). This repository holds the CLI, the packaging and the evaluation record.

- `src/prophound/` is the package. The prophage detector is in `src/prophound/models/prophage/`, including the bundled model in `model/` and the dev-only evaluation scripts in `scripts/`.
- `docs/` holds the living topic docs; the [README's documentation list](../README.md#documentation) says what each covers. The rules the project works by are in [`project-principles.md`](project-principles.md).
- `pixi run test` runs the fast tests, and `pixi run smoke` runs the example end to end.
