# About prophound

What prophound is, how it works, what has been measured, and what it does not claim. Installing and running it is in the [README](../README.md); every figure here, with its interval and source, is in [`release-text-0.1.0.md`](release-text-0.1.0.md).

**Status: version 0.1.0, not yet released.** This is research software. Read [What it does not claim](#what-it-does-not-claim) before using its calls.

## What it is

`prophound` is a recall-first prophage screen for bacterial genomes. It locates candidate integrated elements directly from nucleotide sequence, with no gene calling or reference database. It is tuned to recover known elements and returns everything else as candidates, at a follow-up cost the user chooses. It runs alongside annotation tools to surface what they may not call. [`related-tools.md`](related-tools.md) compares it with neighboring tools, and [`project-principles.md`](project-principles.md) gives the framing.

The name is a prophage hound, and a homophone of "profound". The project was called `mge-detect` until 2026-09-29. Other mobile-element types may follow, hence the `models/prophage/` layout.

## How it works

A convolutional network that reads both strands identically scores 10 kb windows along each contig, 1 kb apart. Three training replicates run as an ensemble and ship with the package (16.5 MB). Calls are runs of high-scoring windows, merged across dips inside large elements. The model was trained on genomes labeled by geNomad 1.9.4. It never sees geNomad's inputs (gene calls, marker HMMs), and nothing is looked up at run time.

**Boundaries.**

- **Att site.** Where an element is flanked by an att site (a direct repeat left by integration), `edges` is `att site` and the ends are exact to the base. About 15% of known elements have one. Over half of these att cores sit on a tRNA or tmRNA gene, the classic integration sites, which the search never looks at ([`att-boundaries.md`](att-boundaries.md) "Are the att sites real?").
- **Estimated.** Elsewhere `edges` is `estimated`: the model's edge, moved to where the composition shifts from host-like to element-like. The median error is about 2 kb, and each end has a range. The ranges held the true end in 88–92% of cases on development samples and about 78% on held-out genomes, so treat them as a guide.

Neither element in the bundled example has an att site. The calls' ends are 1–4 kb from the elements' (58,837–73,073 and 119,962–154,275, 1-based), inside the ranges; [`examples/README.md`](../examples/README.md) walks through the output.

## Modes, measured

Each mode is priced in calls per genome, the follow-up budget. Figures are from the held-out test of the shipped pipeline: 150 complete RefSeq genomes from genera absent from training and every earlier test, 186 geNomad elements, preregistered and scored once ([`lockbox2.md`](lockbox2.md)).

| mode | rule | geNomad elements covered | calls per genome |
|---|---|---|---|
| `exploratory` | window score 0.90; every call kept | 0.984 | 37.4 |
| **`default`** | 0.95; kept only if all three checkpoints make it | **0.989** (95% interval 0.962–0.997) | 12.4 |
| `conservative` | 0.995; all three checkpoints | 0.903 | 5.0 |

"Covered" means calls cover at least half of an element (covered50). On the 77 held-out genomes with no geNomad element, Default still makes 9.2 calls per genome: the floor cost of a run. Exploratory touches every known element, yet covers one fewer than Default by half, because its looser composition gate trims two partial detections.

In Exploratory mode, start triage with `seed_support`: calls only one checkpoint makes are almost never known elements.

## What has been measured

Every comparison with other tools is in [`tool-comparisons.md`](tool-comparisons.md).

- **Agreement with geNomad, held out:** Default covers 0.989 of 186 elements (0.962–0.997) on a preregistered test scored once ([`lockbox2.md`](lockbox2.md)). The lower bound clears 0.95, so by the rule fixed in advance, at least 95% of known elements are covered held out. This is agreement with geNomad, not accuracy. The pre-boundary defaults read 0.960 (0.924–0.980) on 150 other genomes ([`lockbox.md`](lockbox.md)).
- **Expert-curated truth (LAMBDA, 39 genomes):** Default covers 0.933 of 149 curated prophages; geNomad, PHASTER and PIDE about 0.82. It makes 1.7–4.1× their calls per genome, and after adjusting for chance its paired differences to all three span zero.
- **Alongside geNomad:** adding Default's calls to geNomad's recovers 16 more curated LAMBDA prophages and loses none, at about 15 extra calls per genome. Twelve of the 16 come from two genomes.
- **LAMBDA's own candidates:** prophound touches all 19 phages LAMBDA's reviewers judged likely but left out of its curated set, and covers at least half of 16 (0.75 above chance; 18 in the preregistered scoring, before the boundary stages trimmed two calls). geNomad covers 15, and annotation tools call most of them.
- **What the extra calls are made of** (held out): about 4 in 10 Default calls outside known elements carry a mobile-element or phage signal (1 in 10 for random host sequence of the same length); 6 in 10 show none. Long calls carry signal far more often (release text §8b).
- **Off-target calls** carry mobile-element features 3.1× as often as matched background, mostly conjugation and integrase/transposase genes. Many look like other integrative elements.

## What it does not claim

- **No validated novel elements.** No call has been confirmed in the lab. Calls are candidates for follow-up.
- **Not more accurate than** geNomad or any other tool. It finds more at a higher cost.
- **Every agreement figure is with geNomad,** which also supplied the training labels. The model inherits its blind spots.
- **Small elements:** recall under 10 kb is low (0.36 on implanted inoviruses in Default) and partly measured on synthetic data.
- **Fragmented input:** contigs of 3–10 kb are scanned; contigs under 3 kb are not. On MAG-like assemblies Default keeps 0.92–0.94 of its complete-genome recall, which does not beat geNomad.
- **Several elements may fall under one call** (5% of held-out elements).
- **Edge ranges are not 90% ranges:** held out, they held the true end in about 78% of cases.
- **Strand symmetry:** on sequence with IUPAC ambiguity codes, a contig and its reverse complement can give slightly different calls ([`scientific-validity.md`](scientific-validity.md) RC-3).
- **Platforms:** verified on Apple M4 (CPU and GPU), an NVIDIA V100S and Linux x86 CPU. Other platforms are expected to agree but are untested.

## Run time

On a MacBook Pro (Apple M4, 16 GB), Default mode, the 0.1.0 defaults (timed 2026-11-02):

| input | time | memory |
|---|---|---|
| a 5 Mb genome, Apple GPU (picked automatically) | under a minute (~48 s) | ~0.5 GB |
| the same genome as a 374-contig assembly | ~25 s | ~0.5 GB |
| a 5 Mb genome, CPU only | ~45 min | ~1.4 GB |

On an NVIDIA V100S a median 3.6 Mb genome takes ~27 s. Details: [`benchmarks.md`](benchmarks.md).

## Development

Training and large-scale evaluation run on the University of Texas at San Antonio's HPC cluster (ARC). This repository holds the CLI, the packaging and the evaluation record.

- `src/prophound/` is the package. The prophage detector is in `src/prophound/models/prophage/`, with the bundled model in `model/` and dev-only evaluation scripts in `scripts/`.
- `docs/` holds the living topic docs; the [README's documentation list](../README.md#documentation) says what each covers. The project's working rules are in [`project-principles.md`](project-principles.md).
- `pixi run test` runs the fast tests, and `pixi run smoke` runs the example end to end.
