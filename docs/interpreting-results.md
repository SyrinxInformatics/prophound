# How to read prophound's results

prophound looks at a bacterial genome's DNA and marks regions that look like **integrated prophages**: phage genomes sitting in the chromosome. It uses no reference database and no gene annotation, so it can flag elements that databases don't contain. Every region it reports is a **candidate to check**, not a confirmed prophage.

This guide is for reading the results. Every file and column is defined in `docs/output-files.md`.

## Start with `regions.tsv`

One row per candidate, best first. Coordinates are **1-based and inclusive**, as GenBank, GFF3 and genome browsers show them: `start` is the first base of the region and `end` its last.

| column | what it tells you |
|---|---|
| `rank` | the order to look at candidates in: 1 is the one most like known prophages (see "Rank and support") |
| `call_id` | the candidate's name, the same in every output file. Numbered by position, not by rank |
| `contig` | the sequence it is on |
| `start`, `end` | where it is |
| `length_kb` | its length |
| `edges` | how exactly its ends are known: `att site` or `estimated` (next section) |
| `start_range`, `end_range` | where each true end lies. A single number for an att-site edge; a range for an estimated one |
| `support` | how many of prophound's three models find it on their own (`3/3` is strongest) |
| `note` | anything unusual: an att site's repeat, a region running off the end of a contig, a lot of N |

## The edges

### `att site`: exact to the base, if the repeat is real

When a phage integrates, the host's and the phage's attachment sites recombine, leaving the prophage between two copies of the same short sequence, the att core (attL on one side, attR on the other; typically tens of bp). When prophound finds such a direct repeat, one copy near each end of a candidate, it reports the region from the **first base of the left copy to the last base of the right copy**. `note` gives the repeat's length, and `regions.detailed.tsv` its sequence (`att_core`).

- **How often it is wrong.** A repeat that long can turn up by chance. On host sequence with no prophage, the search accepts one at a rate that implies about **12–17% of att-site edges in Default mode are chance repeats** (4–6% in Conservative, 27–32% in Exploratory). These figures come from development samples. Treat an att-site edge as a strong lead, not a fact.
- **How often it is found.** About 15% of known prophages have an att pair prophound can see. Many prophages are degraded, or have a core too short to tell from chance. A candidate without one is not therefore less likely to be a prophage.
- **Against curated annotations,** att-site edges were a median 22 bp from published coordinates, on a small literature sample. Most differences were exactly one core copy (next point).
- **A convention to know.** prophound's region includes both core copies. Some annotations leave one out. To match them, start at `start + core length`, or end at `end − core length` (the core length is in `note`).

### `estimated`: the best available guess

Most candidates have no att site. Their ends are placed in two steps:

1. **The model's score** along the genome. It changes gradually over several kb around a real end, so on its own it places an end only roughly, and usually a few kb past the element.
2. **The sequence's composition.** Prophages usually differ from their host in base composition, so each end is moved to where the sequence turns from host-like to element-like, the way tools for genomic islands find edges. In our tests this halved the typical error, to about **2 kb** (median), without losing any element in Default mode. Exploratory mode applies it a little more often (it accepts host and element compositions that differ less), at the cost of one partially detected element in our tests.

`start_range` and `end_range` say where the true ends lie: each true end fell inside its range in **about 78% of cases on held-out genomes** (the second lockbox), and in 87–91% on the development samples the ranges were fitted to or checked on. So a true end outside its range is not rare. The ranges are wide, often 10–30 kb: the typical error is small, but a minority of ends are far off, mostly where only part of an element was detected. Most of the width lies inside the region. Two practical readings:

- **Sequence that is almost certainly inside the element** lies between the upper end of `start_range` and the lower end of `end_range`.
- **To be sure a primer pair or a cloned fragment spans the whole element,** place it outside the lower end of `start_range` and the upper end of `end_range`.

**Ends on a round number** (a multiple of 1,000) are ends nothing could place: the model's score stays high right up to them, and composition did not move them either. `note` says so ("end not refined"). They are the least certain estimated ends; their ranges held less often in our tests (about 70–85%). In Default mode they are about 1 in 20 ends; in Exploratory mode, under 1 in 100.

A gene-based look at the region (an integrase gene, a tRNA gene at one end, phage structural genes) narrows an estimated edge much further than the range does.

## Rank and support

- **`rank`** sorts candidates by how strongly the model scores them (its highest window score, then its average). In our tests this put candidates that match known prophages first. That also means **a truly novel element can rank low**. If you are looking for something new, read past the top of the list. A low-ranked att-site candidate is worth a look: the repeat is evidence of integration that does not depend on the score.
- **`support`** counts the models that find the candidate independently. Default and Conservative mode keep only `3/3` candidates. A candidate split from a longer one shows the support of the call it was split from. In Exploratory mode, `1/3` candidates very rarely match known elements: start triage with those last.

## What most candidates are

On 150 held-out genomes, we sorted every Default candidate that touches no known element by what its genes and sequence show:

| | share | what it means |
|---|---|---|
| phage-like | about 4% | two or more phage structural genes |
| other mobile element | about 24% | conjugation, plasmid or virus hallmark genes, an integrase or transposase, or an att site |
| island of unknown function | about 12% | composition unlike the rest of the genome, nothing else |
| no detectable signal | about 60% | none of the above |

- **Expect most candidates outside known prophages to show no detectable signal.** That is the price of finding as much as possible. Some of those may still be real, because these checks only recognize what databases already know.
- **Long candidates are worth checking first.** Under 10 kb, about 7 in 10 show no signal. At 60 kb or more, about 6 in 10 carry a mobile-element or phage signal.
- **A known false alarm:** clusters of ribosomal-protein genes, whose unusual composition can look element-like. If a candidate is mostly ribosomal-protein genes, set it aside.

Details: `docs/extra-calls.md`.

## Choosing a mode

The mode sets how many candidates you get, which is the cost of following them up.

| mode | candidates per genome (typical) | use it when |
|---|---|---|
| `exploratory` | ~30 | you would rather check many candidates than miss one |
| `default` | ~11 | most uses |
| `conservative` | ~5 | follow-up is expensive and you want only the strongest |

## What a candidate does not tell you

- **It is not proof.** No candidate has been confirmed in a lab by us. A candidate can be another kind of mobile element (an integrative conjugative element, a genomic island) or a false call.
- **What prophound agrees with.** Its tests count agreement with geNomad, a widely used annotation tool that also supplied its training labels, plus a curated benchmark (LAMBDA). It inherits some of geNomad's blind spots.
- **Small elements** (under about 10 kb) are often missed.
- **Fragmented assemblies.** On contigs under 3 kb nothing is scanned. A candidate that runs off the end of a contig (`note`) may continue on another one.
- **Several prophages next to each other** can come out as one candidate.

## Checking a candidate in the lab

These are standard approaches, listed to show how prophound's outputs feed into them.

- **Att-site candidates.** If the prophage can excise, a PCR with primers just outside both ends (outside `start` and `end`) gives a short product from the empty host site (attB) in the cells that lost it. Primers inside the element near each end, facing outward, give a product only from the circular excised phage (attP). Sequencing either product shows the core and confirms the boundary. Inducing the culture (for example with mitomycin C) raises the excised fraction.
- **Estimated candidates.** Annotate the region first (genes, integrase, tRNA) to place the ends, then design junction primers as above. Primers placed from the ranges alone can miss.
- **Any candidate.** Look at the score track around it (`windows.tsv`, plotted along the contig) and at the genes it contains, and compare with annotation tools. A candidate none of them calls, carrying phage-like genes, is the kind prophound exists to surface.

## Other files

| file | when you need it |
|---|---|
| `regions.detailed.tsv` | every column behind `regions.tsv` (scores, the model's own edges, the att core), 0-based half-open like BED |
| `regions.bed`, `regions.gff3` | to view candidates in a genome browser (IGV, JBrowse) next to your annotation |
| `regions.filtered.tsv` | candidates set aside by a filter, with the reason: a place to look if you think something was missed |
| `windows.tsv`, `fine_windows.tsv` | the score along the genome, for plotting |
| `run_metadata.json` | the settings, model and input checksums: keep it with your results |
