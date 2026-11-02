# How to read prophound's results

prophound looks at a bacterial genome's DNA and marks regions that look like **integrated prophages**: phage genomes sitting in the chromosome. It uses no reference database and no gene annotation, so it can flag elements databases don't contain. Every region it reports is a **candidate to check**.

This guide is for reading the results. Every file and column is defined in [`docs/output-files.md`](output-files.md).

## Start with `regions.tsv`

One row per candidate, best first. Coordinates are **1-based and inclusive**, as GenBank, GFF3 and genome browsers show them: `start` is the region's first base and `end` its last.

| column | what it tells you |
|---|---|
| `rank` | the order to look at candidates in: 1 is the most like known prophages (see "Rank and support") |
| `call_id` | the candidate's name, the same in every output file; numbered by position |
| `contig` | the sequence it is on |
| `start`, `end` | where it is |
| `length_kb` | its length |
| `edges` | how exactly its ends are known: `att site` or `estimated` (next section) |
| `start_range`, `end_range` | where each true end lies: a single number for an att-site edge, a range for an estimated one |
| `support` | how many of prophound's three models find it on their own (`3/3` is strongest) |
| `note` | anything unusual: an att site's repeat, a region running off a contig end, a lot of N |

## The edges

### `att site`: exact to the base, if the repeat is real

When a phage integrates, its attachment site recombines with the host's, leaving the prophage between two copies of the same short sequence, the att core (attL on one side, attR on the other; typically tens of bp). When prophound finds such a direct repeat, one copy near each end, it reports the region from the **first base of the left copy to the last base of the right copy**. `note` gives the repeat's length, and `regions.detailed.tsv` its sequence (`att_core`).

- **How often it is wrong.** A repeat that long can occur by chance. From the rate on host sequence with no prophage, about **12–17% of att-site edges in Default mode are chance repeats** (4–6% in Conservative, 27–32% in Exploratory; development samples). Treat an att-site edge as a strong lead.
- **How often it is found.** About 15% of known prophages have an att pair prophound can see. Many are degraded or have a core too short to tell from chance, so a candidate without one is no less likely to be a prophage.
- **Against curated annotations,** att-site edges sat a median 22 bp from published coordinates (small literature sample), mostly differing by exactly one core copy.
- **A convention to know.** prophound's region includes both core copies; some annotations leave one out. To match them, start at `start + core length` or end at `end − core length` (the core length is in `note`).

### `estimated`: the best available guess

Most candidates have no att site. Their ends are placed in two steps:

1. **The model's score** along the genome. It changes gradually over several kb around a real end, so it places an end only roughly, usually a few kb past the element.
2. **The sequence's composition.** Prophages usually differ from their host in base composition, so each end is moved to where the sequence shifts from host-like to element-like, as tools for genomic islands do. This halved the typical error, to about **2 kb** (median), and lost no element in Default mode. Exploratory mode accepts smaller composition differences, so it moves ends more often, at the cost of one partially detected element in our tests.

`start_range` and `end_range` say where the true ends lie. Each true end fell inside its range in **about 78% of cases on held-out genomes** (the second lockbox) and 87–91% on development samples, so an end outside its range does happen. The ranges are often 10–30 kb wide, mostly inside the region: the typical error is small, but a minority of ends are far off, mostly where only part of an element was detected. Two practical readings:

- **Sequence almost certainly inside the element** lies between the upper end of `start_range` and the lower end of `end_range`.
- **To be sure a primer pair or cloned fragment spans the whole element,** place it outside the lower end of `start_range` and the upper end of `end_range`.

**Ends on a round number** (a multiple of 1,000) are ends nothing could place: the score stays high right up to them, and composition did not move them. `note` says so ("end not refined"). They are the least certain estimated ends; their ranges held less often (about 70–85%). They are about 1 in 20 ends in Default mode and under 1 in 100 in Exploratory.

A gene-based look at the region (an integrase gene, a tRNA gene at one end, phage structural genes) narrows an estimated edge far more than the range does.

## Rank and support

- **`rank`** sorts candidates by score (highest window score, then average). In our tests this put matches to known prophages first, so **a truly novel element can rank low**. If you are looking for something new, read past the top of the list. A low-ranked att-site candidate is worth a look: the repeat is evidence of integration independent of the score.
- **`support`** counts the models that find the candidate independently. Default and Conservative mode keep only `3/3` candidates. A candidate split from a longer one shows the support of the call it came from. In Exploratory mode, `1/3` candidates very rarely match known elements, so triage them last.

## What most candidates are

On 150 held-out genomes, every Default candidate that touches no known element was sorted by what its genes and sequence show:

| | share | what it means |
|---|---|---|
| phage-like | about 4% | two or more phage structural genes |
| other mobile element | about 24% | conjugation, plasmid or virus hallmark genes, an integrase or transposase, or an att site |
| island of unknown function | about 12% | composition unlike the rest of the genome, nothing else |
| no detectable signal | about 60% | none of the above |

- **Expect most candidates outside known prophages to show no detectable signal.** That is the price of finding as much as possible. Some may still be real, since these checks only recognize what databases already know.
- **Check long candidates first.** Under 10 kb, about 7 in 10 show no signal; at 60 kb or more, about 6 in 10 carry a mobile-element or phage signal.
- **Input should be bacterial.** On human sequence, over half of a 10 Mb slice was called, often in large blocks. Archaeal genomes are untested, and on some a third or more of the genome is called. Remove host or eukaryotic contigs from a metagenome or assembly first.
- **Plasmids are often called whole.** On held-out genomes, a third of plasmids had one call covering most of the plasmid. Such a call describes the plasmid, which may itself be mobile, and locates nothing inside it. Check the contig before reading it as a prophage.
- **A known false alarm:** clusters of ribosomal-protein genes, whose unusual composition can look element-like. Set aside a candidate made mostly of them.

Details: [`docs/extra-calls.md`](extra-calls.md).

## Choosing a mode

The mode sets how many candidates you get, and so the cost of following them up. Counts are from 150 held-out genomes.

| mode | candidates per genome (typical) | use it when |
|---|---|---|
| `exploratory` | ~37 | you would rather check many candidates than miss one |
| `default` | ~12 | most uses |
| `conservative` | ~5 | follow-up is expensive and you want only the strongest |

## What a candidate does not tell you

- **A high score is not evidence by itself.** Genomes shuffled so that no genes are left still get calls, and on high-GC genomes most of the shuffled sequence is called. Real bacterial genomes show no such effect, but the score may partly reflect how unlike ordinary host sequence a stretch is (untested). Look at what a candidate contains before acting on it.
- **It is not proof.** We have confirmed no candidate in the lab. A candidate can be another kind of mobile element (an integrative conjugative element, a genomic island) or a false call.
- **What prophound agrees with.** Its tests count agreement with geNomad, a widely used annotation tool that also supplied the training labels, plus a curated benchmark (LAMBDA). It inherits some of geNomad's blind spots.
- **Small elements** (under about 10 kb) are often missed.
- **Fragmented assemblies.** Contigs under 3 kb are never scanned. A candidate running off a contig end (`note`) may continue on another contig.
- **Adjacent prophages** can come out as one candidate.

## Checking a candidate in the lab

These are standard approaches, listed to show how prophound's outputs feed into them.

- **Att-site candidates.** If the prophage can excise, PCR with primers just outside both ends (outside `start` and `end`) gives a short product from the empty host site (attB) in cells that lost it. Primers inside the element near each end, facing outward, give a product only from the circular excised phage (attP). Sequencing either product shows the core and confirms the boundary. Inducing the culture (for example with mitomycin C) raises the excised fraction.
- **Estimated candidates.** Annotate the region first (genes, integrase, tRNA) to place the ends, then design junction primers as above. Primers placed from the ranges alone can miss.
- **Long candidates** (roughly 60 kb or more) may hold several elements side by side. Before designing junction primers, look for a dip in the score track inside the call and for more than one integrase gene. If there are several elements, test each one's junctions.
- **Any candidate.** Plot the score track around it (`windows.tsv`), look at its genes, and compare with annotation tools. A candidate none of them calls, carrying phage-like genes, is what prophound exists to surface.

## Other files

| file | when you need it |
|---|---|
| `regions.detailed.tsv` | every column behind `regions.tsv` (scores, the model's own edges, the att core), 0-based half-open like BED |
| `regions.bed`, `regions.gff3` | to view candidates in a genome browser (IGV, JBrowse) next to your annotation |
| `regions.filtered.tsv` | candidates a filter set aside, with the reason: where to look if you think something was missed |
| `windows.tsv`, `fine_windows.tsv` | the score along the genome, for plotting |
| `run_metadata.json` | settings, model and input checksums: keep it with your results |
