# Comparisons with established tools

Every head-to-head comparison of prophound with other prophage tools, with the source of each result. Figures are for the 0.1.0 pipeline unless marked † (measured before the boundary stages were added).

**How to read these comparisons.**

- **prophound makes more calls than the other tools,** and more calls find more by chance. So recall is always given with calls per genome and, where measured, recall above chance (calls rotated to random positions).
- **Truth comes from outside prophound:** expert curation (LAMBDA), the literature, or another tool's calls. Every tool is scored by the same rule (covered50: calls cover at least half of an element).
- **Bars were fixed before the runs** for LAMBDA (G1–G3), G5, G5b and the LAMBDA candidates check ([Claims and evidence](claims-and-evidence.md)).
- **Samples are small** (LAMBDA 39 genomes, literature 7, G5b 9). Read the paired intervals; most span zero.

## Summary

| comparison | tools | result | verdict |
|---|---|---|---|
| recall on expert-curated prophages (LAMBDA, 39 genomes) | geNomad, PHASTER, PIDE, VIBRANT, genomic language models | Default 0.933 against 0.82 for geNomad, PHASTER and PIDE, at 1.7–4.1× their calls per genome | highest raw recall; **not ahead per unit of cost** (paired intervals span zero) |
| run alongside geNomad (G1) | geNomad | +16 curated prophages, −0, at +14.9 calls per genome | **cleared** the preregistered bar (>+5) |
| matched cost against sequence-only models (G2) | EVO2 4k (best genomic language model on LAMBDA) | 0.792 against 0.785 at EVO2 4k's call count | cleared by one element; trails geNomad and PHASTER at their counts |
| literature-named prophages (7 genomes) | geNomad | 0.912 against 0.807; level on intact prophages, ahead on remnants | an example (7 genomes) |
| worked genomes (Sakai, SL1344) | PHASTEST, geNomad | more curated elements covered, at more calls; ends level on Sakai, behind on SL1344 | an example |
| MAG-like fragmented input (G5, G5b) † | geNomad, Jaeger | keeps 0.92–0.94 of its complete-genome recall, against geNomad 0.96–1.02 and Jaeger 0.88–0.99 | **no comparative claim** (G5 refuted, G5b a tie) |
| LAMBDA's reviewed candidates outside its curated set | geNomad, PIDE, PhageBoost | covers 16 of 19 likely phages; geNomad 15 | the annotation tools reach most of them too |
| held-out extra calls with phage head genes | geNomad, Jaeger, PHASTEST | 4 of 14 called by none of the three | candidates for follow-up, none validated |
| run time and footprint | Jaeger, geNomad | M4 GPU 46.4 s per genome against Jaeger 56–207 s (CPU); model 16.5 MB against geNomad's 1.47 GB database | context only (different devices) |

## Curated recall: LAMBDA

LAMBDA's clean subset: 39 genomes with expert-curated prophages, none in prophound's training or validation splits. The other tools' predictions ship with the benchmark and were rescored by the same rule ([Claims and evidence](claims-and-evidence.md) "LAMBDA result", "Amendment 2026-10-23"; release text §2).

| | prophound Default | geNomad | PHASTER | PIDE | EVO2 4k |
|---|---|---|---|---|---|
| covered50, all 39 genomes (149 elements) | **0.933** (0.874–0.983) | 0.822 | 0.824 | 0.824 | 0.785 |
| calls per genome | 14.7 | 4.6 | 3.9 | 8.4 | 4.8 |
| share of each chromosome called | 7.9% | 3.7% | | | |
| recall above chance | **0.822** | 0.770 | 0.772 | 0.776 | |

- On the 34 genomes not in prophound's corpus (the primary row): 0.920 (113 elements; 0.859–0.983).
- **Paired, Default minus geNomad is +0.048 (−0.012 to +0.109)**, minus PHASTER +0.048 (−0.007 to +0.110), minus PIDE +0.044 (−0.020 to +0.113). All span zero: per unit of cost, prophound is not shown to be ahead.
- Of 8 pairings (Default and Conservative against four tools), three exclude zero, two of them against VIBRANT. Uncorrected for eight comparisons, so not claimed.
- Under 10 kb (13 elements) prophound covered 0.846 against geNomad's 0.500 †, but n is small.
- The 8 kb composition ceiling was chosen with LAMBDA in view; for that setting LAMBDA is confirmation after use.

## Run alongside geNomad (G1)

Adding prophound Default's calls to geNomad's recovers **16 more curated LAMBDA prophages and loses none**: +0.110 covered50 on the 37 genomes both ran (paired 0.018–0.210), at +14.9 calls per genome. The preregistered bar was more than +5, the best existing pairing with geNomad. Twelve of the 16 come from two genomes (*Photorhabdus* TTO1, 9; *P. aeruginosa* LESB61, 3); all 16 have containment ≤0.16 against training positives (release text §3).

## Matched cost against sequence-only models (G2)

Exploratory's calls, ranked by score and cut to EVO2 4k's per-genome count, cover **0.792 against EVO2 4k's 0.785**: the preregistered bar is cleared by one element, with no paired separation. Cut to geNomad's count, prophound trails geNomad (0.801 against 0.822), and to PHASTER's, it trails further (0.770 against 0.824). Before the boundary stages these read 0.805, a tie, and 0.797 (release text §4).

## Literature-named prophages

57 elements in 7 genomes, curated from the literature, with geNomad scored the same way ([release text](release-text-0.1.0.md) §5; results in [`analysis/litpilot/`](../prophound/models/prophage/analysis/litpilot/)). The pooled figure is flattered: 34 of the 57 elements are at least half contained in training.

| subset | n | prophound Default | geNomad |
|---|---|---|---|
| all | 57 | 0.912 | 0.807 (paired +6 / −0) |
| intact prophages | 28 | 28 of 28 | 27 of 28 |
| novel to training (containment <0.1) | 13 | 8 of 13 | 6 of 13 |
| remnants | 12 | 7 of 12 | 3 of 12 |

Both tools are near ceiling on intact prophages; prophound's margin is in remnants and cryptic prophages. **Its ends are worse than geNomad's:** median IoU 0.73 / 0.62 (test / val) against 0.89 / 0.71.

## Worked genomes against PHASTEST

prophound, PHASTEST (web API) and, for Sakai, geNomad, scored against author-curated elements (results in [`analysis/phastest_cmp/`](../prophound/models/prophage/analysis/phastest_cmp/)).

- ***E. coli* O157:H7 Sakai** (24 elements; in the corpus's validation split, so a worked example): Default covers 21 of 24 with 26 calls; PHASTEST 18 of 24 with 16; geNomad 17 of 24 with 14. Off the curated set, Default calls 0.30 Mb and PHASTEST 0.10 Mb.
- ***S.* Typhimurium SL1344** (7 elements; never used before): every tool covers 7 of 7. Default makes 17 calls, PHASTEST 9.
- **Ends:** level on Sakai (median end error 1.24 against 1.25 kb), PHASTEST ahead on SL1344 (3.78 against 3.19 kb). Where prophound finds an att pair, its ends are within 0.06 kb. It merges tandem elements into one call more often than PHASTEST.

## Fragmented input (G5, G5b) †

Genomes run complete and shredded to MAG-like contigs (N50 20 kb) through prophound, geNomad and Jaeger ([Claims and evidence](claims-and-evidence.md) "G5 result", "G5b result"; release text §6). Retention is recall on the shreds divided by recall on the complete genomes.

| | prophound Default | geNomad | Jaeger |
|---|---|---|---|
| G5 (LAMBDA genomes) | 0.916 (0.844–0.968) | 1.022 | 0.987 |
| G5b (Philympics, held out, preregistered) | 0.938 (0.889–0.983) | 0.958 | 0.881 |

G5 refuted the comparative claim (paired intervals below zero). G5b is a tie: prophound minus geNomad −0.021 (−0.074 to +0.026), minus Jaeger +0.057 (−0.064 to +0.165). **No comparative claim on fragmented input.** The design fact stands: prophound scans contigs down to 3 kb by default; Jaeger's prophage mode scans contigs over 500 kb by default.

## Candidates outside the truth sets

- **LAMBDA's reviewed candidates** (19 regions its reviewers judged likely phages and left out of the curated set; preregistered): Default covers 16 of 19 by half and touches all 19; geNomad covers 15, PIDE 0.789, PhageBoost 0.778. Most are phages the database tools also call (release text §9).
- **Held-out extra calls with phage genes.** Of 1,679 Default calls on the second lockbox that touch no geNomad element, 69 passed the preregistered phage-protein test, and Jaeger and PHASTEST were run on their genomes. Jaeger calls 9 of the 43 tail-only candidates, consistent with most being phage-like injection systems. Of the 14 with head genes, **4 are called by none of geNomad, Jaeger and PHASTEST**, and PHASTEST calls 3 others "intact" that geNomad and Jaeger missed. These were picked after the results, so they are quoted with their denominator and only as candidates ([release text](release-text-0.1.0.md) §8b; results in [`analysis/candidate_homology/`](../prophound/models/prophage/analysis/candidate_homology/)).

## Run time and footprint

On the 10 G5 benchmark genomes (median ~5.1 Mb; [Benchmarks](benchmarks.md)):

| tool | hardware | median per genome |
|---|---|---|
| prophound 0.1.0, Default | Apple M4, GPU | 46.4 s |
| Jaeger 1.26.5, `default` model | Apple M4, CPU, 8 workers | 55.5 s |
| Jaeger 1.26.5, `57341_1.5M` model | Apple M4, CPU, 8 workers | 206.5 s |
| geNomad 1.9.4 | cluster CPU node, 4 threads | 342 s |

Jaeger ran on CPU (its TensorFlow Metal plugin does not load) and geNomad on different hardware, so these timings are context only. A same-machine comparison with geNomad (G6) has not been run.

Run-time data: prophound's bundled model is 16.5 MB; geNomad's database v1.9 is 1.47 GB unpacked, about 89× more. This counts only alongside the recall above.

## Design

Design differences (database, gene calling, training context, boundaries, minimum contig length, sensitivity settings) are in [Related tools](related-tools.md), with sources. Database-free and nucleotide-level detection, discovery framing and sensitivity settings all exist already (PhageBoost, PIDE, Jaeger, Seeker, LAMBDA's models). No tool of this style was found to combine a recall target measured on a preregistered held-out set with modes priced in follow-up cost.

## Not compared

- geNomad timed on the same machine as prophound (G6).
- PhageBoost and DeepVirFinder were not run head to head. PhageBoost's figures come from LAMBDA's files; DeepVirFinder classifies whole contigs and gives no prophage coordinates.
- Installed-environment size.
- No comparison uses wet-lab-confirmed elements that no tool's training data contained.
