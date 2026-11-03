# Project principles

The framing behind prophound and the rules its development and evaluation follow. The other docs apply them.

## What prophound is for

**prophound is a discovery tool.** It is built to work without any database, and to be reliable enough that a region it flags which no database contains is worth taking to the lab. It is not built for maximum precision or for agreement with existing databases. Four things follow for how its results are read:

- **Calls outside known elements are the product.** A tool that agreed perfectly with geNomad would score 100% precision and add nothing new. The extra calls are where anything new would be.
- **Precision against geNomad is a cost, not an accuracy.** It is reported as candidates per genome, the follow-up budget a mode asks for. Read as accuracy, it would count every potential discovery as an error.
- **The recall target is a floor.** Recovering known elements shows the tool does not miss what is already findable. It says nothing about whether it finds new ones.
- **Novelty is untested.** It is the largest claim a tool like this could make, and no measurement so far tests it, so the release makes no novelty claim ([release text](release-text-0.1.0.md)).

**Exploratory mode puts recall first.** Its job is proposing candidates a person can confirm in the lab. A missed element cannot be recovered later, while a false candidate costs one check. So the first release establishes a consistent, well-measured baseline, without tuning the model for precision. Expanding the training data and re-tuning come after it.

## The acceptance criterion

**At least 95% of known prophage regions recovered, with extra calls acceptable.** Recovery is measured as covered50: calls cover at least half of an element. Truth is geNomad's calls on complete RefSeq genomes. Every figure is reported with its interval and with the role of its sample (tuning, confirmation or held out; see [Research integrity](research-integrity.md)). Whether a held-out test establishes the criterion is decided by a rule fixed before scoring ([the held-out test](lockbox2.md)).

## Changing a method

A new or alternative method is adopted only after it is compared with what ships and shown not to hurt:

1. **Measure it on the metrics that select it and on the costs those metrics cannot see.** Covered50 and false calls per Mb count elements and calls, not size, so also check:
   - base pairs called, and base pairs called off known elements;
   - boundary overlap (IoU) and candidates per genome;
   - complete genomes and genomes with no prophage;
   - every mode, on both the tuning and the confirmation samples.
2. **If anything gets worse, look for a compromise first:** a smaller parameter, a narrower scope, or a combination with the existing path. The trade-off table decides.
3. **Never choose by an endpoint that cannot see the cost.** Covered50, for example, rewards longer calls.
4. **After implementing, confirm that the shipped code reproduces the measured variant on the same data,** including every offline scorer that rebuilds calls, so the evaluation measures the code that ships.

**The shipped modules are hash-checked.** A change to them must reproduce a fixed set of reference outputs byte for byte on CPU (`scripts/reference_outputs.py compare`: 30 cases covering every mode and option), or explain each difference. Each such change after a held-out test is recorded as an amendment to that test ([the held-out test](lockbox2.md), "Amendments").

## Standing checks

- **Both strands, and gaps.** Every analysis considers the reverse complement. The model and the steps after it must cope with runs of N and with IUPAC ambiguity codes; `scripts/rc_n_resilience.py` asserts the window-level properties.
- **Three replicates per training change,** each with its own seed. A single training run falls inside the spread between seeds.
- **Small jobs run locally;** the computing cluster is for large data, heavy compute and large downloads.

## Records and wording

- **A preregistration is never edited.** Anything changed after it is written as a dated amendment beneath it, so a reader can check what was fixed in advance ([the held-out test](lockbox2.md); the LAMBDA checks in [Claims and evidence](claims-and-evidence.md)).
- **Wording.** Never "first" and never "guarantee". Write "may not call" rather than "does not call", and call the output "candidates for follow-up". Agreement with geNomad is reported as agreement, not accuracy ([Claims and evidence](claims-and-evidence.md)).
