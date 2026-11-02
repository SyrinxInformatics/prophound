# Project principles

The framing and the working rules the other docs cite. Living doc, updated in place. Each rule carries the date it was set; the user sets them, and later docs apply them.

## Product framing

**prophound is a discovery tool** (2026-09-21). It is not designed for absolute precision, or to align reliably with existing databases. It is designed to be database-independent, and reliable enough that a region it flags which no database has is worth taking to the wet lab. Four things follow for how results are read:
- **Off-target calls are the product, not the error term.** A tool that agreed perfectly with geNomad would score 100% precision and have no discovery value.
- **Precision against geNomad is a cost proxy**, best expressed as calls per genome (the follow-up budget). It is not an accuracy figure, and headlining it as one would mislead.
- **The recall target is a floor.** Recovering known elements means the tool does not miss what is already findable. It says nothing about novelty.
- **Novelty is the biggest and most exposed claim, and nothing measured so far tests it** (`docs/scientific-validity.md`). The release makes no novelty claim (`docs/release-text-0.1.0.md`).

**Exploratory mode favors recall** (2026-09-15). Its job is proposing potentially novel sequences that a person can confirm in the lab. A missed novel element cannot be recovered, while a false positive costs one wasted check. So recall comes before precision, and a checkpoint is not tuned for precision before the first release. The aim is a scientifically consistent baseline first, then a larger dataset and re-tuning (`docs/retraining-infrastructure.md`).

## Acceptance criterion

At least 95% region recall on known regions, with extra calls acceptable (2026-09-17). Recall is measured as covered50 against geNomad truth, and is reported with intervals and the sample's role (`docs/research-integrity.md`). Whether a held-out set establishes it is decided by a rule fixed in advance (`docs/lockbox.md`, `docs/lockbox2.md`).

## Method-change procedure

Before adopting a new or alternative method (2026-10-05), compare it against what ships and confirm it does not hurt. If it does, look for a hybrid or compromise first.
1. **Measure it on the metrics that select it and on the costs those metrics cannot see.** Covered50 and false calls per Mb count elements and calls, not size, so also check:
   - bp called and bp off-element;
   - boundary IoU and calls per genome;
   - complete-genome and phage-free regressions;
   - every mode, on both the tuning and the check samples.
2. **If anything gets worse, sweep for a compromise:** a smaller parameter, a narrower scope, or combining it with the existing path. Put the trade-off table to the user.
3. **Never choose by an endpoint that cannot see the cost.** For example, covered50 rewards longer calls.
4. **After implementing, confirm that the shipped code reproduces the measured variant on the same data.** Check every offline scorer that rebuilds calls, so the evaluation actually uses the new code.

Worked example: edge snap (`docs/claims-and-evidence.md` "Amendment 2026-10-05"). Any change to the hash-checked modules must also reproduce the CPU reference outputs byte for byte (`scripts/reference_outputs.py compare`) or explain each difference, and it gets a lockbox amendment.

## Standing checks

- **Reverse complement and N regions** (2026-09-22). Every analysis considers the reverse complement, and the model and the pipeline after it must be resilient to N regions and ambiguity codes (`docs/scientific-validity.md` "Input-robustness standing check").
- **Replicates** (2026-09-17). Three replicates per training-config change, each with its own `run_name` and seed. Single runs fall inside the seed variance (`docs/model-reengineering.md`).
- **Local first** (2026-09-29). Small jobs run locally. The cluster is for large data, heavy compute and large downloads (`docs/environment.md`).

## Records

- **Preregistrations are never edited.** Every change after one is written as a dated amendment (`docs/lockbox.md`, `docs/lockbox2.md`, the LAMBDA checks in `docs/claims-and-evidence.md`).
- **`src/prophound/models/prophage/notes/` is a frozen historical snapshot** of the earlier development cycle. It is not edited, so where it conflicts with `docs/`, `docs/` is current.
- **Wording** (`docs/claims-and-evidence.md`): never "first" and never "guarantee". Write "may not call" rather than "do not call", and call the output "candidates for follow-up".
