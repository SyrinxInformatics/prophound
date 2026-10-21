# Example input

`hin_rd_kw20_region.fna` is 220 kb of *Haemophilus influenzae* Rd KW20 (RefSeq `NC_000907.1`, positions 1,440,001–1,660,000). It contains the genome's two prophage regions that geNomad 1.9.4 calls, with flanking host sequence on both sides. It also carries 25 IUPAC ambiguity codes, as some assemblies do. Those are handled by the default `--ambiguity-policy impute`.

## Run it

```bash
prophound predict --input examples/hin_rd_kw20_region.fna --outdir example_out
```

This takes a few seconds on a GPU (Apple M4: 1.6 s) and about 2 minutes on CPU only (`--device cpu`). It uses the bundled model, Default mode.

## Expected output

The run ends by listing the candidates, best first. `regions.tsv` holds them as a table:

| rank | call | start–end (1-based) | length | edges | start range | end range | geNomad element (1-based) | support |
|---|---|---|---|---|---|---|---|---|
| 1 | `prophound_2` | 116,738–156,000 | 39.3 kb | estimated | 114,579–136,134 | 136,604–158,159 | 119,962–154,275 (34.3 kb) | 3/3 |
| 2 | `prophound_1` | 53,843–81,000 | 27.2 kb | estimated | 52,349–67,259 | 67,584–82,494 | 58,837–73,073 (14.2 kb) | 3/3 |

How to read it (the full guide is [`docs/interpreting-results.md`](../docs/interpreting-results.md)):

- **Both elements are found, and nothing else is called.** Each call covers its whole element, so both count as found under covered50, the release's recall metric (at least half of an element's length covered).
- **The ends are estimates** (`edges` = `estimated`): neither element has a findable att site, a direct repeat left at both ends by integration, which about 15% of known elements do. The calls overshoot the elements by 3–8 kb here. The start and end ranges say where each true end lies, about 90% of the time, and both elements' ends fall inside them. Both calls' ends sit on round numbers (81,000 and 156,000), and `note` says why: the fine pass could not move them, because the score stays high up to them. Those are the least certain ends; here they overshoot by 7.9 and 1.7 kb. Where an att site is found, `edges` is `att site` and the ends are exact to the base.
- **`support` 3/3** means all three bundled models make the call on their own. Default mode keeps only those. In Exploratory mode, calls with support 1/3 are the first to discount.
- **`rank`** puts the call most like known prophages first. Here both score 1.0 at their best window, and the larger one has the higher average.

Files: `expected_regions.tsv` is `regions.tsv` from a CPU run, `expected_regions.detailed.tsv` is `regions.detailed.tsv` (every column, 0-based like BED), and `expected_calls.tsv` holds the columns the smoke test compares.

On another device the window scores can differ in the last digits, because floating-point order differs. The calls are the same: CPU and Apple GPU agree here, and on all 30 reference cases (`docs/benchmarks.md`).

## Use it as a smoke test

```bash
pixi run smoke        # or: pytest -m slow tests/test_example.py
```

This runs the example and checks the calls against `expected_calls.tsv`. Every column, including scores, is in `expected_regions.tsv` for a byte-level check on CPU.
