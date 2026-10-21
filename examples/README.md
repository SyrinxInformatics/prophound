# Example input

`hin_rd_kw20_region.fna` is 220 kb of *Haemophilus influenzae* Rd KW20 (RefSeq `NC_000907.1`, positions 1,440,001–1,660,000). It contains the genome's two prophage regions that geNomad 1.9.4 calls, with flanking host sequence on both sides. It also carries 25 IUPAC ambiguity codes, as some assemblies do. Those are handled by the default `--ambiguity-policy impute`.

## Run it

```bash
prophound predict --input examples/hin_rd_kw20_region.fna --outdir example_out
```

This takes a few seconds on a GPU (Apple M4: 1.6 s) and about 2 minutes on CPU only (`--device cpu`). It uses the bundled model, Default mode.

## Expected output

`expected_calls.tsv` holds the calls. `expected_regions.tsv` is the full `regions.tsv` from a CPU run.

| call | start–end (0-based, half-open) | length | `boundary` | start range | end range | geNomad element (same coordinates) | seed_support |
|---|---|---|---|---|---|---|---|
| `prophound_1` | 53,842–81,000 | 27.2 kb | model | 52,348–67,258 | 67,584–82,494 | 58,836–73,073 (14.2 kb) | 3 |
| `prophound_2` | 116,737–156,000 | 39.3 kb | model | 114,578–136,133 | 136,604–158,159 | 119,961–154,275 (34.3 kb) | 3 |

How to read the calls:

- **Both elements are found, and nothing else is called.** Each call covers its whole element, so both count as found under covered50, the release's recall metric (at least half of an element's length covered).
- **These boundaries are model estimates** (`boundary` = `model`): neither element has a findable att site, a direct repeat left at both edges by integration, which about 15% of known elements do. The fine pass places the edges, which still overshoot the elements by 3–8 kb here. The start and end ranges say where each true edge lies, 90% of the time; both elements' edges fall inside them. Where an att site is found, `boundary` is `att` and the edges are exact to the base.
- **`seed_support` 3** means all three bundled checkpoints make the call on their own. Default mode keeps only those. In Exploratory mode, calls with support 1 are the first to discount.
- `mean_prob` (0.68 and 0.70) is the mean window score across the call. It ranks calls; it is not a probability that the call is a prophage.

On another device the window scores can differ in the last digits, because floating-point order differs. The calls are the same: CPU and Apple GPU agree here, and on all 30 reference cases (`docs/benchmarks.md`).

## Use it as a smoke test

```bash
pixi run smoke        # or: pytest -m slow tests/test_example.py
```

This runs the example and checks the calls against `expected_calls.tsv`. Every column, including scores, is in `expected_regions.tsv` for a byte-level check on CPU.
