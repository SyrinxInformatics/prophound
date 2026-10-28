# Golden regression set

Two small genomes from the second lockbox, with the `regions.tsv` and `regions.detailed.tsv` that prophound 0.1.0 makes on them (Default mode, CPU; made 2026-10-28, and every column matches the second lockbox's V100 run `851599`). Between them they exercise most of the call pipeline. The example genome (`examples/`) adds ambiguity codes.

| genome | covers |
|---|---|
| `GCF_053667655.1` (*Holzapfeliella* sp. JNUCC 80, 1.34 Mb) | 7 contigs: an att-site call, 3 contig-edge calls, 2 short-contig calls, composition edges; 8 calls set aside by the vote |
| `GCF_052828685.1` (*Daquilactobacillus linshuiensis*, 1.40 Mb) | 2 split calls; composition edges on 6 of 7 calls; 9 calls set aside by the vote |

The FASTAs are not stored here. `genomes.tsv` gives each one's NCBI URL and the SHA-256 of the uncompressed file.

```
pixi run python scripts/golden_check.py fetch --fasta_dir <dir>
pixi run python scripts/golden_check.py check --fasta_dir <dir> --out <dir>                # CPU: byte for byte
pixi run python scripts/golden_check.py check --fasta_dir <dir> --out <dir> --device mps   # other devices: calls only
```

On CPU the outputs must match byte for byte. On other devices only the calls are compared, because scores differ in the last digits between devices (`docs/benchmarks.md`). The check ends with `[golden] PASS` or `[golden] FAIL`. The CPU check takes about 25 minutes on an Apple M4.
