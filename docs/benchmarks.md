# Benchmarks: run time, memory and cross-device agreement

What to expect when you run prophound: how long a genome takes, how much memory it needs, and whether different hardware gives the same calls. The headline machine is a current laptop, since most users will run prophound on something similar. Cluster GPU figures are secondary.

## In short

| input | Apple M4 GPU | Apple M4 CPU only | NVIDIA V100S |
|---|---|---|---|
| a typical 5 Mb complete genome | **under a minute** (47.6 s), ~0.5 GB | **about an hour** (62 min), ~1.1 GB | ~27 s for a median 3.6 Mb genome |
| the same genome as a 374-contig assembly | ~25 s, ~0.5 GB | | |
| 100 genomes of that size | roughly 80 minutes | | |

The model's forward pass is nearly all of the time, so a GPU matters far more than anything else. The calls are the same on every device tested.

## The machine

| | |
|---|---|
| model | MacBook Pro, Mac16,1 |
| chip | Apple M4: 10 cores (4 performance + 6 efficiency), integrated GPU |
| memory | 16 GB |
| OS | macOS, native arm64 |
| power | mains (AC attached), otherwise idle (load under 2) |
| software | torch 2.13.0 from conda-forge (pixi), Python 3.14, 4 torch CPU threads |
| prophound | 0.1.0 (`447294a`; the result files record it as `a9f1851`, its hash before a history rewrite), Default mode, the bundled three-checkpoint model, every default |

## How it was measured

- **Each run is a fresh `prophound predict` process,** so wall time includes start-up and model loading. Peak memory is the process's maximum resident set (`/usr/bin/time -l`).
- **GPU:** one warm-up run, then three timed runs per input; the median is reported. **CPU:** two timed runs, no warm-up (a warm-up is small next to an hour). **Batch:** ten genomes, one run each.
- **Inputs** are the G5 benchmark genomes, checked against their recorded MD5s ([`analysis/g5/input_md5.tsv`](../prophound/models/prophage/analysis/g5/input_md5.tsv)).
- Every result file records the machine, power source, torch build and code commit ([`analysis/benchmarks/`](../prophound/models/prophage/analysis/benchmarks/), `final_*.json`).
- **Every repeat of an input wrote byte-identical `regions.tsv`.**

## Results

| input | size | median wall time (range) | per Mb | peak memory |
|---|---|---|---|---|
| typical complete genome (`GCF_000146165.2`), Apple GPU | 4.97 Mb, 1 contig | **47.6 s** (46.8–48.7) | 9.6 s | 497 MB |
| the same genome as a MAG-like shred, Apple GPU | 4.95 Mb, 374 contigs | **25.2 s** (25.1–25.2) | 5.1 s | 545 MB |
| 10 genomes spanning the training corpus's size range, Apple GPU | 1.8–7.0 Mb, 48.9 Mb in total | **7.2 min** in total; median 46.4 s per genome | 3.4 s fixed + 8.15 s per Mb | ≤522 MB |
| typical complete genome, CPU only (4 threads) | 4.97 Mb | **62.0 min** (61.7 and 62.3) | 12.5 min | 1.14 GB |

- **Fragmented assemblies are faster per Mb** than complete genomes: fewer full 10 kb windows fit, and contigs under 3 kb are never scanned.
- **CPU takes about 78 times as long as the laptop GPU** and about twice the memory. The conda-forge torch build is generic arm64 with oneDNN off. Whether another torch build or more threads would help is untested.

## Where the time goes

- **The forward pass** is about 97% of the time on the GPU: 3 checkpoints × about 5,000 windows for a 5 Mb genome, each read on both strands. On CPU it costs about 0.19 s per window per checkpoint.
- **Reading and encoding windows** is done once and shared by the three checkpoints; it takes about half a second per genome.
- **The boundary stages** (refinement, att-site search, composition edges) add about 36% on the GPU and 38% on CPU. That was measured against the code without them, on the same machine and genome back to back (35.0 s against 47.6 s; `final_control_p3_mps.json`). Refinement alone is +10–23%.

## Other hardware: NVIDIA V100S

From the 450 runs of the first held-out test on a cluster V100S (150 genomes × 3 modes, one process per genome and mode, including start-up; [`analysis/lockbox/prophound/wall_seconds.tsv`](../prophound/models/prophage/analysis/lockbox/prophound/wall_seconds.tsv)). That code predates the boundary stages, so expect about a third more today.

| | |
|---|---|
| genome sizes | median 3.6 Mb (0.23–13.4) |
| median wall time per genome | 27 s (90th percentile 42 s), the same in every mode, since all modes score the same windows |
| linear fit | 4.5 s fixed + 6.1 s per Mb |

The M4 GPU is about 1.45 times slower per Mb than a V100S.

## Against other tools on the same genomes

The 10 benchmark genomes are the G5 genomes, so the comparison tools' G5 timings apply to exactly these inputs ([`analysis/g5/`](../prophound/models/prophage/analysis/g5/), `{jaeger,genomad}_wall_seconds.tsv`). Medians over the 10:

| tool | hardware | median per genome (~5.1 Mb) |
|---|---|---|
| **prophound 0.1.0** (3 checkpoints, Default) | Apple M4, GPU | **46.4 s** |
| Jaeger 1.26.5, `default` model | Apple M4, CPU, 8 workers | 55.5 s |
| Jaeger 1.26.5, `57341_1.5M` model | Apple M4, CPU, 8 workers | 206.5 s |
| geNomad 1.9.4 | cluster CPU node, 4 threads | 342 s |

- **Jaeger ran on CPU** because its TensorFlow Metal plugin does not load (TF 2.21), so this is not Jaeger's best case on a Mac.
- **geNomad ran on different hardware,** with its 1.47 GB database, so its row is context only. A same-machine comparison with geNomad has not been run.

## Cross-device agreement

**Outputs are not byte-identical across devices, and cannot be in general.** Each device's convolution kernels add numbers in a different order, so scores differ in the last few significant figures. What matters is whether calls differ.

| comparison | data | largest score difference | calls |
|---|---|---|---|
| Apple M4 CPU vs GPU, 0.1.0 | the typical 5 Mb genome | at floating-point level | `regions.tsv`, `regions.detailed.tsv` and `regions.filtered.tsv` byte-identical (21 calls) |
| Apple M4 CPU vs GPU | the 30-case reference matrix (every mode and option) | 1.7e-5 | identical in 30 of 30 cases, with identical messages |
| Apple M4 GPU vs NVIDIA V100S | 45 runs: 9 genomes, complete and four shredded conditions; 134,070 windows × 4 scores | 1.9e-5 | identical in every mode |
| Linux x86 CPU vs NVIDIA V100S | a clean install: the example and one held-out genome | | identical |

- **A call can differ only when a score lands within about 1e-5 of a threshold.** Across the 45 GPU-vs-V100S runs, exactly one internal decision differed: one checkpoint scored one window 0.98999995 on CUDA and 0.99000007 on the Apple GPU, either side of Conservative's 0.99 extend threshold. It lengthened that checkpoint's own call by 1 kb, and the released calls were identical, since a single checkpoint's call only votes. That is about 1 decision in 500,000 scores ([`analysis/benchmarks/cross_device_g5b_mps.json`](../prophound/models/prophage/analysis/benchmarks/cross_device_g5b_mps.json); `scripts/cross_device_g5b.py`).
- **Verified platforms:** Apple M4 CPU and GPU, NVIDIA V100S (CUDA) and Linux x86 CPU.
- **Untested:** other NVIDIA generations, AMD CPUs and GPUs (ROCm), Intel GPUs (XPU) and Windows. They are expected to agree, but on Apple GPUs a pooling bug once gave wrong scores until it was found and fixed (0.1.0 includes the fix; see the [changelog](../CHANGELOG.md)), so "it should be the same" is not evidence.

## Re-running these

```bash
pixi run python prophound/models/prophage/scripts/benchmark.py timing --code . --device mps \
    --input typical=<genome.fna> --workdir <scratch> --out <result>.json
pixi run python prophound/models/prophage/scripts/benchmark.py profile --code . --device mps \
    --input <genome.fna> --out <result>_profile.json
```

`--device cpu --repeats 2 --no-warmup` gives the CPU protocol. The bundled model is used by default. Run on a quiet machine on mains power.
