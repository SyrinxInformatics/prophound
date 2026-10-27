# Benchmarks: run time, memory and cross-device agreement

Written 2026-10-19. What a user should expect when they run prophound, measured on a current laptop (user, 2026-10-19: "this laptop has hardware that most other users will likely have similar or slightly more / less computing power, so users know what they are in for"). The cluster GPU figures are secondary: most users will not have a V100.

The first run (2026-10-19) measured the code at `b4848b9`: the lockbox code plus the MPS pooling fix (`docs/tech-debt-audit.md` S10). Phase 3 of the refactor targeted speed. **The single-genome Apple GPU rows were re-run after it on 2026-10-20** ("Phase 3 speed work" below), and the release text quotes them. The batch row and the CPU figures are from the first run; Phase 3 leaves CPU time unchanged (profiled).

## The machine

| | |
|---|---|
| model | MacBook Pro, Mac16,1 |
| chip | Apple M4: 10 cores (4 performance + 6 efficiency), integrated GPU |
| memory | 16 GB |
| OS | macOS 27.0.1, native arm64 (not Rosetta) |
| power | mains (AC attached), otherwise idle |
| software | torch 2.13.0 from conda-forge (pixi), Python 3.14, 4 torch CPU threads |
| prophound | `b4848b9`, Default mode, 3-checkpoint ensemble, shipped defaults |

## Headline: Apple M4 GPU (`--device mps`)

Protocol: one warm-up run, then 3 timed runs per input, each a fresh `prophound predict` process. Wall time includes start-up and model loading; peak memory is the process's maximum resident set (`/usr/bin/time -l`). Every repeat wrote byte-identical `regions.tsv`. Scripts and raw results: `scripts/benchmark.py`, `analysis/benchmarks/mps_{typical,batch,profile}.json`.

| input | size | median wall time (range) | per Mb | peak memory |
|---|---|---|---|---|
| typical complete genome (`GCF_000146165.2`) | 4.97 Mb, 1 contig | **33.5 s** (33.3–33.6), after Phase 3; was 43.5 s | 6.7 s | 512 MB |
| the same genome as a MAG-like shred (`mag_s0`) | 4.95 Mb, 374 contigs | **19.3 s** (19.3–19.3), after Phase 3; was 25.6 s | 3.9 s | 538 MB |
| 10 genomes spanning the corpus size range (before Phase 3, as 10 separate runs; not re-measured) | 1.8–7.0 Mb, 48.9 Mb total | **7.4 min** in total; median 46.5 s per genome | 8.9 s, plus 1.0 s fixed per genome | ≤620 MB |

**In plain terms: about 35 seconds for a typical 5 Mb bacterial genome, about 7 seconds per Mb, and about half a GB of memory.** A batch of 100 genomes of that size is roughly 75 minutes.

Fragmented assemblies are faster per Mb than complete genomes. Fewer full 10 kb windows fit, and fragments under 3 kb are never scanned.

**Where the time goes** (typical genome, profiled in process, `mps_profile.json`):

| stage | time | share |
|---|---|---|
| model forward pass on the GPU (3 checkpoints × 4,960 windows; 2.2 ms per window per checkpoint) | 33.3 s | 76% |
| reading and one-hot encoding windows on the CPU (0.68 ms per window, done once per checkpoint, so 3×) | 10.2 s | 23% |
| loading checkpoints, moving batches to the GPU | 0.2 s | <1% |

So the GPU does most of the work, and encoding is the one cheap target. A vectorised encoder (13× faster on its own, `docs/tech-debt-audit.md` S9), run once rather than once per checkpoint, would save about 9 s, to roughly 34 s per 5 Mb genome. That is Phase 3.

### Phase 3 speed work (2026-10-20)

Two changes, both with byte-identical outputs (`docs/tech-debt-audit.md` Phase 3): a vectorised one-hot encoder, and one scoring pass that reads and encodes each window once and runs it through all three checkpoints. Before, every window was encoded once per checkpoint. Measured back to back on the same idle machine, timing the code before Phase 3 (`f3d526e`) and after it. Protocol as above; 1 warm-up plus 3 timed runs. `analysis/benchmarks/p3_{before,after}_{mps,mps_profile,cpu_profile}.json`.

| | before | after | change |
|---|---|---|---|
| typical genome, MPS wall time (median) | 43.2 s | **33.5 s** | −22% |
| MAG shred, MPS wall time (median) | 24.9 s | **19.3 s** | −22% |
| encoding, as the CLI does it (profile) | 10.3 s (0.69 ms per window, 3 passes) | **0.49 s** (0.099 ms, 1 pass) | −95% |
| GPU forward pass (profile) | 33.2 s | 32.5 s | unchanged |
| peak memory | 511 MB | 512 MB | unchanged |

Every run before and after wrote the same `regions.tsv`. The forward pass is now about 97% of the profiled time, so further MPS gains need a faster model or a different device. Peak memory was already about 100 MB below the 2026-10-19 figure (610 MB) before Phase 3. The likely cause is that checkpoints are now loaded once (UI pass, step 1).

**CPU (300-window profile, 4 threads):** the forward pass costs 186–193 ms per window per checkpoint, and encoding 0.1–0.7 ms. Encoding was about 0.3% of the CPU run time, so Phase 3 barely changes it: projected from the profile, a typical genome still takes about 45 minutes. The open CPU questions (pip's torch build with oneDNN, more threads) are untouched.

## CPU only (`--device cpu`)

For users without an Apple Silicon or NVIDIA GPU. Protocol: 2 timed runs of the typical genome, no warm-up; a warm-up is small next to the run time.

Measured 2026-10-19 evening on the same machine, idle and on mains power; code identical to the MPS run (`cpu_typical.json`).

| input | size | wall time (2 runs) | per Mb | peak memory |
|---|---|---|---|---|
| typical complete genome (`GCF_000146165.2`) | 4.97 Mb | **45.0 min** (2,701.8 and 2,700.3 s) | 9.1 min | 1.37 GB |

**In plain terms: on CPU alone, about 45 minutes for a typical 5 Mb genome, 62× slower than the same laptop's GPU, and about twice the memory.** Both runs wrote byte-identical `regions.tsv`. Against the MPS run, the 20 calls are identical (coordinates, `seed_support`, contig-edge flags), with scores within 6e-7. The first attempt (afternoon) was stopped at 16:00 before finishing, and nothing from it is used.

Earlier measurement on the same machine, for scale: `GCF_000027305.1` (1.83 Mb) took 1,053 s, 17.6 minutes (about 9.6 minutes per Mb) on 4 threads. Two causes, neither fixed yet:
- **The torch build.** The pixi/conda-forge build is generic arm64 with oneDNN off (`torch.backends.mkldnn.is_available()` is False).
- **The model.** Its 10 kb convolutions run twice per window (forward and reverse complement), for three checkpoints. Measured in process: about 0.2 s per window per checkpoint on CPU.

Whether pip's torch wheels or more threads help is still open: Phase 3 did not change CPU time (see the CPU profile below).

## Cluster GPU: NVIDIA V100S (ARC, secondary)

From the lockbox run (`841208`, `analysis/lockbox/prophound/wall_seconds.tsv`): 150 genomes × 3 modes, one process per genome and mode, including start-up.

| | |
|---|---|
| genome sizes | median 3.6 Mb (0.23–13.4) |
| median wall time per genome | 27 s (90th percentile 42 s); the same in every mode, since all modes score the same windows |
| linear fit | 4.5 s fixed + 6.1 s per Mb |

The M4 GPU is about 1.45× slower per Mb than a V100S.

## Against other tools on the same genomes

The 10 benchmark genomes are G5 genomes, so the G5 timings apply to exactly these inputs (`analysis/g5/{jaeger,genomad}_wall_seconds.tsv`, `docs/claims-and-evidence.md` "G5"). Medians over the 10:

| tool | hardware | median per genome (~5.1 Mb) |
|---|---|---|
| **prophound** (3 checkpoints, Default) | this M4, GPU | **46.5 s** |
| Jaeger 1.26.5, `57341_1.5M` (the G5 primary model) | this M4, CPU, 8 workers | 206.5 s |
| Jaeger 1.26.5, `default` model | this M4, CPU, 8 workers | 55.5 s |
| geNomad 1.9.4 | ARC CPU node, 4 threads (different hardware) | 342 s |

Caveats:
- **Jaeger ran on CPU here** because its TensorFlow Metal plugin fails to load under TF 2.21 (`metadata/g5-comparators/pixi.toml`). So this is not Jaeger's best case on a Mac.
- **geNomad was timed on different hardware** and with a 1.47 GB database, so its row is context, not a comparison. The G6 comparison needs both tools on one machine (release text, "Not yet measured").
- The Jaeger runs date from 2026-09-30/10-01, not today.

## Cross-device agreement

**Outputs are not byte-identical across devices, and cannot be in general.** Each device's convolution kernels accumulate in a different order, so scores differ in the last few significant figures. What was measured:

| comparison | data | largest score difference | calls |
|---|---|---|---|
| CPU vs CUDA (V100) | 200 windows × 3 checkpoints from G5b, recorded CUDA scores | 8e-6 | no window crosses 0.95 |
| CPU vs MPS, after the pooling fix | Phase 0 reference matrix, 30 cases | 1.7e-5 | identical in 30 of 30 cases, with identical messages |
| CPU vs MPS, before the fix | the same | 0.34 | different in 8 of 29 cases |

A call can differ between devices only when a decision value lands within about 1e-5 of a threshold. At a margin of 1e-4 (6× the largest difference seen), 289 of 536,280 window scores in the 45 G5b runs sit near a threshold (0.90 / 0.95 / 0.99 / 0.995), as do 42 of 153,583 bin means near the 0.1 bridge gate. Most runs have a few such values. Whether any of them decides a call is the next measurement:

**Measured 2026-10-19 (`scripts/cross_device_g5b.py`, `analysis/benchmarks/cross_device_g5b_mps.json`).** All 45 G5b runs (9 genomes: complete plus four shredded conditions; inputs MD5-identical to ARC's) were re-run on the M4 GPU and compared with the scores recorded on the V100S. The comparison covers 134,070 windows × 4 scores (ensemble and three checkpoints) and every bin mean.
- **Largest score difference:** 1.9e-5.
- **Decision flips:** one. Checkpoint 1 scored one window 0.98999995 on CUDA and 0.99000007 on MPS, 1.2e-7 either side of Conservative's 0.99 extend threshold. No bin mean flipped at the bridge gate.
- **Effect of that flip:** rebuilding the calls from both tracks with the shipped builder and vote, it extends checkpoint 1's own Conservative call by 1 kb. **The released calls are identical in every mode**: a single checkpoint's call only votes, by overlap.
- **CPU vs MPS on the 4.97 Mb benchmark genome:** 20 of 20 calls identical, scores within 6e-7.

So the calls are the same on all three verified devices, and the measured rate of even an internal decision differing is about 1 in 500,000 scores.

**What can be said about other platforms.** Verified so far: Apple M4 CPU, Apple M4 GPU (MPS, after the fix) and NVIDIA V100S (CUDA). Not verified:
- x86 CPUs (Intel and AMD, with MKL or oneDNN);
- other NVIDIA generations;
- AMD ROCm and Intel XPU;
- Windows.

The S10 bug shows why "it should be the same" is not evidence. Proposed coverage is in the next section.

## Proposed cross-platform coverage (pending the user)

1. **ARC, one GPU job and one CPU job.** Run `reference_outputs.py run` on the Phase 0 inputs:
   - on a V100 (CUDA);
   - on an x86 CPU node: Intel or AMD, Linux, the conda-forge or pip torch the ARC env has.

   Then `compare --across-devices` each against the M4 CPU references. This covers CUDA and x86 CPU on identical inputs. The CPU job is slow, so a 6-case subset is enough there.
2. **CI on Linux x86** once the model is bundled (§5): a small subset on every release candidate, against stored references with `--across-devices`.
3. **Release wording until then:** "Verified to give the same calls on Apple M4 (CPU and GPU) and NVIDIA V100S; scores agree to about 1e-5 across these devices; other platforms are expected to agree but are untested."

## Re-running these

```
pixi run python src/prophound/models/prophage/scripts/benchmark.py timing --code . --device mps \
    --input typical=<genome.fna> --workdir <scratch> --out analysis/benchmarks/<name>.json
pixi run python src/prophound/models/prophage/scripts/benchmark.py profile --code . --device mps \
    --input <genome.fna> --out analysis/benchmarks/<name>_profile.json
```

Needs the three checkpoints under `src/prophound/models/prophage/runs/` (gitignored). Run on a quiet machine on mains power. Every result records the machine, power source, torch build and code commit.
