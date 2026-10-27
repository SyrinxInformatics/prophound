# Output files

What `prophound predict` writes, column by column. This is the current output contract (since the UI pass, step 2, 2026-10-20; the researcher table and `regions.detailed.tsv` since 2026-10-21). The frozen `notes/V1_OUTPUT_FILES.md` describes an older format. Update this file whenever an output changes. **How to read the results** is a separate guide for users: `docs/interpreting-results.md`.

**Coordinates.** `regions.tsv`, the researcher table, is **1-based, inclusive**, as GenBank, GFF3 and genome browsers show positions. Every other TSV and BED is **0-based, half-open** (`start` included, `end` excluded), as in BED. GFF3 and the FASTA headers are 1-based, inclusive, as those formats require. A call's 1-based `start` is its 0-based `start` + 1; `end` is the same number in both.

## Where outputs go

- Everything goes into `--outdir`, which is created if missing.
- **An outdir that already holds prophound outputs is refused.** The error names the next free directory (`out_1`, `out_2`, …). `--force` replaces the outputs, after removing every file prophound can write there, so no file from an earlier run survives. Files that aren't prophound's are left alone.
- **Writes are atomic.** Each file is written under a temporary name in the outdir, and all of them are renamed once every file is complete. A failed or interrupted run leaves no partial set of outputs. Earlier outputs stay intact, even with `--force`, unless the run reaches the final rename.

## Files

| file | written | contents |
|---|---|---|
| `regions.tsv` | always | the calls for a researcher: ranked, plain columns, 1-based |
| `regions.detailed.tsv` | always | the same calls with every column, 0-based; what scripts should read |
| `regions.filtered.tsv` | always | candidates set aside by the seed vote or the N filter, with the reason |
| `regions.bed` | always | the calls as BED6 |
| `regions.gff3` | always | the calls as GFF3 |
| `windows.tsv` | always | every scored window and its score |
| `run_metadata.json` | always | settings, input and checkpoint hashes, counts, timing, platform |
| `regions.fna` | `--write-region-fasta` | the calls' sequences |
| `fine_windows.tsv` | refinement on (the default; not with `--refine-window-len 0`) | the refinement pass's windows |

## Batches

A single `--input` file writes the layout above. A **batch** (`--input` given more than once, a directory, or `--input-list`) writes one directory per genome, `outdir/<genome_id>/`, each holding that genome's files exactly as a single run would write them. The top of the outdir holds:

| file | contents |
|---|---|
| `regions.tsv` | every genome's calls as the researcher table, with a leading `genome` column; ranked within each genome |
| `regions.detailed.tsv`, `regions.filtered.tsv` | every genome's calls (or set-aside candidates), in input order. `genome_id` tells them apart; `call_id` is numbered per genome |
| `regions.bed`, `regions.gff3` | the combined calls, named `genome_id:call_id` so names are unique. Contig IDs can repeat across genomes |
| `batch_summary.tsv` | one row per genome: `genome_id`, `input`, `status` (`ok` or `failed`), `format`, `contigs`, `bp`, `calls`, `filtered`, `elapsed_s`, `error` |
| `run_metadata.json` | `"batch": true`, the shared settings, checkpoint hashes, counts (`n_genomes`, `n_ok`, `n_failed`), one entry per genome, timing, command and platform. Each genome's own `run_metadata.json` has the full detail |

- **Before any work**, a batch checks that every input exists and that genome IDs are unique (the file name without extension; `x.fna` and `x.fa` collide).
- **A directory** contributes its genome files: `.fa`, `.fna`, `.fasta`, `.fas`, `.fsa`, `.gb`, `.gbk`, `.gbff`, `.genbank`, `.embl`, `.gff`, `.gff3`, optionally `.gz`, `.bz2` or `.xz`. It is not recursive, and hidden files are skipped.
- **`--input-list`** takes one path per line, relative to the list file, with `#` comments.
- **A failing genome** is recorded in `batch_summary.tsv` with its error, and the batch continues. The exit code is 1 if any genome failed.
- **Checkpoints are read once** per batch.
- **`--force`** over an earlier batch also removes the prophound files in its genome directories that the new batch does not rewrite, including a genome that fails this time.

## `regions.tsv`

The researcher table (`src/report.py`): one row per call, **sorted by `rank`**, 1-based inclusive.

| column | meaning |
|---|---|
| `rank` | 1 … N within the genome: `max_prob` (4 dp), ties broken by `mean_prob`, then position. On the tuning sample this ordered calls on known elements first better than `mean_prob` alone, in every mode, and it held on the check sample (`docs/att-boundaries.md` "Ranking") |
| `call_id` | as in `regions.detailed.tsv` (numbered by position) |
| `contig` | the sequence ID (`contig_id` in the detailed table; the different name makes a script that expects 0-based `contig_id` fail instead of reading 1-based coordinates) |
| `start`, `end` | 1-based, inclusive |
| `length_kb` | `(end − start + 1) / 1000`, 1 dp |
| `edges` | `att site` (`boundary` = `att`) or `estimated` (`model`) |
| `start_range`, `end_range` | the edge ranges, 1-based: `lo-hi`, or one number when the edge is exact (an att site or a contig end) |
| `support` | `seed_support` as `k/n` checkpoints (`1/1` with one checkpoint) |
| `note` | plain-language flags, `; `-separated: `att site: L bp repeat at both edges`; `runs off the contig start/end/both ends; the element may continue`; `P% N` (from 10% N); `on a short contig`; `start` / `end` / `start and end not refined: …` (an edge whose source is `coarse`: neither refinement nor composition moved it, because the fine score stays high past it) |

The terminal prints the top 10 the same way, with a two-line legend and a link to `docs/interpreting-results.md`.

## `regions.detailed.tsv`

One row per call, sorted by contig (in input order), then by start. 0-based, half-open. Before 2026-10-21 this table was `regions.tsv`, without `rank`.

| column | meaning |
|---|---|
| `call_id` | `prophound_1` … `prophound_N`; the same name in BED, GFF3 and FASTA |
| `rank` | as in `regions.tsv` |
| `genome_id` | the input file name without its extension (and without `.gz`/`.bz2`/`.xz`) |
| `contig_id` | the sequence ID: the first word of the FASTA header, or the versioned accession (`VERSION`) of a GenBank/EMBL record |
| `start`, `end` | the call, 0-based half-open. With `boundary` = `att`, from the first base of attL to the last base of attR (both copies of the att core inside); otherwise estimated edges (refinement, then composition; `start_source`/`end_source` say which) |
| `boundary` | `att`: both edges are an att site, a direct repeat with one copy at each edge left by integration, exact to the base. `model`: the edges are estimated, by the model's score and composition (docs/att-boundaries.md). The same values are `att site` and `estimated` in `regions.tsv`'s `edges` column and the terminal list; `att` and `model` are kept in this file and the GFF3, which scripts read |
| `start_lo`, `start_hi`, `end_lo`, `end_hi` | where each true edge lies. Equal to the edge for an att edge and for an edge on a contig end. For a model edge, a range relative to the call's length L (fitted as a 90% range; it held 0.782 of Default edges on the held-out second lockbox), clipped to the contig: from 0.143·L outside to 0.463·L inside the reported edge, in every mode (fitted on Default and Conservative: coverage 0.898 tuning, 0.891 and 0.883 on two check samples; Exploratory 0.909 tuning, 0.907 and 0.886 on the check samples, 0.804 on LAMBDA). Most of each range is inward, because model edges more often overshoot |
| `n_windows` | windows overlapping the model's call |
| `region_len_bp` | `end − start` |
| `mean_prob` | mean score of those windows (the fine windows when the call was refined), 4 dp. The rank's tie-break |
| `max_prob` | highest score among them, 4 dp; the rank's main key. Its median is about 1, but it separates calls on known elements from the rest better than `mean_prob` (an earlier version of this doc said the opposite, unmeasured) |
| `window_bp` | window length that made the call: `--window-len` (10,000 by default), or 5,000 / 3,000 on a contig shorter than that |
| `contig_edge` | `none`, `start`, `end` or `both`: an edge on a contig end. That edge is the end of the sequence, not a detected element boundary, and the element may continue on another contig |
| `n_frac` | fraction of the call that is N, 4 dp |
| `seed_support` | with several checkpoints: how many make an overlapping call of their own (absent with one checkpoint). Counted on the call as the vote saw it (`coarse_start`/`coarse_end` when refinement ran), so a split piece reports its parent call's support and never reads below its mode's vote |
| `refine_status`, `coarse_start`, `coarse_end`, `n_pieces` | from the refinement pass (on by default; absent with `--refine-window-len 0`). `refine_status` is `refined`, `split` (one coarse call became `n_pieces`), `no_fine_support` (kept unchanged) or `short_contig` (not refined). `coarse_*` is the call before refinement; the att search uses it |
| `model_start`, `model_end` | the model's edges after refinement, before att and composition |
| `start_source`, `end_source` | what placed each edge: `att`, `composition` (`--composition-edges`), `fine` (the refinement pass), `coarse` (nothing moved it off the 1 kb grid) or `contig_end` |
| `att_len`, `att_core` | the att core's length and sequence (empty for a model call). Authors differ on whether one core copy belongs to the element: with this core, `start + att_len` or `end − att_len` gives the other convention |

## `regions.filtered.tsv`

The same columns as `regions.detailed.tsv` (0-based), numbered `filtered_1` … `filtered_N`, plus:

| column | meaning |
|---|---|
| `filter` | why the call was set aside: `seed_vote` (fewer checkpoints than the mode's vote made it), `n_gap` (at least 50% N), or `seed_vote;n_gap` |

These calls are not in BED, GFF3 or FASTA. They are kept so a discovery user can recover them without a re-run. Refinement columns are empty for calls set aside by the vote, because only kept calls are refined. They are searched for att sites like kept calls.

## `regions.bed`, `regions.gff3`, `regions.fna`

- **BED6:** `contig_id`, `start`, `end`, `call_id`, score, `.`. The score is `round(1000 × mean_prob)`, for display; `rank` is the ranking.
- **GFF3:** type `prophage` (Sequence Ontology), source `prophound`, score `mean_prob`. Attributes: `ID` (= `call_id`), `note=candidate`, `boundary` (`att` or `model`, shown as `att site` / `estimated` in `regions.tsv`), `att_core` (att calls only), `mean_prob`, `max_prob`, `seed_support`, `n_frac`, `contig_edge`, `window_bp`.
- **FASTA:** `>call_id contig:start-end` (1-based, inclusive), 80 bases per line.

## `windows.tsv`

One row per scored window: full-length windows first, then short-contig windows.

| column | meaning |
|---|---|
| `genome_id`, `contig_id`, `start`, `end` | the window, 0-based half-open |
| `prob` | the ensemble score: the highest of the checkpoints' scores with `--ensemble max` |
| `pred_label` | 1 if `prob` reaches the score threshold: the mode's (`--threshold`) for full-length windows, the short-contig threshold for short ones |
| `prob_ckpt1` … | each checkpoint's own score (`--write-seed-scores`) |
| `region_pred_label` | 1 if the window overlaps a call |

Scores are written at full precision, because offline scorers rebuild calls from them. Read them with `pandas.read_csv(..., float_precision="round_trip")`; the default parser can be one unit in the last place off.

`fine_windows.tsv` has the same layout for the refinement windows (`genome_id`, `contig_id`, `start`, `end`, `prob`, `region_pred_label`).

## `run_metadata.json`

| key | meaning |
|---|---|
| `prophound_version`, `run_name` | package version, and the checkpoint's training run name |
| `bundled_model`, `config`, `checkpoint`, `checkpoint_sha256` | whether the bundled model ran; config and checkpoint paths, and each checkpoint's SHA-256 (the bundled ones are listed in `docs/lockbox.md`) |
| `input_fasta`, `input_sha256`, `input` | the input as given, the SHA-256 of that file, and a summary (contigs, bp, ambiguity codes) |
| `input_format` | detected `format` (`fasta`, `genbank`, `embl`, `gff3`), `compression`, `u_read_as_t` (RNA `U` bases read as `T`), and `contig_names` (GenBank/EMBL record name per contig ID, where they differ) |
| `mode`, `threshold`, `seed_threshold`, `extend_threshold`, `ensemble` | the thresholds actually used |
| `window_len`, `stride`, `trained_window_len`, `trained_stride`, `grid_offset`, `ambiguity_policy` | window settings |
| `region_builder`, `coverage_frac`, `fragment_merge`, `refine`, `boundary_refinement` | call-building settings (`refine` only with refinement); `boundary_refinement` names the edge stages that ran (`fine_pass_then_att_then_composition` by default in every mode) |
| `boundaries` | `att_search` and `composition_edges` (whether each ran), `calls_with_att`, `edges_by_source`, and the `edge_range` fractions and coverage |
| `short_contigs` | contigs under the window length: counts, and the thresholds used to scan them |
| `seed_support` | the vote: checkpoints, minimum support, calls set aside, calls by support |
| `edge_snap`, `n_gap_filter` | edge snapping and the N filter: settings and counts |
| `n_windows`, `n_regions`, `outputs` | counts, and every file written |
| `device`, `started_utc`, `elapsed_s`, `command`, `platform` | where and when it ran, the command as typed, and Python, OS and torch versions |

Every path shows the home directory as `~`, so a metadata file can be shared without revealing a user name.
