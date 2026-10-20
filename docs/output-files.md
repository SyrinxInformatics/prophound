# Output files

What `prophound predict` writes, column by column. This is the current output contract (since the UI pass, step 2, 2026-10-20). The frozen `notes/V1_OUTPUT_FILES.md` describes an older format. Update this file whenever an output changes.

All coordinates in the TSV files and BED are **0-based, half-open** (`start` included, `end` excluded), as in BED. GFF3 and the FASTA headers are 1-based, inclusive, as those formats require.

## Where outputs go

- Everything goes into `--outdir`, which is created if missing.
- **An outdir that already holds prophound outputs is refused.** The error names the next free directory (`out_1`, `out_2`, …). `--force` replaces the outputs, after removing every file prophound can write there, so no file from an earlier run survives. Files that aren't prophound's are left alone.
- **Writes are atomic.** Each file is written under a temporary name in the outdir, and all of them are renamed once every file is complete. A failed or interrupted run leaves no partial set of outputs. Earlier outputs stay intact, even with `--force`, unless the run reaches the final rename.

## Files

| file | written | contents |
|---|---|---|
| `regions.tsv` | always | the calls: candidate prophage regions to follow up |
| `regions.filtered.tsv` | always | candidates set aside by the seed vote or the N filter, with the reason |
| `regions.bed` | always | the calls as BED6 |
| `regions.gff3` | always | the calls as GFF3 |
| `windows.tsv` | always | every scored window and its score |
| `run_metadata.json` | always | settings, input and checkpoint hashes, counts, timing, platform |
| `regions.fna` | `--write-region-fasta` | the calls' sequences |
| `fine_windows.tsv` | `--refine-window-len` > 0 | the refinement pass's windows |

## Batches

A single `--input` file writes the layout above. A **batch** (`--input` given more than once, a directory, or `--input-list`) writes one directory per genome, `outdir/<genome_id>/`, each holding that genome's files exactly as a single run would write them. The top of the outdir holds:

| file | contents |
|---|---|
| `regions.tsv`, `regions.filtered.tsv` | every genome's calls (or set-aside candidates), in input order. `genome_id` tells them apart; `call_id` is numbered per genome |
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

One row per call, sorted by contig (in input order), then by start.

| column | meaning |
|---|---|
| `call_id` | `prophound_1` … `prophound_N`; the same name in BED, GFF3 and FASTA |
| `genome_id` | the input file name without its extension (and without `.gz`/`.bz2`/`.xz`) |
| `contig_id` | the sequence ID: the first word of the FASTA header, or the versioned accession (`VERSION`) of a GenBank/EMBL record |
| `start`, `end` | the call, 0-based half-open |
| `n_windows` | windows overlapping the call |
| `region_len_bp` | `end − start` |
| `mean_prob` | mean score of those windows, 4 dp. **Use this to rank calls**: it separates calls that match known elements from those that don't better than `max_prob` |
| `max_prob` | highest score among them, 4 dp; about 1 for almost every call |
| `window_bp` | window length that made the call: `--window-len` (10,000 by default), or 5,000 / 3,000 on a contig shorter than that |
| `contig_edge` | `none`, `start`, `end` or `both`: an edge on a contig end. That edge is the end of the sequence, not a detected element boundary, and the element may continue on another contig |
| `n_frac` | fraction of the call that is N, 4 dp |
| `seed_support` | with several checkpoints: how many make an overlapping call of their own (absent with one checkpoint) |
| `refine_status`, `coarse_start`, `coarse_end`, `n_pieces` | only with `--refine-window-len`. `refine_status` is `refined`, `split` (one coarse call became `n_pieces`), `no_fine_support` (kept unchanged) or `short_contig` (not refined). `coarse_*` is the call before refinement |

## `regions.filtered.tsv`

The same columns as `regions.tsv`, numbered `filtered_1` … `filtered_N`, plus:

| column | meaning |
|---|---|
| `filter` | why the call was set aside: `seed_vote` (fewer checkpoints than the mode's vote made it), `n_gap` (at least 50% N), or `seed_vote;n_gap` |

These calls are not in BED, GFF3 or FASTA. They are kept so a discovery user can recover them without a re-run. Refinement columns are empty for calls set aside by the vote, because only kept calls are refined.

## `regions.bed`, `regions.gff3`, `regions.fna`

- **BED6:** `contig_id`, `start`, `end`, `call_id`, score, `.`. The score is `round(1000 × mean_prob)`, so it ranks calls.
- **GFF3:** type `prophage` (Sequence Ontology), source `prophound`, score `mean_prob`. Attributes: `ID` (= `call_id`), `note=candidate`, `mean_prob`, `max_prob`, `seed_support`, `n_frac`, `contig_edge`, `window_bp`.
- **FASTA:** `>call_id contig:start-end` (1-based, inclusive), 80 bases per line.

## `windows.tsv`

One row per scored window: full-length windows first, then short-contig windows.

| column | meaning |
|---|---|
| `genome_id`, `contig_id`, `start`, `end` | the window, 0-based half-open |
| `prob` | the ensemble score: the highest of the checkpoints' scores with `--ensemble max` |
| `pred_label` | 1 if `prob` reaches the score threshold: the mode's (`--threshold`) for full-length windows, the short-contig threshold for short ones |
| `prob_ckpt1` … | each checkpoint's own score (`--write-seed-scores`) |
| `region_pred_label` | 1 if the window overlaps a call in `regions.tsv` |

Scores are written at full precision, because offline scorers rebuild calls from them. Read them with `pandas.read_csv(..., float_precision="round_trip")`; the default parser can be one unit in the last place off.

`fine_windows.tsv` has the same layout for the refinement windows (`genome_id`, `contig_id`, `start`, `end`, `prob`, `region_pred_label`).

## `run_metadata.json`

| key | meaning |
|---|---|
| `prophound_version`, `run_name` | package version, and the checkpoint's training run name |
| `config`, `checkpoint`, `checkpoint_sha256` | config and checkpoint paths, and each checkpoint's SHA-256 |
| `input_fasta`, `input_sha256`, `input` | the input as given, the SHA-256 of that file, and a summary (contigs, bp, ambiguity codes) |
| `input_format` | detected `format` (`fasta`, `genbank`, `embl`, `gff3`), `compression`, `u_read_as_t` (RNA `U` bases read as `T`), and `contig_names` (GenBank/EMBL record name per contig ID, where they differ) |
| `mode`, `threshold`, `seed_threshold`, `extend_threshold`, `ensemble` | the thresholds actually used |
| `window_len`, `stride`, `trained_window_len`, `trained_stride`, `grid_offset`, `ambiguity_policy` | window settings |
| `region_builder`, `coverage_frac`, `fragment_merge`, `refine` | call-building settings (`refine` only with refinement) |
| `short_contigs` | contigs under the window length: counts, and the thresholds used to scan them |
| `seed_support` | the vote: checkpoints, minimum support, calls set aside, calls by support |
| `edge_snap`, `n_gap_filter` | edge snapping and the N filter: settings and counts |
| `n_windows`, `n_regions`, `outputs` | counts, and every file written |
| `device`, `started_utc`, `elapsed_s`, `command`, `platform` | where and when it ran, the command as typed, and Python, OS and torch versions |

Every path shows the home directory as `~`, so a metadata file can be shared without revealing a user name.
