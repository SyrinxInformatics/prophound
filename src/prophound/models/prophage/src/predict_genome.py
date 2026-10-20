# SPDX-License-Identifier: Apache-2.0
#
# Copyright 2025 Jacob Alford <jalford0000@gmail.com>
#                and Benjamin Alford <ssolar5.dev@gmail.com>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""`prophound predict`: score a genome's windows with the checkpoints and write the calls.

The stages are in `run_predict`. Call building is in `calls.py`, shared with the offline
scorers; options and defaults are in `options.py`.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import re
import shlex
import shutil
import sys
import tempfile
import time
import traceback
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd
import torch
import yaml
from pyfaidx import Fasta
from torch.utils.data import DataLoader, Dataset
from tqdm import tqdm

from .calls import (  # noqa: F401 -- re-exported: scripts and tests import these from here
    REGION_COLUMNS,
    CallRules,
    PassCalls,
    apply_seed_vote,
    build_pass,
    build_region_table_coverage,
    concat_nonempty,
    finalize,
    label_windows,
    seed_support,
    snap_contig_edges,
)
from .config import DataConfig, load_cfg  # noqa: F401 -- DataConfig is used by scripts/benchmark.py
from .devices import resolve_device
from .encoding import collate_batch, encode_window, onehot_encode, stable_sigmoid  # noqa: F401 -- re-exported
from .inputs import PreparedInput, prepare_input
from .io_checks import package_version, sha256, validate_fasta, write_bed, write_gff3, write_region_fasta
from .model_rc_cnn import build_model
from .multiscale import plan_fine_windows, plan_neighborhoods, refine_regions
from .options import (  # noqa: F401 -- re-exported: scripts import the defaults and tables from here
    AMBIGUITY_POLICIES,
    DEFAULT_AMBIGUITY_POLICY,
    DEFAULT_BRIDGE_PROB,
    DEFAULT_COVERAGE_FRAC,
    DEFAULT_EDGE_SNAP_FRAC,
    DEFAULT_MAX_GAP_BP,
    MAX_CALL_N_FRAC,
    MODE_MIN_SEED_SUPPORT,
    MODE_THRESHOLDS,
    PREDICT_EPILOG,
    SHORT_CONTIG_SCALES,
    SHORT_CONTIG_STRIDE,
    SHORT_CONTIG_THRESHOLDS,
    PredictError,
    PredictOptions,
    add_predict_arguments,
)

# every file `predict` can write; --force removes all of them before writing
OUTPUT_FILES = ("windows.tsv", "regions.tsv", "regions.filtered.tsv", "run_metadata.json", "regions.bed",
                "regions.gff3", "regions.fna", "fine_windows.tsv")
WINDOW_COLUMNS = ["genome_id", "contig_id", "start", "end", "fasta_path", "prob", "pred_label", "region_pred_label"]


class Log:
    """stderr output at the run's verbosity: warnings always, progress lines unless
    --quiet, `[info]` stage details only with --verbose."""

    def __init__(self, verbosity: int = 0):
        self.verbosity = verbosity

    def warn(self, msg: str) -> None:
        print(f"[warn] {msg}", file=sys.stderr)

    def info(self, msg: str) -> None:
        if self.verbosity >= 0:
            print(msg, file=sys.stderr)

    def detail(self, msg: str) -> None:
        if self.verbosity > 0:
            print(f"[info] {msg}", file=sys.stderr)


class Progress:
    """One progress bar for every scoring pass of a run (hidden with --quiet or off a terminal)."""

    def __init__(self, total: int, quiet: bool):
        self.bar = tqdm(total=total, unit="window", unit_scale=True, desc="scoring", leave=False,
                        file=sys.stderr, dynamic_ncols=True, disable=True if quiet else None)

    def add(self, n: int) -> None:
        """Grow the total once a pass's size is known (the refinement pass)."""
        self.bar.total += n
        self.bar.refresh()

    def update(self, n: int) -> None:
        self.bar.update(n)

    def close(self) -> None:
        self.bar.close()


class GenomeWindowDataset(Dataset):
    """Window dataset for direct inference on a FASTA file."""

    def __init__(self, windows_df: pd.DataFrame, window_len: int, ambiguity_policy: str = DEFAULT_AMBIGUITY_POLICY):
        self.windows_df = windows_df.reset_index(drop=True)
        self.window_len = window_len
        self.ambiguity_policy = ambiguity_policy
        fasta_paths = self.windows_df["fasta_path"].drop_duplicates().tolist()
        if len(fasta_paths) != 1:
            raise ValueError("GenomeWindowDataset expects exactly one FASTA path.")
        self.fasta_path = fasta_paths[0]
        self.fasta = Fasta(self.fasta_path, as_raw=True, sequence_always_upper=True)

    def __len__(self) -> int:
        return len(self.windows_df)

    def __getitem__(self, index: int) -> Tuple[torch.Tensor, torch.Tensor]:
        row = self.windows_df.iloc[index]
        contig_id = row["contig_id"]
        start = int(row["start"])
        end = int(row["end"])
        sequence = self.fasta[contig_id][start:end]
        if len(sequence) < self.window_len:
            # build_windows emits only full-length windows; padding with N could
            # silently zero the score (docs/scientific-validity.md "1a"), so fail instead
            raise RuntimeError(
                f"short window {contig_id}:{start}-{end} ({len(sequence)} < {self.window_len} bp); "
                "padding with N would silently zero its score"
            )
        elif len(sequence) > self.window_len:
            sequence = sequence[: self.window_len]
        return encode_window(sequence, self.ambiguity_policy), torch.tensor(0.0, dtype=torch.float32)

    def close(self) -> None:
        try:
            self.fasta.close()
        except Exception:
            pass

    def __del__(self):
        self.close()


@torch.no_grad()
def predict_probs(
    model: torch.nn.Module,
    loader: DataLoader,
    device: torch.device,
    on_batch=None,
) -> np.ndarray:
    """Score every batch of `loader`; `on_batch(n)` is called after each batch of n windows."""
    model.eval()
    all_logits: List[np.ndarray] = []
    for batch_x, _ in loader:
        batch_x = batch_x.to(device, non_blocking=True)
        logits = model(batch_x).view(-1)
        all_logits.append(logits.detach().cpu().numpy())
        if on_batch is not None:
            on_batch(len(batch_x))
    if not all_logits:
        return np.array([], dtype=np.float32)
    return stable_sigmoid(np.concatenate(all_logits))


def iter_fasta_lengths(fasta_path: str) -> Iterable[Tuple[str, int]]:
    fasta = Fasta(fasta_path, as_raw=True, sequence_always_upper=True)
    try:
        for contig_id in fasta.keys():
            yield contig_id, len(fasta[contig_id])
    finally:
        fasta.close()


def build_windows(
    fasta_path: str,
    genome_id: str,
    window_len: int,
    stride: int,
    offset: int = 0,
    contig_lengths: Dict[str, int] | None = None,
) -> pd.DataFrame:
    """Full-length windows over every contig at least `window_len` long."""
    lengths = contig_lengths if contig_lengths is not None else dict(iter_fasta_lengths(fasta_path))
    rows: List[Dict[str, object]] = []
    for contig_id, contig_length in lengths.items():
        if contig_length < window_len:
            continue
        for start in range(offset, contig_length - window_len + 1, stride):
            rows.append({"genome_id": genome_id, "contig_id": contig_id, "start": int(start),
                         "end": int(start + window_len), "fasta_path": fasta_path})
    return pd.DataFrame(rows)


def build_short_contig_windows(
    fasta_path: str,
    genome_id: str,
    window_len: int,
    scales: Tuple[int, ...] = SHORT_CONTIG_SCALES,
    stride: int = SHORT_CONTIG_STRIDE,
    contig_lengths: Dict[str, int] | None = None,
) -> Tuple[pd.DataFrame, Dict[str, int]]:
    """Windows for contigs shorter than `window_len`, plus a report of them all.

    `build_windows` skips these contigs. Each one at least as long as the
    smallest scale is tiled with the largest scale that fits (column `scale`);
    shorter ones stay unscanned. The report always covers every short contig, so
    the caller can warn even when scanning is off.
    """
    rows: List[Dict[str, object]] = []
    report = {"n_short": 0, "bp_short": 0, "n_scannable": 0, "n_too_short": 0, "bp_too_short": 0}
    lengths = contig_lengths if contig_lengths is not None else dict(iter_fasta_lengths(fasta_path))
    for contig_id, length in lengths.items():
        if length >= window_len:
            continue
        report["n_short"] += 1
        report["bp_short"] += int(length)
        fitting = [sc for sc in scales if sc <= length and sc < window_len]
        if not fitting:
            report["n_too_short"] += 1
            report["bp_too_short"] += int(length)
            continue
        report["n_scannable"] += 1
        scale = max(fitting)
        for start in range(0, length - scale + 1, stride):
            rows.append({"genome_id": genome_id, "contig_id": contig_id, "start": int(start),
                         "end": int(start + scale), "fasta_path": fasta_path, "scale": int(scale)})
    cols = ["genome_id", "contig_id", "start", "end", "fasta_path", "scale"]
    return pd.DataFrame(rows, columns=cols), report


def call_n_fraction(regions_df: pd.DataFrame, fasta_path: str) -> np.ndarray:
    """Fraction of each call's bases that are N."""
    if regions_df.empty:
        return np.zeros(0)
    fa = Fasta(fasta_path, as_raw=True, sequence_always_upper=True)
    try:
        return np.array([fa[r.contig_id][int(r.start):int(r.end)].count("N") / max(1, int(r.end) - int(r.start))
                         for r in regions_df.itertuples(index=False)])
    finally:
        fa.close()


def load_checkpoints(model: torch.nn.Module, paths: Tuple[str, ...], device: torch.device,
                     config_path: str) -> List[dict]:
    """Read each checkpoint's weights once, checking that they fit the model the config describes."""
    states = []
    for path in paths:
        try:
            state = torch.load(path, map_location=device, weights_only=False)["model"]
        except (KeyError, TypeError, IndexError):
            raise PredictError(f"{path} is not a prophound checkpoint (it has no model weights)") from None
        except Exception as e:  # unreadable, truncated, not a torch file
            raise PredictError(f"cannot read checkpoint {path}: {e}") from None
        try:
            model.load_state_dict(state)
        except RuntimeError as e:
            # torch's first line is generic; the next names the first mismatched weight
            detail = next((ln.strip() for ln in str(e).splitlines()[1:] if ln.strip()), str(e))
            raise PredictError(f"checkpoint {path} does not fit the model described by {config_path} "
                               f"({detail.rstrip('.')})") from None
        states.append(state)
    return states


class Scorer:
    """The model with every checkpoint's weights, and the run's progress bar."""

    def __init__(self, model: torch.nn.Module, states: List[dict], device: torch.device, opts: PredictOptions,
                 progress: Progress):
        self.model, self.states, self.device, self.opts, self.progress = model, states, device, opts, progress

    def score(self, windows_df: pd.DataFrame, window_len: int) -> Tuple[np.ndarray, List[np.ndarray]]:
        """Score `windows_df` with every checkpoint: the ensemble track and each checkpoint's own."""
        opts, device = self.opts, self.device
        dataset = GenomeWindowDataset(windows_df, window_len=window_len, ambiguity_policy=opts.ambiguity_policy)
        loader = DataLoader(
            dataset,
            batch_size=opts.batch_size,
            shuffle=False,
            num_workers=opts.num_workers,
            collate_fn=collate_batch,
            pin_memory=(device.type == "cuda"),
        )
        per_checkpoint = []
        for state in self.states:
            self.model.load_state_dict(state)
            probs = predict_probs(self.model, loader, device, on_batch=self.progress.update)
            per_checkpoint.append(np.asarray(probs, dtype=np.float64))
        dataset.close()

        if len(per_checkpoint) == 1:
            return per_checkpoint[0], per_checkpoint
        stacked = np.stack(per_checkpoint)
        # max, not mean: averaging pulls borderline seeds under the seed threshold and
        # costs recall (docs/cli-and-packaging.md "Defaults changed + seed ensembling added")
        ensemble = stacked.max(axis=0) if opts.ensemble == "max" else stacked.mean(axis=0)
        return ensemble, per_checkpoint


def scored(windows_df: pd.DataFrame, ensemble: np.ndarray, per_checkpoint: List[np.ndarray],
           threshold: float, opts: PredictOptions) -> pd.DataFrame:
    """Windows with their scores: prob, pred_label, and prob_ckpt<k> with --write-seed-scores."""
    out = windows_df.copy()
    out["prob"] = ensemble
    out["pred_label"] = (out["prob"] >= threshold).astype(int)
    if opts.write_seed_scores:
        for k, p in enumerate(per_checkpoint):
            out[f"prob_ckpt{k + 1}"] = p
    return out


def call_coarse(scorer: Scorer, windows_df: pd.DataFrame, opts: PredictOptions, rules: CallRules) -> PassCalls:
    """The full-length window pass: score, build calls with fragment merging, vote."""
    if windows_df.empty:   # every contig is short
        empty = pd.DataFrame(columns=REGION_COLUMNS)
        return PassCalls(pd.DataFrame(columns=WINDOW_COLUMNS), empty, empty.assign(filter=pd.Series(dtype=str)), [])
    ensemble, per_checkpoint = scorer.score(windows_df, opts.window_len)
    return build_pass(scored(windows_df, ensemble, per_checkpoint, opts.threshold, opts), per_checkpoint,
                      opts.seed_threshold, opts.extend_threshold, rules,
                      max_gap_bp=opts.max_gap_bp, bridge_prob=opts.bridge_prob)


def refine(scorer: Scorer, calls: pd.DataFrame, contig_lengths: Dict[str, int], genome_id: str, fasta: str,
           opts: PredictOptions) -> Tuple[pd.DataFrame, pd.DataFrame | None, dict | None]:
    """Fine-scale pass over the coarse calls (off unless --refine-window-len). Never drops a call.

    Returns the refined calls, the scored fine windows, and the refine metadata.
    """
    if not opts.refine_window_len:
        return calls, None, None
    fine_windows_df, counts = None, {}
    if not calls.empty:
        neighborhoods = plan_neighborhoods(calls, contig_lengths, opts.refine_margin)
        fine_windows_df = plan_fine_windows(
            neighborhoods, contig_lengths, genome_id, fasta,
            window_len=opts.refine_window_len, stride=opts.refine_stride,
        )
        if not fine_windows_df.empty:
            fine_windows_df = fine_windows_df.copy()
            scorer.progress.add(len(fine_windows_df) * len(opts.checkpoints))
            fine_windows_df["prob"] = scorer.score(fine_windows_df, opts.refine_window_len)[0]
            # one threshold for seed and extend at the fine scale; no merging
            # here -- merging belongs to the coarse track (merge coarse, split fine)
            fine_windows_df, fine_regions_df = build_region_table_coverage(
                fine_windows_df,
                seed_threshold=opts.refine_threshold,
                extend_threshold=opts.refine_threshold,
                min_coverage_frac=opts.coverage_frac,
                max_gap_bp=0,
            )
        else:
            fine_regions_df = pd.DataFrame(columns=["contig_id", "start", "end"])
        ramp = opts.edge_alpha_ramp
        calls = refine_regions(
            calls,
            fine_windows_df,
            fine_regions_df,
            edge_alpha=opts.edge_alpha,
            split_min_gap_bp=opts.split_min_gap_bp,
            split_max_prob=opts.split_max_prob,
            kind_suffix=f"_refine{opts.refine_window_len}_a{opts.edge_alpha:g}"
            + (f"_ramp{ramp[0]}-{ramp[1]}" if ramp else ""),
            alpha_ramp_bp=ramp,
        )
        counts = {k: int(v) for k, v in calls["refine_status"].value_counts().items()}
    meta = {
        "window_len": opts.refine_window_len, "stride": opts.refine_stride,
        "threshold": opts.refine_threshold, "margin_bp": opts.refine_margin,
        "edge_alpha": opts.edge_alpha,
        "edge_alpha_ramp_bp": list(opts.edge_alpha_ramp) if opts.edge_alpha_ramp else None,
        "split_min_gap_bp": opts.split_min_gap_bp, "split_max_prob": opts.split_max_prob,
        "n_fine_windows": int(0 if fine_windows_df is None else len(fine_windows_df)),
        "status_counts": counts,
    }
    return calls, fine_windows_df, meta


def call_short_contigs(scorer: Scorer, short_windows_df: pd.DataFrame, opts: PredictOptions,
                       rules: CallRules, log: Log) -> Tuple[List[PassCalls], dict]:
    """Contigs shorter than --window-len, one pass per scale (5 or 3 kb), no merging.

    Returns the passes and the metadata additions (thresholds used, etc.); no passes
    when scanning is off or nothing is scannable.
    """
    if opts.short_contigs != "scan" or short_windows_df.empty:
        return [], {}
    if opts.short_contig_threshold is not None:
        log.warn(f"--short-contig-threshold {float(opts.short_contig_threshold):g} overrides the calibrated "
             "short-contig thresholds.")
    passes, used = [], {}
    for scale, sw in short_windows_df.groupby("scale"):
        t = opts.short_threshold(int(scale))
        used[str(int(scale))] = {"seed": t, "extend": t}
        sw = sw.drop(columns="scale").reset_index(drop=True)
        ensemble, per_checkpoint = scorer.score(sw, int(scale))
        # per-contig windows share one length, so the coverage builder's
        # per-contig window_len/stride inference holds
        passes.append(build_pass(scored(sw, ensemble, per_checkpoint, t, opts), per_checkpoint, t, t, rules,
                                 max_gap_bp=0, kind_suffix=f"_short{int(scale)}"))
    meta = {"thresholds": used, "scales": list(SHORT_CONTIG_SCALES), "stride": SHORT_CONTIG_STRIDE,
            "n_regions": int(sum(len(p.calls) for p in passes)), "calibrated": opts.short_contig_threshold is None}
    return passes, meta


def combine_passes(coarse: PassCalls, coarse_calls: pd.DataFrame,
                   short: List[PassCalls]) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Windows, calls and filtered candidates of every pass, coarse first."""
    windows_df, calls, filtered = coarse.windows, coarse_calls, coarse.filtered
    if short:
        short_calls = concat_nonempty([p.calls for p in short])
        if "refine_status" in calls.columns:
            short_calls = short_calls.assign(refine_status="short_contig", coarse_start=short_calls["start"],
                                             coarse_end=short_calls["end"], n_pieces=1)
        # concatenating onto an empty frame would change the column dtypes
        calls = concat_nonempty([calls, short_calls]) if not calls.empty else short_calls
        short_windows = [p.windows for p in short]
        windows_df = pd.concat([windows_df, *short_windows], ignore_index=True) \
            if not windows_df.empty else pd.concat(short_windows, ignore_index=True)
        filtered = concat_nonempty([filtered, *(p.filtered for p in short)])
    return windows_df, calls, filtered


def call_rules(opts: PredictOptions) -> CallRules:
    return CallRules(
        seed_threshold=opts.seed_threshold, extend_threshold=opts.extend_threshold,
        coverage_frac=opts.coverage_frac, max_gap_bp=opts.max_gap_bp, bridge_prob=opts.bridge_prob,
        min_seed_support=opts.min_seed_support,
        short_thresholds={sc: opts.short_threshold(sc) for sc in SHORT_CONTIG_SCALES},
        snap_bp=opts.snap_bp, max_call_n_frac=MAX_CALL_N_FRAC,
    )


def _load_config(path: str):
    try:
        return load_cfg(path)
    except FileNotFoundError:
        raise PredictError(f"config not found: {path}") from None
    except ValueError as e:   # load_cfg's messages start with the path
        raise PredictError(f"config {e}") from None
    except (OSError, yaml.YAMLError) as e:
        raise PredictError(f"cannot read config {path}: {e}") from None


def _check_paths(opts: PredictOptions) -> None:
    """Fail before any work on a missing checkpoint or an outdir that cannot be written."""
    missing = [c for c in opts.checkpoints if not Path(c).is_file()]
    if missing:
        raise PredictError(f"checkpoint not found: {', '.join(missing)}")
    out = Path(opts.outdir)
    if out.exists() and not out.is_dir():
        raise PredictError(f"--outdir {out} exists and is not a directory")
    existing = next(p for p in (out, *out.parents) if p.exists())
    if not os.access(existing, os.W_OK | os.X_OK):
        raise PredictError(f"cannot write to --outdir {out}: no write permission on {existing}")
    held = _held_outputs(out)
    if held and not opts.force:
        raise PredictError(f"--outdir {out} already holds prophound outputs ({', '.join(held)}); "
                           f"use a new directory such as {_next_free(out)}, or --force to replace them")


def _held_outputs(d: Path) -> List[str]:
    return [name for name in OUTPUT_FILES if (d / name).exists()] if d.is_dir() else []


def _next_free(out: Path) -> Path:
    """The first of out_1, out_2, ... that holds no outputs."""
    k = 1
    while _held_outputs(out.with_name(f"{out.name}_{k}")):
        k += 1
    return out.with_name(f"{out.name}_{k}")


def _nothing_to_scan(opts: PredictOptions, contig_lengths: Dict[str, int]) -> PredictError:
    longest = max(contig_lengths.values())
    if opts.short_contigs == "skip":
        return PredictError(f"nothing to scan: every contig is shorter than --window-len ({opts.window_len:,} bp), "
                            f"and --short-contigs skip leaves those unscanned (the longest is {longest:,} bp)")
    return PredictError(f"nothing to scan: every contig is shorter than {min(SHORT_CONTIG_SCALES):,} bp, the "
                        f"smallest window (the longest is {longest:,} bp)")


def run_predict(args: argparse.Namespace) -> None:
    """Run inference end to end from parsed `predict` arguments."""
    started = time.perf_counter()
    started_utc = datetime.now(timezone.utc).isoformat(timespec="seconds")
    cfg = _load_config(args.config)
    device = resolve_device(args.device)
    opts = PredictOptions.from_args(args, cfg.data.window_len, cfg.data.stride)
    log = Log(opts.verbosity)
    _check_paths(opts)
    if opts.window_len != opts.trained_window_len:
        # the model accepts any length, but its scores on negatives rise as windows
        # shrink, so the thresholds don't transfer (docs/adaptive-multiscale-detection.md)
        log.warn(f"--window-len {opts.window_len:,} differs from the trained {opts.trained_window_len:,}; the "
                 f"mode thresholds were calibrated at {opts.trained_window_len:,} bp and are not validated "
                 "at this length.")

    with input_workspace(Path(opts.outdir)) as work:
        prepared = prepare_input(opts.input_fasta, work)
        predict_prepared(opts, cfg, device, log, prepared, started, started_utc)


@contextmanager
def input_workspace(outdir: Path):
    """A hidden work directory in the outdir for the normalised input and its index,
    removed afterwards (with the outdir too, if this run created it and wrote nothing)."""
    created = not outdir.exists()
    outdir.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix=".prophound-input-", dir=outdir))
    try:
        yield work
    finally:
        shutil.rmtree(work, ignore_errors=True)
        if created and not any(outdir.iterdir()):
            outdir.rmdir()


def predict_prepared(opts: PredictOptions, cfg, device: torch.device, log: Log, prepared: PreparedInput,
                     started: float, started_utc: str) -> None:
    """Everything after the input is ready: windows, scoring, calls and outputs."""
    fasta = prepared.fasta
    if prepared.u_as_t:
        log.warn(f"{_plural(prepared.u_as_t, 'U base')} read as T (RNA sequence?).")
    input_summary = validate_fasta(fasta, label=opts.input_fasta)
    contig_lengths = dict(iter_fasta_lengths(fasta))
    genome_id = prepared.genome_id
    windows_df = build_windows(fasta, genome_id, opts.window_len, opts.stride,
                               offset=opts.grid_offset, contig_lengths=contig_lengths)
    short_windows_df, short_report = build_short_contig_windows(
        fasta, genome_id, window_len=opts.window_len, contig_lengths=contig_lengths)
    if short_report["n_short"]:
        # never silent (docs/cli-and-packaging.md "Input robustness")
        verb = "scanned with shorter windows" if opts.short_contigs == "scan" else "NOT scanned (--short-contigs skip)"
        log.warn(f"{_plural(short_report['n_short'], 'contig')} shorter than {opts.window_len:,} bp "
                 f"({short_report['bp_short']:,} bp): {short_report['n_scannable']:,} {verb}; "
                 f"{short_report['n_too_short']:,} under {min(SHORT_CONTIG_SCALES):,} bp never scanned "
                 f"({short_report['bp_too_short']:,} bp).")
    scan_short = opts.short_contigs == "scan" and not short_windows_df.empty
    if windows_df.empty and not scan_short:
        raise _nothing_to_scan(opts, contig_lengths)

    model = build_model(cfg, device)
    states = load_checkpoints(model, opts.checkpoints, device, opts.config)
    if opts.requested_min_seed_support > 1 and opts.min_seed_support < opts.requested_min_seed_support:
        # never silent: the mode thresholds were set with the vote over 3 checkpoints
        k = opts.min_seed_support
        applied = f"voting at {k}" if k > 1 else "no vote applied"
        log.warn(f"seed vote of {opts.requested_min_seed_support} requested but only "
                 f"{_plural(opts.n_checkpoints, 'checkpoint')}: {applied}. The '{opts.mode}' "
                 "thresholds were calibrated with the vote over 3 checkpoints; expect more calls.")
    rules = call_rules(opts)
    n_short = len(short_windows_df) if scan_short else 0
    log.detail(f"{len(windows_df):,} windows of {opts.window_len:,} bp and {n_short:,} short-contig windows, "
               f"{_plural(opts.n_checkpoints, 'checkpoint')}, on {device}")

    progress = Progress((len(windows_df) + n_short) * opts.n_checkpoints, quiet=opts.verbosity < 0)
    scorer = Scorer(model, states, device, opts, progress)
    try:
        t = time.perf_counter()
        coarse = call_coarse(scorer, windows_df, opts, rules)
        log.detail(f"full-length windows: {_plural(len(coarse.calls), 'call')} in {_duration(time.perf_counter() - t)}")
        t = time.perf_counter()
        coarse_calls, fine_windows_df, refine_meta = refine(scorer, coarse.calls, contig_lengths, genome_id,
                                                            fasta, opts)
        if refine_meta:
            log.detail(f"refinement: {_plural(len(coarse_calls), 'call')} in {_duration(time.perf_counter() - t)}")
        t = time.perf_counter()
        short, short_scan_meta = call_short_contigs(scorer, short_windows_df, opts, rules, log)
        if short:
            log.detail(f"short contigs: {_plural(sum(len(p.calls) for p in short), 'call')} "
                       f"in {_duration(time.perf_counter() - t)}")
    finally:
        progress.close()
    windows_df, calls, filtered = combine_passes(coarse, coarse_calls, short)

    per_seed = [[p.per_seed[k] for p in (coarse, *short) if p.per_seed] for k in range(opts.n_checkpoints)] \
        if opts.n_checkpoints > 1 else []
    final = finalize(calls, filtered, per_seed, contig_lengths,
                     lambda df: call_n_fraction(df, fasta), rules, opts.edge_snap_frac,
                     requested_min_seed_support=opts.requested_min_seed_support)
    gap_meta = final.meta["n_gap_filter"]
    if gap_meta["calls_dropped"]:
        log.warn(f"{_plural(gap_meta['calls_dropped'], 'call')} set aside in regions.filtered.tsv: at least "
                 f"{MAX_CALL_N_FRAC:.0%} N (imputed gap sequence).")
    windows_df = label_windows(windows_df, final.calls)

    short_meta = {"mode": opts.short_contigs, **short_report, **short_scan_meta}
    contig_order = {c: i for i, c in enumerate(contig_lengths)}
    calls_out = call_table(final.calls, contig_order, "prophound_", opts.window_len)
    filtered_out = call_table(final.filtered, contig_order, "filtered_", opts.window_len)
    filtered_out = filtered_out.reindex(columns=[*calls_out.columns, "filter"])
    run_info = {"started_utc": started_utc, "elapsed_s": round(time.perf_counter() - started, 1)}
    written = write_outputs(opts, cfg.run_name, str(device), prepared, input_summary, windows_df, calls_out,
                            filtered_out, fine_windows_df, short_meta, final.meta, refine_meta, run_info)
    log.info(f"[done] wrote {', '.join(written)} to {opts.outdir}")
    print_run_summary(opts, prepared, input_summary, calls_out, short_meta, final.meta.get("seed_support"), gap_meta,
                      str(device), time.perf_counter() - started)


_SHORT_SCALE = re.compile(r"_short(\d+)")


def call_table(calls: pd.DataFrame, contig_order: Dict[str, int], prefix: str, window_len: int) -> pd.DataFrame:
    """Calls as written: sorted by contig (input order) then start, numbered `<prefix>1..N`
    in a leading call_id column, probabilities to 4 dp, and the internal `boundary_kind`
    replaced by `window_bp`, the window length that made the call."""
    out = calls.copy()
    if "boundary_kind" in out.columns:
        kinds = out["boundary_kind"].astype(str)
        window_bp = [int(m.group(1)) if (m := _SHORT_SCALE.search(k)) else window_len for k in kinds]
        out.insert(out.columns.get_loc("boundary_kind"), "window_bp", window_bp)
        out = out.drop(columns="boundary_kind")
    out = (out.assign(_order=out["contig_id"].map(contig_order))
              .sort_values(["_order", "start", "end"], kind="stable")
              .drop(columns="_order").reset_index(drop=True))
    out.insert(0, "call_id", [f"{prefix}{i}" for i in range(1, len(out) + 1)])
    for col in ("mean_prob", "max_prob"):
        if col in out.columns:
            out[col] = out[col].astype(float).round(4)
    return out


def _shown(path: str) -> str:
    """A path for run_metadata.json, with the home directory written as ~."""
    p, home = str(path), str(Path.home())
    return "~" + p[len(home):] if p == home or p.startswith(home + os.sep) else p


def _command_line(argv: List[str]) -> str:
    """The command as a user would type it, home directory as ~ (left unquoted so it expands)."""
    def arg(a: str) -> str:
        shown = _shown(a)
        return "~/" + shlex.quote(shown[2:]) if shown.startswith("~/") else shlex.quote(shown)
    program = "prophound" if argv[:1] == ["predict"] else "python -m src.predict_genome"
    return " ".join([program, *map(arg, argv)])


def run_metadata(opts: PredictOptions, run_name: str, device: str, prepared: PreparedInput, input_summary: dict,
                 windows_df: pd.DataFrame,
                 calls: pd.DataFrame, short_meta: dict, call_meta: dict, refine_meta: dict | None,
                 run_info: dict, outputs: List[str]) -> dict:
    """The run_metadata.json contents. Paths show the home directory as ~."""
    many = opts.n_checkpoints > 1
    checkpoints = [_shown(c) for c in opts.checkpoints]
    metadata = {
        "prophound_version": package_version(),
        "run_name": run_name,
        "config": _shown(opts.config),
        "checkpoint": checkpoints if many else checkpoints[0],
        "ensemble": opts.ensemble if many else "single",
        "input_fasta": _shown(opts.input_fasta),
        "input_sha256": prepared.sha256,
        "input_format": {"format": prepared.format, "compression": prepared.compression,
                         "u_read_as_t": prepared.u_as_t,
                         **({"contig_names": prepared.contig_names} if prepared.contig_names else {})},
        "mode": opts.mode,
        "threshold": opts.threshold,
        "seed_threshold": opts.seed_threshold,
        "extend_threshold": opts.extend_threshold,
        "device": device,
        "window_len": opts.window_len,
        "stride": opts.stride,
        "trained_window_len": opts.trained_window_len,
        "trained_stride": opts.trained_stride,
        "ambiguity_policy": opts.ambiguity_policy,
        "grid_offset": opts.grid_offset,
        "n_windows": int(len(windows_df)),
        "n_regions": int(len(calls)),
        "input": input_summary,
        "n_gap_filter": call_meta["n_gap_filter"],
        "edge_snap": call_meta["edge_snap"],
        "checkpoint_sha256": {_shown(c): sha256(c) for c in opts.checkpoints},
        "outputs": outputs,
        "short_contigs": short_meta,
    }
    if "seed_support" in call_meta:
        metadata["seed_support"] = call_meta["seed_support"]
    metadata["region_builder"] = "coverage_fraction"
    metadata["boundary_refinement"] = "none_bin_edges"
    metadata["coverage_frac"] = opts.coverage_frac
    if opts.max_gap_bp > 0:
        metadata["fragment_merge"] = {"max_gap_bp": opts.max_gap_bp, "bridge_prob": opts.bridge_prob}
    if refine_meta:
        metadata["refine"] = refine_meta
    metadata.update(run_info)
    metadata["command"] = _command_line(sys.argv[1:])
    metadata["platform"] = {"python": platform.python_version(), "system": platform.system(),
                            "release": platform.release(), "machine": platform.machine(), "torch": torch.__version__}
    return metadata


class StagedOutputs:
    """Write each output under a temporary name in the outdir, then rename them all at
    once, so a failed run leaves no partial set (and, with --force, no stale file)."""

    def __init__(self, outdir: Path, force: bool):
        self.outdir, self.force = outdir, force
        self.staged: Dict[str, Path] = {}

    def path(self, name: str) -> Path:
        self.staged[name] = self.outdir / f".{name}.prophound-tmp"
        return self.staged[name]

    def commit(self) -> None:
        if self.force:
            for name in OUTPUT_FILES:
                (self.outdir / name).unlink(missing_ok=True)
        for name, tmp in self.staged.items():
            os.replace(tmp, self.outdir / name)

    def discard(self) -> None:
        for tmp in self.staged.values():
            tmp.unlink(missing_ok=True)


def write_outputs(opts: PredictOptions, run_name: str, device: str, prepared: PreparedInput, input_summary: dict,
                  windows_df: pd.DataFrame,
                  calls: pd.DataFrame, filtered: pd.DataFrame, fine_windows_df: pd.DataFrame | None,
                  short_meta: dict, call_meta: dict, refine_meta: dict | None, run_info: dict) -> List[str]:
    """Write every output file atomically; returns their names."""
    outdir = Path(opts.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    stage = StagedOutputs(outdir, opts.force)
    try:
        windows_df.drop(columns="fasta_path").to_csv(stage.path("windows.tsv"), sep="\t", index=False)
        calls.to_csv(stage.path("regions.tsv"), sep="\t", index=False)
        filtered.to_csv(stage.path("regions.filtered.tsv"), sep="\t", index=False)
        if fine_windows_df is not None:
            fine_windows_df.drop(columns="fasta_path", errors="ignore").to_csv(
                stage.path("fine_windows.tsv"), sep="\t", index=False)
        write_bed(calls, stage.path("regions.bed"))
        write_gff3(calls, stage.path("regions.gff3"), opts.mode, package_version())
        if opts.write_region_fasta:
            write_region_fasta(calls, prepared.fasta, stage.path("regions.fna"))
        outputs = [*stage.staged, "run_metadata.json"]
        metadata = run_metadata(opts, run_name, device, prepared, input_summary, windows_df, calls, short_meta,
                                call_meta, refine_meta, run_info, outputs)
        stage.path("run_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
        stage.commit()
    except BaseException:
        stage.discard()
        raise
    return outputs


def _plural(n: int, word: str, plural: str | None = None) -> str:
    return f"{n:,} {word if n == 1 else (plural or word + 's')}"


def _duration(seconds: float) -> str:
    if seconds < 60:
        return f"{seconds:.1f} s"
    minutes, sec = divmod(int(round(seconds)), 60)
    if minutes < 60:
        return f"{minutes} min {sec} s"
    hours, minutes = divmod(minutes, 60)
    return f"{hours} h {minutes} min"


def print_run_summary(opts: PredictOptions, prepared: PreparedInput, inp: dict, regions_df: pd.DataFrame,
                      short_meta: dict, support_meta: dict | None, gap_meta: dict, device: str,
                      elapsed: float) -> None:
    """One human-readable block on stderr: what went in, what came out, what was skipped."""
    lines = [f"prophound {package_version()} | mode {opts.mode} | {_plural(opts.n_checkpoints, 'checkpoint')}",
             f"input: {_input_kind(prepared)}, {_plural(inp['n_contigs'], 'contig')}, {inp['total_bp']:,} bp; "
             f"{_plural(inp['ambiguous_bases'], 'ambiguity code')} ({inp['ambiguous_ppm']:,g} ppm), "
             f"policy {opts.ambiguity_policy}"]
    if short_meta.get("n_short"):
        lines.append(f"short contigs: {short_meta['n_short']:,} under the window length "
                     f"({short_meta['bp_short']:,} bp); {short_meta['n_too_short']:,} never scanned "
                     f"({short_meta['bp_too_short']:,} bp)")
    n = len(regions_df)
    bp = int((regions_df.end - regions_df.start).sum()) if n else 0
    lines.append(f"calls: {n:,} ({bp:,} bp, {100 * bp / max(1, inp['total_bp']):.1f}% of input)")
    if support_meta:
        k, total = support_meta["min_seed_support"], support_meta["n_checkpoints"]
        dropped = support_meta["calls_dropped_by_vote"]
        if k > 1:
            lines.append(f"seed vote: kept calls made by at least {k} of {total} checkpoints; "
                         f"{_plural(dropped, 'call')} set aside in regions.filtered.tsv")
        else:
            lines.append("seed vote: off, every call kept")
        by_support = support_meta["calls_by_support"]
        if len(by_support) > 1:
            lines.append("calls by checkpoint support: " + ", ".join(
                f"{s} of {total}: {c:,}" for s, c in sorted(by_support.items(), key=lambda kv: -int(kv[0]))))
    if gap_meta["calls_dropped"] or gap_meta["calls_with_any_n"]:
        parts = []
        if gap_meta["calls_dropped"]:
            parts.append(f"{_plural(gap_meta['calls_dropped'], 'call')} set aside as at least "
                         f"{gap_meta['max_call_n_frac']:.0%} N")
        if gap_meta["calls_with_any_n"]:
            n_any = gap_meta["calls_with_any_n"]
            parts.append(f"{_plural(n_any, 'kept call')} {'contains' if n_any == 1 else 'contain'} some N "
                         "(see n_frac)")
        lines.append("N gaps: " + "; ".join(parts))
    lines.append("calls are candidates for follow-up, not annotations")
    lines.append(f"ran on {device} in {_duration(elapsed)}")
    print("\n".join("[summary] " + line for line in lines), file=sys.stderr)


def _input_kind(prepared: PreparedInput) -> str:
    kind = {"fasta": "FASTA", "genbank": "GenBank", "embl": "EMBL", "gff3": "GFF3 with sequence"}[prepared.format]
    return f"{kind} ({prepared.compression})" if prepared.compression else kind


def _out_of_memory(e: BaseException) -> bool:
    return isinstance(e, getattr(torch.cuda, "OutOfMemoryError", ())) or "out of memory" in str(e).lower()


def run_cli(args: argparse.Namespace) -> int:
    """Run `predict` and turn failures into one `error:` line; returns the exit code."""
    verbose = bool(getattr(args, "verbose", False))
    try:
        run_predict(args)
        return 0
    except KeyboardInterrupt:
        print("error: interrupted", file=sys.stderr)
        return 130
    except PredictError as e:
        print(f"error: {e}", file=sys.stderr)
        return 1
    except Exception as e:
        if verbose:
            traceback.print_exc()
        if _out_of_memory(e):
            print(f"error: the device ran out of memory; try a smaller --batch-size (it was {args.batch_size})",
                  file=sys.stderr)
        else:
            print(f"error: unexpected {type(e).__name__}: {e}"
                  + ("" if verbose else " (rerun with --verbose for the traceback)"), file=sys.stderr)
        return 1


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run prophage inference on a genome FASTA (same as `prophound predict`).",
        epilog=PREDICT_EPILOG, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_predict_arguments(parser)
    return parser


def main(argv: List[str] | None = None) -> None:
    # `python -m src.predict_genome` launchers spell the input `--input_fasta`
    argv = sys.argv[1:] if argv is None else list(argv)
    argv = ["--input" + a[len("--input_fasta"):] if a == "--input_fasta" or a.startswith("--input_fasta=") else a
            for a in argv]
    sys.exit(run_cli(build_arg_parser().parse_args(argv)))


if __name__ == "__main__":
    main()
