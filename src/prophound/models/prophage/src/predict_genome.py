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
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd
import torch
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
    SHORT_CONTIG_SCALES,
    SHORT_CONTIG_STRIDE,
    SHORT_CONTIG_THRESHOLDS,
    PredictOptions,
    add_predict_arguments,
)

WINDOW_COLUMNS = ["genome_id", "contig_id", "start", "end", "fasta_path", "prob", "pred_label", "region_pred_label"]


def warn(msg: str) -> None:
    print(f"[warn] {msg}", file=sys.stderr)


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
) -> np.ndarray:
    model.eval()
    all_logits: List[np.ndarray] = []
    for batch_x, _ in tqdm(loader, desc="predict", leave=False):
        batch_x = batch_x.to(device, non_blocking=True)
        logits = model(batch_x).view(-1)
        all_logits.append(logits.detach().cpu().numpy())
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


def score_windows(
    model: torch.nn.Module,
    windows_df: pd.DataFrame,
    window_len: int,
    opts: PredictOptions,
    device: torch.device,
) -> Tuple[np.ndarray, List[np.ndarray]]:
    """Score `windows_df` with every checkpoint: the ensemble track and each checkpoint's own."""
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
    for ckpt_path in opts.checkpoints:
        checkpoint = torch.load(ckpt_path, map_location=device, weights_only=False)
        model.load_state_dict(checkpoint["model"])
        per_checkpoint.append(np.asarray(predict_probs(model, loader, device), dtype=np.float64))
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


def call_coarse(model, device, windows_df: pd.DataFrame, opts: PredictOptions, rules: CallRules) -> PassCalls:
    """The full-length window pass: score, build calls with fragment merging, vote."""
    if windows_df.empty:   # every contig is short
        empty = pd.DataFrame(columns=REGION_COLUMNS)
        return PassCalls(pd.DataFrame(columns=WINDOW_COLUMNS), empty, empty.assign(filter=pd.Series(dtype=str)), [])
    ensemble, per_checkpoint = score_windows(model, windows_df, opts.window_len, opts, device)
    return build_pass(scored(windows_df, ensemble, per_checkpoint, opts.threshold, opts), per_checkpoint,
                      opts.seed_threshold, opts.extend_threshold, rules,
                      max_gap_bp=opts.max_gap_bp, bridge_prob=opts.bridge_prob)


def refine(model, device, calls: pd.DataFrame, contig_lengths: Dict[str, int], genome_id: str,
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
            neighborhoods, contig_lengths, genome_id, opts.input_fasta,
            window_len=opts.refine_window_len, stride=opts.refine_stride,
        )
        if not fine_windows_df.empty:
            fine_windows_df = fine_windows_df.copy()
            fine_windows_df["prob"] = score_windows(model, fine_windows_df, opts.refine_window_len, opts, device)[0]
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


def call_short_contigs(model, device, short_windows_df: pd.DataFrame, opts: PredictOptions,
                       rules: CallRules) -> Tuple[List[PassCalls], dict]:
    """Contigs shorter than --window-len, one pass per scale (5 or 3 kb), no merging.

    Returns the passes and the metadata additions (thresholds used, etc.); no passes
    when scanning is off or nothing is scannable.
    """
    if opts.short_contigs != "scan" or short_windows_df.empty:
        return [], {}
    if opts.short_contig_threshold is not None:
        warn(f"--short-contig-threshold {float(opts.short_contig_threshold):g} overrides the calibrated "
             "short-contig thresholds.")
    passes, used = [], {}
    for scale, sw in short_windows_df.groupby("scale"):
        t = opts.short_threshold(int(scale))
        used[str(int(scale))] = {"seed": t, "extend": t}
        sw = sw.drop(columns="scale").reset_index(drop=True)
        ensemble, per_checkpoint = score_windows(model, sw, int(scale), opts, device)
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
        short_calls = pd.concat([p.calls for p in short], ignore_index=True)
        if "refine_status" in calls.columns:
            short_calls = short_calls.assign(refine_status="short_contig", coarse_start=short_calls["start"],
                                             coarse_end=short_calls["end"], n_pieces=1)
        # concatenating onto an empty frame would change the column dtypes
        calls = pd.concat([calls, short_calls], ignore_index=True) if not calls.empty else short_calls
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


def run_predict(args: argparse.Namespace) -> None:
    """Run inference end to end from parsed `predict` arguments."""
    cfg = load_cfg(args.config)
    device = resolve_device(args.device)
    opts = PredictOptions.from_args(args, cfg.data.window_len, cfg.data.stride)
    if opts.window_len != opts.trained_window_len:
        # the model accepts any length, but its scores on negatives rise as windows
        # shrink, so the thresholds don't transfer (docs/adaptive-multiscale-detection.md)
        warn(f"--window-len {opts.window_len} differs from the trained {opts.trained_window_len}; the mode "
             f"thresholds were calibrated at {opts.trained_window_len} bp and are not validated at this length.")

    input_summary = validate_fasta(opts.input_fasta)   # fails early with a clear message
    contig_lengths = dict(iter_fasta_lengths(opts.input_fasta))
    genome_id = Path(opts.input_fasta).stem
    windows_df = build_windows(opts.input_fasta, genome_id, opts.window_len, opts.stride,
                               offset=opts.grid_offset, contig_lengths=contig_lengths)
    short_windows_df, short_report = build_short_contig_windows(
        opts.input_fasta, genome_id, window_len=opts.window_len, contig_lengths=contig_lengths)
    if short_report["n_short"]:
        # never silent (docs/cli-and-packaging.md "Input robustness")
        verb = "scanned at reduced scale" if opts.short_contigs == "scan" else "NOT scanned (--short-contigs skip)"
        warn(f"{short_report['n_short']} contig(s) shorter than {opts.window_len} bp "
             f"({short_report['bp_short']} bp): {short_report['n_scannable']} {verb}; "
             f"{short_report['n_too_short']} under {min(SHORT_CONTIG_SCALES)} bp never scanned "
             f"({short_report['bp_too_short']} bp).")
    if windows_df.empty and (opts.short_contigs == "skip" or short_windows_df.empty):
        raise RuntimeError("No inference windows were generated from the input FASTA.")

    model = build_model(cfg, device)
    if opts.requested_min_seed_support > 1 and opts.min_seed_support < opts.requested_min_seed_support:
        # never silent: the mode thresholds were set with the vote over 3 checkpoints
        k = opts.min_seed_support
        warn(f"seed vote of {opts.requested_min_seed_support} requested but only {opts.n_checkpoints} "
             f"checkpoint(s): {'voting at ' + str(k) if k > 1 else 'no vote applied'}. The '{opts.mode}' "
             "thresholds were calibrated with the vote over 3 checkpoints; expect more calls.")
    rules = call_rules(opts)

    coarse = call_coarse(model, device, windows_df, opts, rules)
    coarse_calls, fine_windows_df, refine_meta = refine(model, device, coarse.calls, contig_lengths, genome_id, opts)
    short, short_scan_meta = call_short_contigs(model, device, short_windows_df, opts, rules)
    windows_df, calls, filtered = combine_passes(coarse, coarse_calls, short)

    per_seed = [[p.per_seed[k] for p in (coarse, *short) if p.per_seed] for k in range(opts.n_checkpoints)] \
        if opts.n_checkpoints > 1 else []
    final = finalize(calls, filtered, per_seed, contig_lengths,
                     lambda df: call_n_fraction(df, opts.input_fasta), rules, opts.edge_snap_frac,
                     requested_min_seed_support=opts.requested_min_seed_support)
    gap_meta = final.meta["n_gap_filter"]
    if gap_meta["calls_dropped"]:
        warn(f"{gap_meta['calls_dropped']} call(s) dropped: at least {MAX_CALL_N_FRAC:.0%} N "
             "(imputed gap sequence).")
    windows_df = label_windows(windows_df, final.calls)

    short_meta = {"mode": opts.short_contigs, **short_report, **short_scan_meta}
    write_outputs(opts, cfg.run_name, str(device), input_summary, windows_df, final.calls, fine_windows_df,
                  short_meta, final.meta, refine_meta)
    print_run_summary(opts, input_summary, final.calls, short_meta, final.meta.get("seed_support"), gap_meta)


def run_metadata(opts: PredictOptions, run_name: str, device: str, input_summary: dict, windows_df: pd.DataFrame,
                 calls: pd.DataFrame, short_meta: dict, call_meta: dict, refine_meta: dict | None) -> dict:
    """The run_metadata.json contents."""
    many = opts.n_checkpoints > 1
    metadata = {
        "prophound_version": package_version(),
        "run_name": run_name,
        "checkpoint": list(opts.checkpoints) if many else opts.checkpoints[0],
        "ensemble": opts.ensemble if many else "single",
        "input_fasta": opts.input_fasta,
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
        "checkpoint_sha256": {c: sha256(c) for c in opts.checkpoints},
        "outputs": ["windows.tsv", "regions.tsv", "regions.bed", "regions.gff3", "run_metadata.json"]
        + (["regions.fna"] if opts.write_region_fasta else []),
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
    return metadata


def write_outputs(opts: PredictOptions, run_name: str, device: str, input_summary: dict, windows_df: pd.DataFrame,
                  calls: pd.DataFrame, fine_windows_df: pd.DataFrame | None, short_meta: dict, call_meta: dict,
                  refine_meta: dict | None) -> None:
    outdir = Path(opts.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    windows_df.to_csv(outdir / "windows.tsv", sep="\t", index=False)
    calls.to_csv(outdir / "regions.tsv", sep="\t", index=False)
    metadata = run_metadata(opts, run_name, device, input_summary, windows_df, calls, short_meta, call_meta,
                            refine_meta)
    (outdir / "run_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"[done] wrote {outdir / 'windows.tsv'}", file=sys.stderr)
    if fine_windows_df is not None:
        fine_windows_df.to_csv(outdir / "fine_windows.tsv", sep="\t", index=False)
        print(f"[done] wrote {outdir / 'fine_windows.tsv'}", file=sys.stderr)
    print(f"[done] wrote {outdir / 'regions.tsv'}", file=sys.stderr)
    print(f"[done] wrote {outdir / 'run_metadata.json'}", file=sys.stderr)
    write_bed(calls, outdir / "regions.bed")
    write_gff3(calls, outdir / "regions.gff3", opts.mode, package_version())
    print(f"[done] wrote {outdir / 'regions.bed'} and {outdir / 'regions.gff3'}", file=sys.stderr)
    if opts.write_region_fasta:
        write_region_fasta(calls, opts.input_fasta, outdir / "regions.fna")
        print(f"[done] wrote {outdir / 'regions.fna'}", file=sys.stderr)


def print_run_summary(opts: PredictOptions, inp: dict, regions_df: pd.DataFrame, short_meta: dict,
                      support_meta: dict | None, gap_meta: dict) -> None:
    """One human-readable block on stderr: what went in, what came out, what was skipped."""
    lines = [f"prophound {package_version()} | mode {opts.mode} | {opts.n_checkpoints} checkpoint(s)",
             f"input: {inp['n_contigs']} contig(s), {inp['total_bp']:,} bp; "
             f"ambiguity codes {inp['ambiguous_bases']:,} ({inp['ambiguous_ppm']:g} ppm), "
             f"policy {opts.ambiguity_policy}"]
    if short_meta.get("n_short"):
        lines.append(f"short contigs: {short_meta['n_short']} under the window length ({short_meta['bp_short']:,} bp); "
                     f"{short_meta['n_too_short']} never scanned ({short_meta['bp_too_short']:,} bp)")
    n = len(regions_df)
    bp = int((regions_df.end - regions_df.start).sum()) if n else 0
    lines.append(f"calls: {n} ({bp:,} bp, {100 * bp / max(1, inp['total_bp']):.1f}% of input)")
    if support_meta:
        lines.append(f"seed vote: min support {support_meta['min_seed_support']}, "
                     f"{support_meta['calls_dropped_by_vote']} call(s) dropped; kept by support "
                     f"{support_meta['calls_by_support']}")
    if gap_meta["calls_dropped"] or gap_meta["calls_with_any_n"]:
        lines.append(f"N gaps: {gap_meta['calls_dropped']} call(s) dropped as >= {gap_meta['max_call_n_frac']:.0%} N; "
                     f"{gap_meta['calls_with_any_n']} kept call(s) contain some N (see n_frac)")
    lines.append("calls are candidates for follow-up, not annotations")
    print("\n".join("[summary] " + line for line in lines), file=sys.stderr)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run prophage inference on a genome FASTA (same as `prophound predict`).")
    add_predict_arguments(parser)
    return parser


def main(argv: List[str] | None = None) -> None:
    # `python -m src.predict_genome` launchers spell the input `--input_fasta`
    argv = sys.argv[1:] if argv is None else list(argv)
    argv = ["--input" + a[len("--input_fasta"):] if a == "--input_fasta" or a.startswith("--input_fasta=") else a
            for a in argv]
    run_predict(build_arg_parser().parse_args(argv))


if __name__ == "__main__":
    main()
