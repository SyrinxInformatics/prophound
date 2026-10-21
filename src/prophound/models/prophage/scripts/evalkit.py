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
"""Shared code for the offline scorers and analyses (docs/tech-debt-audit.md X2).

Every script imports measurement code from here, not from another script, so changing one
analysis cannot silently change another's result. Scripts in one investigation may still
import from each other (the g5/g5dev, lambda, short-contig, N-gap and eval-sample families).

- shipped code: imported once, from the installed package or from the ARC source tree
- tracks and calls: `load_tracks`, `calls`, `calls_at`, `vote`, `shipped_calls` (with `refine_pass`,
  `StoredScorer` and the shipped att boundaries when fine tracks and sequence are given)
- refinement offline: `load_run`, `fine_regions_by_genome`, `refine_all`
- intervals and coverage: `overlap`, `covered_bp`, `covered50`, `cov_fracs`, `cov_frac`,
  `split_frag`, `score_variant`, `element_overlap_bp`
- files: `load_genome_list`, `calls_file`, `load_pred_regions`
- statistics: `wilson`
- k-mers: `canon_kmers`
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# The shipped code, imported once for every script: from the installed package, or from the
# ARC source tree when a script runs as `python scripts/<name>.py` from the model directory.
# The names scripts use are re-exported, so no script needs this fallback itself.
try:
    from prophound.models.prophage.src import predict_genome as pg
    from prophound.models.prophage.src.att import Genome, finish_boundaries
    from prophound.models.prophage.src.calls import (
        WINDOW_COORDS,
        CallRules,
        build_pass,
        concat_nonempty,
        finalize,
    )
    from prophound.models.prophage.src.multiscale import refine_regions
    from prophound.models.prophage.src.predict_genome import (  # noqa: F401 -- re-exported
        DEFAULT_BRIDGE_PROB,
        DEFAULT_COVERAGE_FRAC,
        DEFAULT_EDGE_SNAP_FRAC,
        DEFAULT_MAX_GAP_BP,
        MAX_CALL_N_FRAC,
        MODE_MIN_SEED_SUPPORT,
        MODE_THRESHOLDS,
        SHORT_CONTIG_THRESHOLDS,
        apply_seed_vote,
        build_region_table_coverage,
        call_n_fraction,
        encode_window,
        iter_fasta_lengths,
        seed_support,
        snap_contig_edges,
    )
    from prophound.models.prophage.src.predict_genome import refine as _shipped_refine
except ImportError:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from src.att import Genome, finish_boundaries
    from src.calls import WINDOW_COORDS, CallRules, build_pass, concat_nonempty, finalize
    from src.multiscale import refine_regions
    from src.predict_genome import (  # noqa: F401 -- re-exported
        DEFAULT_BRIDGE_PROB,
        DEFAULT_COVERAGE_FRAC,
        DEFAULT_EDGE_SNAP_FRAC,
        DEFAULT_MAX_GAP_BP,
        MAX_CALL_N_FRAC,
        MODE_MIN_SEED_SUPPORT,
        MODE_THRESHOLDS,
        SHORT_CONTIG_THRESHOLDS,
        apply_seed_vote,
        build_region_table_coverage,
        call_n_fraction,
        encode_window,
        iter_fasta_lengths,
        seed_support,
        snap_contig_edges,
    )
    from src.predict_genome import refine as _shipped_refine

    from src import predict_genome as pg  # noqa: F401 -- re-exported

SEP = "__"          # shred_genomes.py names a fragment <contig>__<offset>
K = 21              # k-mer length for containment
BASES = np.full(256, 4, dtype=np.uint8)
for _i, _b in enumerate(b"ACGT"):
    BASES[_b] = _i
    BASES[_b + 32] = _i


# --- tracks and calls


def load_tracks(d: str, prob_col: str = "prob") -> dict:
    """`prob_col` picks the score column: `prob` (the run's own, ensemble if it
    had several checkpoints) or `prob_ckpt<k>` from --write-seed-scores runs."""
    out = {}
    for g in sorted(Path(d).iterdir()):
        w, m = g / "windows.tsv", g / "run_metadata.json"
        if w.exists() and m.exists():
            wl = int(json.loads(m.read_text())["window_len"])
            # round_trip: the default parser can be 1 ulp off the score the CLI wrote
            t = pd.read_csv(w, sep="\t", usecols=["genome_id", "contig_id", "start", "end", prob_col],
                            float_precision="round_trip")
            out[g.name] = (t.rename(columns={prob_col: "prob"}), wl)
    if not out:
        raise SystemExit(f"no <genome>/windows.tsv + run_metadata.json under {d}")
    return out


def calls(tracks: dict, mode: str, seed: float | None = None) -> pd.DataFrame:
    """Unvoted calls per genome at `mode`'s thresholds (`seed` overrides its seed threshold;
    short contigs keep the mode's own), rebuilt as `predict` builds them."""
    seed = MODE_THRESHOLDS[mode] if seed is None else seed
    extend = min(seed, 0.99)   # resolve_thresholds' default
    regs = []
    for w, wl in tracks.values():
        w = w.assign(fasta_path="x")
        span = w.end - w.start
        if (span == wl).any():
            _, r = build_region_table_coverage(w[span == wl], seed_threshold=seed, extend_threshold=extend,
                                               min_coverage_frac=DEFAULT_COVERAGE_FRAC,
                                               max_gap_bp=DEFAULT_MAX_GAP_BP, bridge_prob=DEFAULT_BRIDGE_PROB)
            regs.append(r.assign(short=False))
        for scale, sw in w[span < wl].groupby(span[span < wl]):
            t = SHORT_CONTIG_THRESHOLDS[mode][int(scale)]
            _, r = build_region_table_coverage(sw, seed_threshold=t, extend_threshold=t,
                                               min_coverage_frac=DEFAULT_COVERAGE_FRAC, max_gap_bp=0)
            regs.append(r.assign(short=True))
    regs = [r for r in regs if len(r)]
    return pd.concat(regs, ignore_index=True) if regs else pd.DataFrame(
        columns=["genome_id", "contig_id", "start", "end", "max_prob", "short"])


def calls_at(tracks: dict, t: float) -> pd.DataFrame:
    """`calls` at seed threshold `t`; short contigs use the thresholds of the mode whose seed
    is the largest not exceeding `t` (exploratory below 0.95)."""
    mode = max((m for m, s in MODE_THRESHOLDS.items() if s <= t), key=MODE_THRESHOLDS.get, default="exploratory")
    return calls(tracks, mode, seed=t)


def vote(reg: pd.DataFrame, seed_regs: list, k: int) -> pd.DataFrame:
    """predict_genome's call-level seed vote, applied genome by genome."""
    if k <= 1 or reg.empty:
        return reg
    kept = [apply_seed_vote(sub.reset_index(drop=True), [s[s.genome_id == g] for s in seed_regs], k)[0]
            for g, sub in reg.groupby("genome_id", sort=False)]
    return pd.concat(kept, ignore_index=True)


_warned_no_fasta = False
_warned_no_fine = False


class _NoProgress:
    def add(self, n: int) -> None:
        pass


class StoredScorer:
    """Stands in for predict's Scorer: fine-window scores come from a stored fine_windows.tsv."""

    def __init__(self, fine: pd.DataFrame):
        self.progress = _NoProgress()
        self.lookup = fine.set_index(["contig_id", "start", "end"]).prob

    def score(self, windows: pd.DataFrame, window_len: int):
        prob = self.lookup.reindex(pd.MultiIndex.from_frame(windows[["contig_id", "start", "end"]])).to_numpy()
        if np.isnan(prob).any():
            raise SystemExit(f"{int(np.isnan(prob).sum())} fine windows missing from the stored run")
        return prob, None


def predict_options(mode: str, extra: list | None = None):
    """PredictOptions from the real `predict` parser, at its defaults (plus `extra` arguments)."""
    import argparse

    try:
        from prophound.models.prophage.src.options import PredictOptions, add_predict_arguments
    except ImportError:
        from src.options import PredictOptions, add_predict_arguments
    p = argparse.ArgumentParser()
    add_predict_arguments(p)
    args = p.parse_args(["--input", "unused.fna", "--outdir", "unused", "--mode", mode, *(extra or [])])
    return PredictOptions.from_args(args, trained_window_len=10000, trained_stride=1000)


def refine_pass(calls: pd.DataFrame, fine: pd.DataFrame, lengths: dict, genome_id: str, fasta: str,
                mode: str, extra: list | None = None) -> pd.DataFrame:
    """The shipped refine() on coarse calls, fed stored fine scores (an Exploratory refinement
    run's fine windows cover every mode's neighbourhoods), at `predict`'s defaults (plus `extra`
    predict arguments, e.g. ["--refine-threshold", "0.9"])."""
    if calls.empty:
        return calls
    return _shipped_refine(StoredScorer(fine), calls, lengths, genome_id, fasta, predict_options(mode, extra))[0]


def shipped_calls(ens_tracks: dict, seed_tracks: list, mode: str, genomes_dir: str | None = None,
                  with_filtered: bool = False, fine_tracks: dict | None = None):
    """Calls as `predict` ships them for `mode`, built by the CLI's own pipeline (src/calls.py).

    The ensemble track is `ens_tracks`; the vote counts the single-seed tracks in
    `seed_tracks`, aligned to it by window. With `genomes_dir` (holding <genome>.fna)
    every step runs, including the edge snap and the N filter, and the calls carry
    `contig_edge`, `n_frac` and `seed_support` as in regions.detailed.tsv. Without it those two
    steps cannot run (they need contig lengths and sequence) and a warning says so;
    the calls are then the CLI's before snapping. A `short` column marks short-contig
    calls. With `with_filtered`, also returns the candidates a filter removed.

    Since 0.1.0 `predict` refines calls and places att boundaries by default
    (docs/att-boundaries.md). `fine_tracks` ({genome: fine windows}, from refinement-on
    runs) adds the refinement, and `genomes_dir` the att boundaries and edge ranges.
    Without `fine_tracks` a warning says the calls are the pre-refinement ones.
    """
    global _warned_no_fasta, _warned_no_fine
    if genomes_dir is None and not _warned_no_fasta:
        print("[mode_recheck] no --genomes_dir: edge snap, the N filter and att boundaries are skipped, so "
              "calls can differ from the CLI's", file=sys.stderr)
        _warned_no_fasta = True
    if fine_tracks is None and not _warned_no_fine:
        print("[evalkit] no fine tracks: calls are rebuilt without the refinement pass, which predict runs "
              "by default since 0.1.0 (docs/att-boundaries.md)", file=sys.stderr)
        _warned_no_fine = True
    kept, removed = [], []
    for g, (w, wl) in ens_tracks.items():
        per_ckpt = []
        for t in seed_tracks:
            sw = w[WINDOW_COORDS].merge(t[g][0][[*WINDOW_COORDS, "prob"]], on=WINDOW_COORDS, how="left")
            if sw.prob.isna().any():
                raise SystemExit(f"{g}: a seed track lacks windows the ensemble track has")
            per_ckpt.append(sw.prob.to_numpy())
        rules = CallRules.for_mode(mode, wl, len(seed_tracks))
        span = (w.end - w.start).to_numpy()
        passes = []
        full = span == wl
        fasta = str(Path(genomes_dir) / f"{g}.fna") if genomes_dir is not None else None
        if full.any():
            coarse = build_pass(w[full].reset_index(drop=True), [p[full] for p in per_ckpt],
                                rules.seed_threshold, rules.extend_threshold, rules,
                                max_gap_bp=rules.max_gap_bp, bridge_prob=rules.bridge_prob)
            if fine_tracks is not None and g in fine_tracks:
                if fasta is None:
                    raise SystemExit("refinement offline needs genomes_dir")
                coarse.calls = refine_pass(coarse.calls, fine_tracks[g], dict(iter_fasta_lengths(fasta)), g,
                                           fasta, mode)
            passes.append((coarse, False))
        for scale in sorted(set(span[span < wl])):   # ascending, as predict groups them
            m = span == scale
            t = rules.short_thresholds[int(scale)]
            passes.append((build_pass(w[m].reset_index(drop=True), [p[m] for p in per_ckpt], t, t, rules,
                                      max_gap_bp=0, kind_suffix=f"_short{int(scale)}"), True))
        if not passes:
            continue
        call_g = concat_nonempty([p.calls.assign(short=short) for p, short in passes])
        filt_g = concat_nonempty([p.filtered.assign(short=short) for p, short in passes])
        if genomes_dir is not None:
            per_seed = [[p.per_seed[k] for p, _ in passes if p.per_seed] for k in range(len(seed_tracks))] \
                if len(seed_tracks) > 1 else []
            lengths = dict(iter_fasta_lengths(fasta))
            final = finalize(call_g, filt_g, per_seed, lengths,
                             lambda df, f=fasta: call_n_fraction(df, f), rules, DEFAULT_EDGE_SNAP_FRAC)
            genome = Genome.from_fasta(fasta)
            call_g = finish_boundaries(final.calls, genome, lengths, mode)
            filt_g = finish_boundaries(final.filtered, genome, lengths, mode)
        kept.append(call_g)
        removed.append(filt_g)
    cols = ["genome_id", "contig_id", "start", "end", "max_prob", "short"]
    out = concat_nonempty(kept) if kept else pd.DataFrame(columns=cols)
    if not with_filtered:
        return out
    return out, (concat_nonempty(removed) if removed else pd.DataFrame(columns=[*cols, "filter"]))


# --- refinement, offline


def load_run(d: str):
    """(coarse regions, fine windows) per genome, rebuilt from the run's metadata."""
    coarse, fine, meta0 = [], [], None
    for gdir in sorted(Path(d).iterdir()):
        w, f, m = gdir / "windows.tsv", gdir / "fine_windows.tsv", gdir / "run_metadata.json"
        if not (w.exists() and m.exists()):
            continue
        meta = json.loads(m.read_text())
        meta0 = meta0 or meta
        wdf = pd.read_csv(w, sep="\t")
        # short-contig windows (--short-contigs scan) have their own thresholds and
        # are never refined; mode_recheck.py scores them
        wdf = wdf[(wdf.end - wdf.start) == int(meta["window_len"])]
        fm = meta.get("fragment_merge", {})
        _, reg = build_region_table_coverage(
            wdf, seed_threshold=meta["seed_threshold"], extend_threshold=meta["extend_threshold"],
            min_coverage_frac=meta["coverage_frac"], max_gap_bp=fm.get("max_gap_bp", 0),
            bridge_prob=fm.get("bridge_prob", 0.0))
        coarse.append(reg)
        if f.exists():
            fine.append(pd.read_csv(f, sep="\t", usecols=["genome_id", "contig_id", "start", "end", "prob"]))
    if not coarse:
        raise SystemExit(f"no <genome>/windows.tsv + run_metadata.json under {d}")
    fine_df = pd.concat(fine, ignore_index=True) if fine else pd.DataFrame(
        columns=["genome_id", "contig_id", "start", "end", "prob"])
    return pd.concat(coarse, ignore_index=True), fine_df, meta0, len(coarse), len(fine)


def fine_regions_by_genome(fine: pd.DataFrame, thr: float, frac: float) -> dict:
    """genome -> (fine windows, fine regions) at one threshold; reused across variants."""
    out = {}
    for g, fg in fine.groupby("genome_id", sort=False):
        out[g] = build_region_table_coverage(fg.assign(fasta_path="x"), seed_threshold=thr,
                                             extend_threshold=thr, min_coverage_frac=frac, max_gap_bp=0)
    return out


def refine_all(coarse: pd.DataFrame, fine_by_genome: dict, alpha: float, gap: int,
               maxp: float, ramp=None) -> pd.DataFrame:
    out = []
    for g, cg in coarse.groupby("genome_id", sort=False):
        if g not in fine_by_genome:
            out.append(cg)
            continue
        fg, fr = fine_by_genome[g]
        out.append(refine_regions(cg, fg, fr, edge_alpha=alpha, split_min_gap_bp=gap,
                                  split_max_prob=maxp, alpha_ramp_bp=ramp))
    return pd.concat(out, ignore_index=True) if out else coarse


# --- intervals and coverage


def overlap(a_start, a_end, b_start, b_end) -> float:
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def covered_bp(t_start: int, t_end: int, calls) -> int:
    """bp of [t_start, t_end) under the union of `calls` (start, end pairs)."""
    ivs = sorted((max(t_start, s), min(t_end, e)) for s, e in calls if min(t_end, e) > max(t_start, s))
    bp, cur_s, cur_e = 0, None, None
    for s, e in ivs:
        if cur_e is None or s > cur_e:
            bp += 0 if cur_e is None else cur_e - cur_s
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    return bp + (0 if cur_e is None else cur_e - cur_s)


def covered50(ev: pd.DataFrame, reg: pd.DataFrame) -> float:
    """Fraction of elements with >= 50% of their bp under some call.

    Covered recall (>= 1 bp) credits grazes -- a neighbouring call's overshoot
    clipping an element's edge. Both losses 817370's refinement showed were
    grazes (docs/downstream-postprocessing.md "Per-element dig"), so every
    variant is also scored on this graze-proof version.
    """
    ok = 0
    for r in ev.itertuples():
        g = reg[(reg.genome_id == r.accession) & (reg.contig_id == r.contig)
                & (reg.end > r.start) & (reg.start < r.end)].sort_values("start")
        bp, cur_s, cur_e = 0, None, None
        for c in g.itertuples():  # union of overlapping calls, clipped to the element
            s, e = max(int(c.start), r.start), min(int(c.end), r.end)
            if cur_e is None or s > cur_e:
                bp += 0 if cur_e is None else cur_e - cur_s
                cur_s, cur_e = s, e
            else:
                cur_e = max(cur_e, e)
        bp += 0 if cur_e is None else cur_e - cur_s
        ok += bp >= 0.5 * (r.end - r.start)
    return ok / max(1, len(ev))


def cov_fracs(ev: pd.DataFrame, reg: pd.DataFrame) -> np.ndarray:
    """Per element, the fraction of its bp under the union of overlapping calls."""
    out = []
    for e in ev.itertuples():
        x = reg[(reg.genome_id == e.accession) & (reg.contig_id == e.contig)
                & (reg.end > e.start) & (reg.start < e.end)]
        iv = sorted((max(int(a), e.start), min(int(b), e.end)) for a, b in zip(x.start, x.end))
        bp, cs, ce = 0, None, None
        for a, b in iv:
            if ce is None or a > ce:
                bp += (ce - cs) if ce is not None else 0
                cs, ce = a, b
            else:
                ce = max(ce, b)
        bp += (ce - cs) if ce is not None else 0
        out.append(bp / (e.end - e.start))
    return np.array(out)


def split_frag(name: str) -> tuple[str, int]:
    """`CP002907__12345` -> (CP002907, 12345); an unshredded name has offset 0."""
    seq, _, off = str(name).partition(SEP)
    return seq, int(off) if off else 0


def cov_frac(truth: pd.DataFrame, calls: pd.DataFrame) -> np.ndarray:
    """Per element, the fraction of its bp under the union of calls on its sequence."""
    out = np.zeros(len(truth))
    by = {k: g.sort_values("start")[["start", "end"]].to_numpy() for k, g in calls.groupby(["assembly", "sequence"])}
    for i, r in enumerate(truth.itertuples()):
        c = by.get((r.assembly, r.sequence))
        if c is None:
            continue
        c = c[(c[:, 1] > r.start) & (c[:, 0] < r.end)]
        bp, cs, ce = 0, None, None
        for s, e in c:
            s, e = max(s, r.start), min(e, r.end)
            if ce is None or s > ce:
                bp += 0 if ce is None else ce - cs
                cs, ce = s, e
            else:
                ce = max(ce, e)
        bp += 0 if ce is None else ce - cs
        out[i] = bp / (r.end - r.start)
    return out


def score_variant(ev: pd.DataFrame, reg: pd.DataFrame) -> dict:
    covered, ious, n_overlaps = [], [], []
    per_call: dict = {}
    for r in ev.itertuples():
        g = reg[(reg.genome_id == r.accession) & (reg.contig_id == r.contig)]
        o = g[(g.end > r.start) & (g.start < r.end)]
        covered.append(len(o) > 0)
        n_overlaps.append(len(o))
        for c in o.itertuples():
            per_call.setdefault((c.genome_id, c.contig_id, c.start, c.end), []).append(
                float(r.kb) if hasattr(r, "kb") else 0.0)
        if len(o):
            inter = np.minimum(o.end, r.end) - np.maximum(o.start, r.start)
            union = np.maximum(o.end, r.end) - np.minimum(o.start, r.start)
            ious.append(float((inter / union).max()))
        else:
            ious.append(0.0)
    n_genomes = ev.accession.nunique()
    # Conflation: one call claiming several distinct elements. Covered recall is
    # blind to this by construction (3 elements under 1 call score 3/3), so it is
    # reported explicitly -- see docs/downstream-postprocessing.md "Conflation".
    conflated = {k: v for k, v in per_call.items() if len(v) >= 2}
    absorbed = sum(len(v) for v in conflated.values())
    n_covered = int(np.sum(covered))
    return {
        "covered_recall": float(np.mean(covered)),
        "n_covered": n_covered,
        # matched recall: covered, minus elements that merely share a call with
        # another. The 1-1 achievable ceiling; the metric that *does* see conflation.
        "matched_recall": float((n_covered - (absorbed - len(conflated))) / max(1, len(ev))),
        "n_elements": int(len(ev)),
        "n_calls": int(len(reg)),
        "calls_per_genome": float(len(reg) / max(1, n_genomes)),
        "median_iou": float(np.median(ious)),
        "iou_ge_50": float(np.mean(np.asarray(ious) >= 0.5)),
        "fragmented_elements": int(np.sum(np.asarray(n_overlaps) > 1)),
        "conflated_calls": len(conflated),
        "absorbed_elements": absorbed,
        "_covered_vec": [bool(c) for c in covered],
    }


def element_overlap_bp(ev: pd.DataFrame, r) -> int:
    x = ev[(ev.accession == r.genome_id) & (ev.contig == r.contig_id) & (ev.end > r.start) & (ev.start < r.end)]
    return int(sum(min(e, r.end) - max(s, r.start) for s, e in zip(x.start, x.end)))


# --- files


def load_genome_list(path: Path) -> list:
    with open(path) as f:
        return [line.strip() for line in f if line.strip()]


def calls_file(genome_dir: Path) -> Path:
    """predict's calls in 0-based coordinates: regions.detailed.tsv since 0.1.0 (its regions.tsv
    is the 1-based researcher table), regions.tsv in runs from before it."""
    detailed = Path(genome_dir) / "regions.detailed.tsv"
    return detailed if detailed.exists() else Path(genome_dir) / "regions.tsv"


def load_pred_regions(predict_dir: Path, accession: str) -> pd.DataFrame:
    tsv = calls_file(predict_dir / accession)
    # predict_genome.py writes a bare newline (no header) when it finds zero
    # regions for a genome -- pd.read_csv chokes on that with EmptyDataError.
    if not tsv.exists() or tsv.stat().st_size <= 1:
        return pd.DataFrame(columns=["contig_id", "start", "end"])
    df = pd.read_csv(tsv, sep="\t")
    if df.empty:
        return pd.DataFrame(columns=["contig_id", "start", "end"])
    return df[["contig_id", "start", "end"]]


# --- statistics


def wilson(k: int, n: int, z: float = 1.959964) -> list:
    if n == 0:
        return [float("nan"), float("nan")]
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [float(c - h), float(c + h)]


# --- k-mers


def canon_kmers(seq: str, k: int = K) -> np.ndarray:
    """Unique canonical k-mer codes (2 bits/base); k-mers spanning non-ACGT are dropped."""
    arr = BASES[np.frombuffer(seq.encode("ascii", "ignore"), dtype=np.uint8)]
    n = arr.size - k + 1
    if n <= 0:
        return np.empty(0, dtype=np.uint64)
    bad = np.convolve((arr >= 4).astype(np.int32), np.ones(k, dtype=np.int32), "valid") > 0
    a = np.where(arr < 4, arr, 0).astype(np.uint64)
    fwd = np.zeros(n, dtype=np.uint64)
    rev = np.zeros(n, dtype=np.uint64)
    for i in range(k):
        fwd = (fwd << np.uint64(2)) | a[i:i + n]
        rev |= (np.uint64(3) - a[i:i + n]) << np.uint64(2 * i)
    return np.unique(np.minimum(fwd, rev)[~bad])
