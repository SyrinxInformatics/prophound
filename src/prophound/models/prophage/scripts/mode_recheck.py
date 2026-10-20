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
"""Operating modes under the shipped pipeline, rebuilt offline from score tracks.

For each run directory (`<genome>/windows.tsv` + `run_metadata.json`) and each
mode, calls are rebuilt the way predict_genome builds them (with runs labelled
r1, r2, r3 and ens, an extra `ens_shipped` row applies each mode's seed vote,
MODE_MIN_SEED_SUPPORT; the other rows are unvoted): full-length windows
with the mode's seed/extend, coverage builder and fragment merging; short-contig
windows (shorter than the run's window_len) per scale at
SHORT_CONTIG_THRESHOLDS[mode], no merging. Refinement is not applied (off by
default).

With --elements: covered and covered50 recall (overall and per stratum), calls
per genome, the on-target fraction of calls (a cost proxy, not precision), and
false calls per Mb of element-free sequence. Without it (a geNomad-negative
sample): calls per genome and per Mb, the fraction of bp called, and the calls'
max_prob -- which cannot distinguish false positives from elements geNomad
missed.

Example:
  python scripts/mode_recheck.py --elements scripts/eval_confirm18_elements.csv \\
      --genomes_dir ../raw/genomes --run ens=analysis/confirm18/ens --out analysis/confirm18/modes.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from compare_variants import score_variant  # noqa: E402
from refine_sweep import covered50  # noqa: E402

try:  # repo layout
    from prophound.models.prophage.src.calls import WINDOW_COORDS, CallRules, build_pass, concat_nonempty, finalize
    from prophound.models.prophage.src.predict_genome import (
        DEFAULT_BRIDGE_PROB,
        DEFAULT_COVERAGE_FRAC,
        DEFAULT_EDGE_SNAP_FRAC,
        DEFAULT_MAX_GAP_BP,
        MODE_MIN_SEED_SUPPORT,  # noqa: F401 -- used as mr.MODE_MIN_SEED_SUPPORT by other scripts
        MODE_THRESHOLDS,
        SHORT_CONTIG_THRESHOLDS,
        apply_seed_vote,
        build_region_table_coverage,
        call_n_fraction,
        iter_fasta_lengths,
    )
except ImportError:  # HPCC data_refactored/ layout
    from src.calls import WINDOW_COORDS, CallRules, build_pass, concat_nonempty, finalize
    from src.predict_genome import (
        DEFAULT_BRIDGE_PROB,
        DEFAULT_COVERAGE_FRAC,
        DEFAULT_EDGE_SNAP_FRAC,
        DEFAULT_MAX_GAP_BP,
        MODE_MIN_SEED_SUPPORT,  # noqa: F401 -- used as mr.MODE_MIN_SEED_SUPPORT by other scripts
        MODE_THRESHOLDS,
        SHORT_CONTIG_THRESHOLDS,
        apply_seed_vote,
        build_region_table_coverage,
        call_n_fraction,
        iter_fasta_lengths,
    )


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


def calls(tracks: dict, mode: str) -> pd.DataFrame:
    seed = MODE_THRESHOLDS[mode]
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


def vote(reg: pd.DataFrame, seed_regs: list, k: int) -> pd.DataFrame:
    """predict_genome's call-level seed vote, applied genome by genome."""
    if k <= 1 or reg.empty:
        return reg
    kept = [apply_seed_vote(sub.reset_index(drop=True), [s[s.genome_id == g] for s in seed_regs], k)[0]
            for g, sub in reg.groupby("genome_id", sort=False)]
    return pd.concat(kept, ignore_index=True)


_warned_no_fasta = False


def shipped_calls(ens_tracks: dict, seed_tracks: list, mode: str, genomes_dir: str | None = None,
                  with_filtered: bool = False):
    """Calls as `predict` ships them for `mode`, built by the CLI's own pipeline (src/calls.py).

    The ensemble track is `ens_tracks`; the vote counts the single-seed tracks in
    `seed_tracks`, aligned to it by window. With `genomes_dir` (holding <genome>.fna)
    every step runs, including the edge snap and the N filter, and the calls carry
    `contig_edge`, `n_frac` and `seed_support` as in regions.tsv. Without it those two
    steps cannot run (they need contig lengths and sequence) and a warning says so;
    the calls are then the CLI's before snapping. A `short` column marks short-contig
    calls. With `with_filtered`, also returns the candidates a filter removed.
    """
    global _warned_no_fasta
    if genomes_dir is None and not _warned_no_fasta:
        print("[mode_recheck] no --genomes_dir: edge snap and the N filter are skipped, so calls can "
              "differ from the CLI's on fragmented or gappy input", file=sys.stderr)
        _warned_no_fasta = True
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
        if full.any():
            passes.append((build_pass(w[full].reset_index(drop=True), [p[full] for p in per_ckpt],
                                      rules.seed_threshold, rules.extend_threshold, rules,
                                      max_gap_bp=rules.max_gap_bp, bridge_prob=rules.bridge_prob), False))
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
            fasta = str(Path(genomes_dir) / f"{g}.fna")
            per_seed = [[p.per_seed[k] for p, _ in passes if p.per_seed] for k in range(len(seed_tracks))] \
                if len(seed_tracks) > 1 else []
            final = finalize(call_g, filt_g, per_seed, dict(iter_fasta_lengths(fasta)),
                             lambda df, f=fasta: call_n_fraction(df, f), rules, DEFAULT_EDGE_SNAP_FRAC)
            call_g, filt_g = final.calls, final.filtered
        kept.append(call_g)
        removed.append(filt_g)
    cols = ["genome_id", "contig_id", "start", "end", "max_prob", "short"]
    out = concat_nonempty(kept) if kept else pd.DataFrame(columns=cols)
    if not with_filtered:
        return out
    return out, (concat_nonempty(removed) if removed else pd.DataFrame(columns=[*cols, "filter"]))


def element_overlap_bp(ev: pd.DataFrame, r) -> int:
    x = ev[(ev.accession == r.genome_id) & (ev.contig == r.contig_id) & (ev.end > r.start) & (ev.start < r.end)]
    return int(sum(min(e, r.end) - max(s, r.start) for s, e in zip(x.start, x.end)))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--elements", default=None, help="omit for a geNomad-negative sample")
    ap.add_argument("--genomes_dir", default=None, help="<genome>.fna, for sequence lengths")
    ap.add_argument("--lengths_csv", default=None, help="alternative: accession,total_bp (e.g. summary.csv)")
    ap.add_argument("--run", action="append", required=True, metavar="LABEL=DIR")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    if not (args.genomes_dir or args.lengths_csv):
        raise SystemExit("need --genomes_dir or --lengths_csv")
    lengths = pd.read_csv(args.lengths_csv).set_index("accession").total_bp if args.lengths_csv else None
    ev = pd.read_csv(args.elements) if args.elements else None
    res, flat = {}, []
    specs = [s.partition("=")[::2] for s in args.run]
    loaded = {label: load_tracks(d) for label, d in specs}
    seeds = [loaded[s] for s in ("r1", "r2", "r3") if s in loaded]
    jobs = [(label, loaded[label], None) for label, _ in specs]
    if "ens" in loaded and len(seeds) == 3:   # the shipped pipeline: ensemble + the mode's vote
        jobs.append(("ens_shipped", loaded["ens"], seeds))
    for label, tracks, vote_seeds in jobs:
        if lengths is not None:
            total_bp = float(lengths.loc[list(tracks)].sum())
        else:
            total_bp = float(sum(n for g in tracks for _, n in iter_fasta_lengths(str(Path(args.genomes_dir) / f"{g}.fna"))))
        el_bp = float((ev.end - ev.start).sum()) if ev is not None else 0.0
        free_bp = total_bp - el_bp
        res[label] = {"n_genomes": len(tracks), "total_bp": total_bp, "modes": {}}
        for mode in MODE_THRESHOLDS:
            reg = (shipped_calls(tracks, vote_seeds, mode, args.genomes_dir) if vote_seeds
                   else calls(tracks, mode))
            n = len(reg)
            m = {"calls_per_genome": n / len(tracks), "short_contig_calls": int(reg["short"].sum()) if n else 0,
                 "calls_per_mb": n / (total_bp / 1e6)}
            if ev is not None:
                ov = [element_overlap_bp(ev, r) for r in reg.itertuples()]
                n_false = sum(o == 0 for o in ov)
                m.update({"covered_recall": score_variant(ev, reg)["covered_recall"] if n else 0.0,
                          "covered50_recall": covered50(ev, reg),
                          "on_target": 1 - n_false / max(1, n),
                          "false_calls_per_mb": n_false / (free_bp / 1e6),
                          "by_stratum": {str(s): {"n": int(len(sub)), "covered50_recall": covered50(sub, reg)}
                                         for s, sub in ev.groupby("stratum")}})
            else:
                bp = float((reg.end - reg.start).sum()) if n else 0.0
                m.update({"frac_bp_called": bp / total_bp,
                          "genomes_with_calls": int(reg.genome_id.nunique()) if n else 0,
                          "max_prob_median": float(reg.max_prob.median()) if n else None,
                          "calls_by_genome": reg.genome_id.value_counts().to_dict() if n else {}})
            res[label]["modes"][mode] = m
            flat.append({"run": label, "mode": mode, **{k: v for k, v in m.items() if not isinstance(v, dict)}})
            print(f"[{label} {mode}] " + " ".join(f"{k} {v:.3f}" if isinstance(v, float) else f"{k} {v}"
                                                  for k, v in m.items() if not isinstance(v, dict)), flush=True)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(res, indent=2) + "\n")
    pd.DataFrame(flat).to_csv(Path(args.out).with_suffix(".csv"), index=False)
    print(f"[done] wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
