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
"""Boundary pipeline variants, rebuilt offline from refinement-on runs (docs/att-boundaries.md).

Input: a `predict --mode exploratory --write-seed-scores --refine-window-len ...` batch
directory (<genome>/windows.tsv, fine_windows.tsv, run_metadata.json). Exploratory
calls contain every other mode's calls, so their fine windows cover every mode's
refinement neighbourhoods. For each mode the calls are rebuilt with the shipped
code: build_pass, then the shipped `refine()` (fed the stored fine scores by a
stand-in scorer, with options from the real `predict` parser), then finalize.

Variants per mode:
  shipped        coarse calls (refinement off, as released)
  refine         + the fine pass
  att            coarse calls, att boundaries where found
  refine+att     the proposal: att where found (searched on each call's coarse
                 coordinates, as tuned), refined edges elsewhere

Writes <out>.csv (metrics per sample, mode and variant) and <out>.calls.tsv (the
refine+att calls with their coarse and refined coordinates, for calibrating and
checking uncertainty ranges).

Example:
  python scripts/boundary_eval.py --runs ../../../../arc_data/att_runs/tune42 \\
      --genomes_dir ../../../../arc_data/raw/genomes --elements scripts/eval_sizestrat_elements.csv \\
      --out analysis/att/tune42_pipeline
"""

from __future__ import annotations

import argparse
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

import att_boundaries as ab
import evalkit as ek

REFINE_ARGS = ["--refine-window-len", "5000", "--refine-threshold", "0.6", "--edge-alpha", "0.6",
               "--edge-alpha-ramp", "15000", "60000"]
ATT_RULE = {"w_in": 10000, "w_out": 5000, "max_len": 300, "max_copies": 2, "min_entropy": 1.5,
            "core_slack": None, "e_max": 1e-6}
MODES = ["exploratory", "default", "conservative"]


class _Progress:
    def add(self, n: int) -> None:
        pass


class StoredScorer:
    """Stands in for predict's Scorer: fine-window scores come from a stored fine_windows.tsv."""

    def __init__(self, fine: pd.DataFrame):
        self.progress = _Progress()
        self.lookup = fine.set_index(["contig_id", "start", "end"]).prob

    def score(self, windows: pd.DataFrame, window_len: int):
        key = pd.MultiIndex.from_frame(windows[["contig_id", "start", "end"]])
        prob = self.lookup.reindex(key).to_numpy()
        if np.isnan(prob).any():
            raise SystemExit(f"{int(np.isnan(prob).sum())} fine windows missing from the stored run")
        return prob, None


def predict_opts(mode: str, fasta: str):
    from prophound.models.prophage.src.options import PredictOptions, add_predict_arguments

    p = argparse.ArgumentParser()
    add_predict_arguments(p)
    args = p.parse_args(["--input", fasta, "--outdir", "unused", "--mode", mode, *REFINE_ARGS])
    return PredictOptions.from_args(args, trained_window_len=10000, trained_stride=1000)


def genome_calls(job: dict) -> list:
    """{mode: {variant: calls}} for one genome, built by the shipped code."""
    from prophound.models.prophage.src.predict_genome import refine

    g, d, fasta = job["genome"], Path(job["dir"]), job["fasta"]
    w = pd.read_csv(d / "windows.tsv", sep="\t", float_precision="round_trip")
    fine = pd.read_csv(d / "fine_windows.tsv", sep="\t", float_precision="round_trip")
    wl = int(json.loads((d / "run_metadata.json").read_text())["window_len"])
    lengths = dict(ek.iter_fasta_lengths(fasta))
    seqs, text, p = ab.load_genome(Path(fasta))
    out = []
    for mode in job["modes"]:
        opts = predict_opts(mode, fasta)
        rules = ek.CallRules.for_mode(mode, wl, 3)
        span = (w.end - w.start).to_numpy()
        per = [w[f"prob_ckpt{k}"].to_numpy() for k in (1, 2, 3)]
        full = span == wl
        ww = w[["genome_id", "contig_id", "start", "end", "prob"]]
        coarse = ek.build_pass(ww[full].reset_index(drop=True), [x[full] for x in per], rules.seed_threshold,
                               rules.extend_threshold, rules, max_gap_bp=rules.max_gap_bp,
                               bridge_prob=rules.bridge_prob)
        shorts = []
        for scale in sorted(set(span[span < wl])):
            m = span == scale
            t = rules.short_thresholds[int(scale)]
            shorts.append(ek.build_pass(ww[m].reset_index(drop=True), [x[m] for x in per], t, t, rules,
                                        max_gap_bp=0, kind_suffix=f"_short{int(scale)}"))
        variants = {}
        for name, do_refine in (("shipped", False), ("refine", True)):
            calls = coarse.calls
            if do_refine and len(calls):
                calls, _, _ = refine(StoredScorer(fine), calls, lengths, g, fasta, opts)
            allc = ek.concat_nonempty([calls.assign(short=False), *(s.calls.assign(short=True) for s in shorts)]) \
                if len(calls) or shorts else calls
            filt = ek.concat_nonempty([coarse.filtered, *(s.filtered for s in shorts)])
            per_seed = [[pp.per_seed[k] for pp in (coarse, *shorts) if pp.per_seed] for k in range(3)]
            final = ek.finalize(allc, filt, per_seed, lengths, lambda df: ek.call_n_fraction(df, fasta), rules,
                                ek.DEFAULT_EDGE_SNAP_FRAC)
            variants[name] = final.calls.assign(genome_id=g)
        # att: searched on each call's coarse coordinates, the frame the rule was tuned in
        for name in ("shipped", "refine"):
            c = variants[name].copy()
            cs = c["coarse_start"] if "coarse_start" in c else c["start"]
            ce = c["coarse_end"] if "coarse_end" in c else c["end"]
            pieces = c["n_pieces"] if "n_pieces" in c else pd.Series(1, index=c.index)
            att_s, att_e, att_l = [], [], []
            for i, r in enumerate(c.itertuples()):
                hit = None
                if not r.short and r.contig_edge == "none" and int(pieces.iloc[i]) == 1:
                    hit = att_for(seqs[r.contig_id], text, int(cs.iloc[i]), int(ce.iloc[i]), p)
                att_s.append(hit["att_start"] if hit else np.nan)
                att_e.append(hit["att_end"] if hit else np.nan)
                att_l.append(hit["att_len"] if hit else np.nan)
            c["att_start"], c["att_end"], c["att_len"] = att_s, att_e, att_l
            c["model_start"], c["model_end"] = c["start"], c["end"]
            has = c.att_start.notna()
            c.loc[has, "start"] = c.loc[has, "att_start"].astype(int)
            c.loc[has, "end"] = c.loc[has, "att_end"].astype(int)
            variants["att" if name == "shipped" else "refine+att"] = c
        out.append((mode, variants))
    return out


def att_for(contig: str, text: str, s: int, e: int, p: float) -> dict | None:
    """The tuned att rule (ATT_RULE) on one call."""
    cands = ab.find_candidates(contig, text, s, e)
    if not cands:
        return None
    item = pd.DataFrame([{"item_id": 0, "searchable": True, "s": s, "e": e, "contig_len": len(contig),
                          "core_lo": (s + e) // 2, "core_hi": (s + e) // 2, "p_match": p}])
    sel = ab.select(item, pd.DataFrame(cands).assign(item_id=0), **ATT_RULE)
    return None if sel.empty else sel.iloc[0].to_dict()


def metrics(ev: pd.DataFrame, reg: pd.DataFrame, n_genomes: int) -> dict:
    v = ek.score_variant(ev, reg)
    err = ab.edge_errors(reg, ev)
    return {"calls_per_genome": len(reg) / max(1, n_genomes), "covered50": ek.covered50(ev, reg),
            "covered": v["covered_recall"], "matched": v["matched_recall"], "median_iou": v["median_iou"],
            "absorbed": v["absorbed_elements"], "fragmented": v["fragmented_elements"], **ab._cost(ev, reg),
            "edge_med_bp": float(np.median(err)) if len(err) else np.nan,
            "edge_le100": float(np.mean(err <= 100)) if len(err) else np.nan,
            "edge_le1k": float(np.mean(err <= 1000)) if len(err) else np.nan}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", required=True)
    ap.add_argument("--genomes_dir", required=True)
    ap.add_argument("--elements", required=True)
    ap.add_argument("--modes", nargs="+", default=MODES)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    dirs = sorted(d for d in Path(args.runs).iterdir() if (d / "fine_windows.tsv").exists())
    jobs = [{"genome": d.name, "dir": str(d), "fasta": str(Path(args.genomes_dir) / f"{d.name}.fna"),
             "modes": args.modes} for d in dirs]
    by = {}
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        for res in ex.map(genome_calls, jobs):
            for mode, variants in res:
                for name, c in variants.items():
                    by.setdefault((mode, name), []).append(c)
    truth = pd.read_csv(args.elements)
    ev = truth[truth.accession.isin([d.name for d in dirs])]
    rows, keep = [], []
    for (mode, name), parts in by.items():
        reg = ek.concat_nonempty(parts)
        row = {"mode": mode, "variant": name, **metrics(ev, reg, len(dirs))}
        if "att_start" in reg:
            row["att_frac_calls"] = float(reg.att_start.notna().mean())
        rows.append(row)
        if name == "refine+att":
            keep.append(reg.assign(mode=mode))
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(f"{args.out}.csv", index=False)
    pd.concat(keep, ignore_index=True).to_csv(f"{args.out}.calls.tsv", sep="\t", index=False)
    print(pd.DataFrame(rows).round(3).to_string(index=False))


if __name__ == "__main__":
    main()
