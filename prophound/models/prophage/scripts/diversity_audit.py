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


"""Audit where one-size-fits-all assumptions bite, across the diversity axes.

Operates on already-scored output
(`<genome>/windows.tsv` + `regions.tsv`), so it costs no inference and can be
re-run against any sweep directory.

Four reports, one per axis the design doc names:

  `background`  per-genome window-score distribution and the fraction clearing a
                global threshold. Establishes whether a single threshold is
                defensible (measured 2026-09-22: it is not -- 47x spread).
  `size`        per-stratum recall, IoU and predicted length, i.e. whether the
                pipeline treats large and small elements equally.
  `context`     how many calls contain >=2 truth elements (the splitting
                opportunity) and how many truth elements are split across >=2
                calls (the merging opportunity). Sizes both before either is
                built.
  `merge`       gap-budget sanity: for each pair of adjacent calls, the gap
                between them relative to their own lengths. Flags where the
                absolute `--max-gap-bp` is large compared with the elements it
                would join -- the defect that makes a 50 kb budget wrong for a
                2 kb element.

Example:
  python scripts/diversity_audit.py background --predict_dir <dir> --summary summary.csv
  python scripts/diversity_audit.py context --predict_dir <dir> --elements scripts/eval_sizestrat_elements.csv
"""

from __future__ import annotations

import argparse
import glob
import os

import numpy as np
import pandas as pd


def load_tsv(d: str, name: str) -> pd.DataFrame:
    dfs = []
    for f in glob.glob(os.path.join(d, "*", name)):
        if name == "regions.tsv" and os.path.exists(os.path.join(os.path.dirname(f), "regions.detailed.tsv")):
            f = os.path.join(os.path.dirname(f), "regions.detailed.tsv")   # 0-based calls since 0.1.0
        try:
            dfs.append(pd.read_csv(f, sep="\t"))
        except pd.errors.EmptyDataError:
            pass
    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()


def cmd_background(args: argparse.Namespace) -> None:
    rows = []
    summ = pd.read_csv(args.summary).set_index("accession") if args.summary else None
    for f in sorted(glob.glob(os.path.join(args.predict_dir, "*", "windows.tsv"))):
        g = os.path.basename(os.path.dirname(f))
        w = pd.read_csv(f, sep="\t")
        p = w.prob.values
        row = dict(genome=g[:28], n_windows=len(p), median=np.median(p),
                   q90=np.quantile(p, .90), q99=np.quantile(p, .99),
                   frac_over=(p >= args.threshold).mean())
        if summ is not None and g in summ.index:
            r = summ.loc[g]
            row["genomad_frac_bp"] = r.positive_bp / r.total_bp
        rows.append(row)
    d = pd.DataFrame(rows).sort_values("frac_over", ascending=False)
    pd.set_option("display.width", 200)
    print(d.head(args.top).round(5).to_string(index=False))
    lo, hi = d["median"].min(), d["median"].max()
    print(f"\ngenomes: {len(d)}")
    print(f"background (median window score) spread: {lo:.2e} -> {hi:.2e} = {hi/max(lo,1e-12):.0f}x")
    print(f"fraction clearing {args.threshold}: min {d.frac_over.min():.4f}, "
          f"median {d.frac_over.median():.4f}, max {d.frac_over.max():.4f} "
          f"({d.frac_over.max()/max(d.frac_over.min(),1e-12):.0f}x)")
    if "genomad_frac_bp" in d:
        r = d.frac_over.corr(d.genomad_frac_bp)
        print(f"correlation(frac_over, geNomad prophage bp fraction) = {r:.3f}")
        print("Read this against the background spread above. A HIGH correlation (measured 0.888 on "
              "the 42-genome sample) means the global threshold's upper tail is mostly tracking real "
              "prophage content, so a genome with many above-threshold windows is not thereby "
              "over-calling -- do not argue for adaptive thresholds from frac_over alone. The case "
              "rests on the BACKGROUND spread (the lower tail), which a fixed threshold ignores, and "
              "on whether real elements sit just under the threshold.")


def _overlaps(reg: pd.DataFrame, acc: str, contig: str, s: int, e: int) -> pd.DataFrame:
    g = reg[(reg.genome_id == acc) & (reg.contig_id == contig)]
    return g[(g.end > s) & (g.start < e)]


def cmd_size(args: argparse.Namespace) -> None:
    ev = pd.read_csv(args.elements)
    reg = load_tsv(args.predict_dir, "regions.tsv")
    rows = []
    for r in ev.itertuples():
        o = _overlaps(reg, r.accession, r.contig, r.start, r.end)
        if len(o):
            inter = np.minimum(o.end, r.end) - np.maximum(o.start, r.start)
            union = np.maximum(o.end, r.end) - np.minimum(o.start, r.start)
            iou = float((inter / union).max())
            best = o.iloc[int((inter / union).values.argmax())]
            called = int(best.end - best.start)
        else:
            iou, called = 0.0, 0
        rows.append(dict(stratum=getattr(r, "stratum", "all"), kb=r.kb,
                         covered=len(o) > 0, iou=iou, called_len=called,
                         n_calls=len(o)))
    d = pd.DataFrame(rows)
    g = d.groupby("stratum").agg(n=("kb", "size"), covered_recall=("covered", "mean"),
                                 median_iou=("iou", "median"), median_true_kb=("kb", "median"),
                                 median_called_kb=("called_len", lambda x: np.median(x) / 1000),
                                 fragmented=("n_calls", lambda x: (x > 1).mean()))
    print(g.round(3).to_string())
    print("\nA flat median_iou column across strata is the goal; a falling one means the "
          "pipeline is tuned for large elements.")


def cmd_context(args: argparse.Namespace) -> None:
    ev = pd.read_csv(args.elements)
    reg = load_tsv(args.predict_dir, "regions.tsv")
    # calls containing >=2 truth elements -> splitting opportunity
    per_call = {}
    for r in ev.itertuples():
        for o in _overlaps(reg, r.accession, r.contig, r.start, r.end).itertuples():
            per_call.setdefault((o.genome_id, o.contig_id, o.start, o.end), []).append(r.kb)
    multi = {k: v for k, v in per_call.items() if len(v) >= 2}
    # truth elements split across >=2 calls -> merging opportunity
    frag = sum(1 for r in ev.itertuples()
               if len(_overlaps(reg, r.accession, r.contig, r.start, r.end)) >= 2)
    print(f"truth elements: {len(ev)}; calls overlapping >=1 truth element: {len(per_call)}")
    print(f"\nSPLIT opportunity -- calls containing >=2 truth elements: {len(multi)}")
    if multi:
        sizes = [len(v) for v in multi.values()]
        print(f"  elements absorbed this way: {sum(sizes)} ({sum(sizes)/len(ev)*100:.1f}% of truth)")
        print(f"  max elements under one call: {max(sizes)}")
        ex = sorted(multi.items(), key=lambda kv: -len(kv[1]))[:5]
        for (g, c, s, e), kbs in ex:
            print(f"    {g[:26]} {c} {s}-{e} ({(e-s)//1000} kb) covers {len(kbs)}: "
                  f"{[round(k,1) for k in sorted(kbs)]} kb")
    print(f"\nMERGE opportunity -- truth elements hit by >=2 calls: {frag} "
          f"({frag/len(ev)*100:.1f}%)")
    print("\nIf the split count is near zero the splitting work in Part 3 is not worth building; "
          "if it is material, the current greedy 1-1 metric is scoring those as single hits.")


def cmd_merge(args: argparse.Namespace) -> None:
    reg = load_tsv(args.predict_dir, "regions.tsv")
    rows = []
    for (g, c), sub in reg.groupby(["genome_id", "contig_id"]):
        sub = sub.sort_values("start")
        prev = None
        for r in sub.itertuples():
            if prev is not None:
                gap = r.start - prev.end
                if gap >= 0:
                    smaller = min(prev.end - prev.start, r.end - r.start)
                    rows.append(dict(gap_bp=gap, smaller_call_bp=smaller,
                                     ratio=gap / max(1, smaller)))
            prev = r
    d = pd.DataFrame(rows)
    if d.empty:
        print("no adjacent call pairs found")
        return
    print(f"adjacent call pairs: {len(d)}")
    print(d[["gap_bp", "smaller_call_bp", "ratio"]].describe(
        percentiles=[.25, .5, .75, .9]).round(2).to_string())
    for budget in (10000, 25000, 50000):
        would = (d.gap_bp <= budget)
        risky = would & (d.ratio > 1.0)
        print(f"  budget {budget:>6} bp: would merge {would.sum():>4} pairs, of which "
              f"{risky.sum():>3} span a gap LARGER than the smaller call itself")
    print("\nA gap wider than the element it is being joined to is the signature of an absolute "
          "budget applied to a small element -- the defect a proportional budget fixes.")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name, fn, needs_elements in (("background", cmd_background, False), ("size", cmd_size, True),
                                     ("context", cmd_context, True), ("merge", cmd_merge, False)):
        p = sub.add_parser(name)
        p.add_argument("--predict_dir", required=True, help="dir of <genome>/{windows,regions}.tsv")
        if needs_elements:
            p.add_argument("--elements", required=True)
        p.set_defaults(fn=fn)
    sub.choices["background"].add_argument("--summary", help="summary.csv, for geNomad prophage density")
    sub.choices["background"].add_argument("--threshold", type=float, default=0.995)
    sub.choices["background"].add_argument("--top", type=int, default=12)
    args = ap.parse_args()
    args.fn(args)


if __name__ == "__main__":
    main()
