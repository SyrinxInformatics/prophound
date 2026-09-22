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


"""Probability-saturation diagnostics for scored genomes (docs/scientific-validity.md item 1a).

Separates two failure modes that both look like "a missed region" in the
region tables but are mechanically different and need different fixes:

* **near-threshold** -- windows inside the true region score high but below the
  operating threshold (0.67-0.99 for the four persistent misses on the
  size-stratified sample). Recoverable by lowering the threshold; Exploratory
  mode recovers all but one.
* **saturated** -- windows score at numerical zero (1e-19 to 0.0, i.e. logits
  around -100), unrecoverable at *any* threshold. The only instance in the
  current evidence base is EDL933 (`GCF_000006665.1`) 892,845-930,014.

`profile` prints the per-window trace around a region so the two can be told
apart by eye; `survey` reports each genome's saturated fraction next to its
geNomad prophage density, which is how the "dense genomes saturate" hypothesis
was refuted (corr -0.09 over 42 held-out genomes).

Examples (from the prophage model dir):
  python scripts/saturation_diagnostics.py profile --windows <dir>/windows.tsv \
      --contig NC_002655.2 --start 892845 --end 930014 --flank 30000
  python scripts/saturation_diagnostics.py survey --predict_dir <dir> --summary summary.csv
"""

from __future__ import annotations

import argparse
import glob
import os

import pandas as pd

SATURATED = 1e-10  # below this a window cannot be recovered by any threshold


def cmd_profile(args: argparse.Namespace) -> None:
    w = pd.read_csv(args.windows, sep="\t")
    w = w[w.contig_id == args.contig]
    inside = w[(w.end > args.start) & (w.start < args.end)]
    near = w[(w.end > args.start - args.flank) & (w.start < args.end + args.flank)]
    if inside.empty:
        raise SystemExit(f"no windows overlap {args.contig}:{args.start}-{args.end}")
    sat = (inside.prob < SATURATED).mean()
    verdict = ("SATURATED -- unrecoverable at any threshold" if sat > 0.5
               else "near-threshold -- recoverable by lowering the threshold" if inside.prob.max() > 0.5
               else "weak signal")
    print(f"{args.contig}:{args.start}-{args.end}  n={len(inside)} windows")
    print(f"  max={inside.prob.max():.3e}  median={inside.prob.median():.3e}  "
          f"saturated_frac={sat:.2f}  ->  {verdict}")
    if args.trace:
        for r in near.itertuples():
            mark = "*" if r.end > args.start and r.start < args.end else " "
            print(f"  {mark} {r.start:>9} {r.end:>9} {r.prob:.6e}")


def cmd_survey(args: argparse.Namespace) -> None:
    summ = pd.read_csv(args.summary).set_index("accession")
    rows = []
    for f in sorted(glob.glob(os.path.join(args.predict_dir, "*", "windows.tsv"))):
        g = os.path.basename(os.path.dirname(f))
        w = pd.read_csv(f, sep="\t")
        row = dict(genome=g[:30], saturated=(w.prob < SATURATED).mean(), n_windows=len(w))
        if g in summ.index:
            r = summ.loc[g]
            row.update(density=r.positive_bp / r.total_bp, n_segments=int(r.positive_segments))
        rows.append(row)
    d = pd.DataFrame(rows).sort_values("saturated", ascending=False)
    pd.set_option("display.width", 200)
    print(d.head(args.top).round(4).to_string(index=False))
    print(f"\ngenomes with any saturated window: {(d.saturated > 0).sum()}/{len(d)}")
    if "density" in d and d.density.notna().sum() > 2:
        print(f"corr(saturated, density) = {d.saturated.corr(d.density):.3f}  "
              f"corr(saturated, n_segments) = {d.saturated.corr(d.n_segments):.3f}  (n={d.density.notna().sum()})")
        print("A near-zero correlation refutes 'prophage-dense genomes saturate'.")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("profile")
    p.add_argument("--windows", required=True)
    p.add_argument("--contig", required=True)
    p.add_argument("--start", type=int, required=True)
    p.add_argument("--end", type=int, required=True)
    p.add_argument("--flank", type=int, default=30000)
    p.add_argument("--trace", action="store_true", help="print every window in the flanked range")
    p.set_defaults(fn=cmd_profile)
    p = sub.add_parser("survey")
    p.add_argument("--predict_dir", required=True, help="dir of <genome>/windows.tsv")
    p.add_argument("--summary", required=True, help="summary.csv (for geNomad prophage density)")
    p.add_argument("--top", type=int, default=12)
    p.set_defaults(fn=cmd_survey)
    args = ap.parse_args()
    args.fn(args)


if __name__ == "__main__":
    main()
