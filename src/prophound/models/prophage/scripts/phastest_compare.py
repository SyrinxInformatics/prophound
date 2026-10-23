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
"""prophound against PHASTEST (and optionally geNomad) on one genome with author-curated
prophages. docs/final-benchmark.md "PHASTEST comparisons".

Inputs: prophound's shipped `regions.detailed.tsv` per mode (0-based half-open), PHASTEST's
API result JSON (its region table, 1-based inclusive), an optional geNomad provirus TSV
(1-based inclusive), and a truth TSV with contig, start, end (1-based inclusive), name, class.
Per tool: calls, covered50 and covered (the literature pilot's rules), covered50 per class,
median best IoU, bp called and bp off every curated element, and the edge error of calls on
exactly one element. PHASTEST is scored on all its regions and on intact + questionable
regions only. Also writes each element's coverage per tool.

Example:
  python /abs/path/scripts/phastest_compare.py --name sl1344 \\
      --truth /abs/path/analysis/phastest_cmp/sl1344_truth.tsv \\
      --prophound default=/abs/path/arc_data/phastest_cmp/sl1344/default ... \\
      --phastest /abs/path/arc_data/phastest_cmp/inputs/sl1344_phastest_ZZ_ce699ad464.json \\
      --out /abs/path/analysis/phastest_cmp/sl1344
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd


def phastest_regions(path: str, contig: str) -> pd.DataFrame:
    """PHASTEST's region table as 0-based half-open calls, with completeness and score."""
    rows = []
    for line in json.loads(Path(path).read_text())["summary"].splitlines():
        m = re.match(r"\s*(\d+)\s+(\S+)\s+([a-z]+)\((\d+)\)\s+\S+\s+(\d+)-(\d+)\s", line)
        if m:
            rows.append({"contig_id": contig, "start": int(m[5]) - 1, "end": int(m[6]),
                         "completeness": m[3], "score": int(m[4])})
    if not rows:
        raise SystemExit(f"no PHASTEST regions parsed from {path}")
    return pd.DataFrame(rows)


def score(calls: pd.DataFrame, t: pd.DataFrame) -> pd.DataFrame:
    """Per element: share covered by the union of calls, and the best single-call IoU."""
    rows = []
    for r in t.itertuples():
        x = calls[calls.contig_id == r.contig]
        L = r.e0 - r.s0
        cs, ce = x.start.to_numpy(np.int64), x.end.to_numpy(np.int64)
        ov = np.clip(np.minimum(ce, r.e0) - np.maximum(cs, r.s0), 0, None)
        if len(x):
            pos = np.zeros(L, bool)            # union coverage, overlapping calls counted once
            for a, b in zip(cs, ce):
                pos[max(0, a - r.s0):max(0, min(L, b - r.s0))] = True
            cov, iou = pos.mean(), float((ov / ((ce - cs) + L - ov)).max())
        else:
            cov, iou = 0.0, 0.0
        rows.append({"name": r.name, "class": r.cls, "len_kb": round(L / 1e3, 1), "cov_frac": round(float(cov), 3),
                     "best_iou": round(iou, 3)})
    return pd.DataFrame(rows)


def edge_errors(calls: pd.DataFrame, t: pd.DataFrame) -> np.ndarray:
    out = []
    for c in calls.itertuples():
        x = t[(t.contig == c.contig_id) & (t.e0 > c.start) & (t.s0 < c.end)]
        if len(x) == 1:
            out += [abs(int(c.start) - int(x.s0.iloc[0])), abs(int(c.end) - int(x.e0.iloc[0]))]
    return np.asarray(out)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--name", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--prophound", action="append", required=True, metavar="MODE=OUTDIR")
    ap.add_argument("--phastest", required=True)
    ap.add_argument("--genomad", default=None, help="geNomad *_provirus.tsv")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    t = pd.read_csv(a.truth, sep="\t").rename(columns={"class": "cls"})
    t["s0"], t["e0"] = t.start - 1, t.end
    contig = t.contig.iloc[0]

    tools = {}
    for spec in a.prophound:
        mode, _, d = spec.partition("=")
        tools[f"prophound {mode}"] = pd.read_csv(Path(d) / "regions.detailed.tsv", sep="\t")
    if a.genomad:
        g = pd.read_csv(a.genomad, sep="\t")
        tools["geNomad"] = pd.DataFrame({"contig_id": g.source_seq, "start": g.start - 1, "end": g.end})
    ph = phastest_regions(a.phastest, contig)
    tools["PHASTEST (all)"] = ph
    tools["PHASTEST (intact+questionable)"] = ph[ph.completeness != "incomplete"]

    rows, per = [], []
    for name, c in tools.items():
        c = c[c.contig_id == contig]
        es = score(c, t)
        per.append(es.assign(tool=name))
        on = sum(int(np.clip(np.minimum(x.end, t.e0) - np.maximum(x.start, t.s0), 0, None).sum())
                 for x in c.itertuples())
        bp = int((c.end - c.start).sum())
        err = edge_errors(c, t)
        row = {"genome": a.name, "tool": name, "calls": len(c),
               "covered50": f"{int((es.cov_frac >= 0.5).sum())}/{len(es)}",
               "covered": f"{int((es.cov_frac > 0).sum())}/{len(es)}"}
        for cl, sub in es.groupby("class"):
            row[f"covered50_{cl}"] = f"{int((sub.cov_frac >= 0.5).sum())}/{len(sub)}"
        row.update({"median_best_iou": float(es.best_iou.median()), "Mb_called": round(bp / 1e6, 3),
                    "Mb_off_curated": round((bp - on) / 1e6, 3),
                    "edge_median_kb": round(float(np.median(err)) / 1e3, 2) if len(err) else np.nan,
                    "edges_within_1kb": f"{int((err <= 1000).sum())}/{len(err)}"})
        rows.append(row)
    summary = pd.DataFrame(rows)
    elements = pd.concat(per)
    wide = elements.pivot_table(index=["name", "class", "len_kb"], columns="tool", values="cov_frac", sort=False)
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(f"{out}.csv", index=False)
    elements.to_csv(f"{out}.elements.csv", index=False)
    ph.to_csv(f"{out}.phastest_regions.csv", index=False)
    pd.set_option("display.width", 250)
    print(summary.to_string(index=False))
    print(wide[list(tools)].to_string())
    print(f"[done] wrote {out}.csv, {out}.elements.csv, {out}.phastest_regions.csv")


if __name__ == "__main__":
    main()
