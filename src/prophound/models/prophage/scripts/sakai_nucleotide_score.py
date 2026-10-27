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
"""Scores the Sakai nucleotide (megablast) test as preregistered (docs/sakai-candidates.md).

Reads the per-region BLAST XML fetched by scripts/ncbi_blast_jobs.py. Per region: the fraction
of its bases covered by the union of its hits (E <= 1e-10), and its longest alignment. A region
passes the virus search at 10% or more. The non-Escherichia search is reported alongside.

  python scripts/sakai_nucleotide_score.py --dir analysis/sakai_novel
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from Bio.Blast import NCBIXML

LOCAL = {"default_novel": "default", "q_exploratory_only": "exploratory_only", "q_positive_control": "positive_control",
         "q_negative_control": "negative_control"}
SETS = {"v_default_novel": "default", "v_q_exploratory_only": "exploratory_only", "v_q_positive_control": "positive_control",
        "v_q_negative_control": "negative_control", "ne_default_novel": "default (non-Escherichia search)"}
E_MAX = 1e-10
MIN_FRAC = 0.10


def union(iv: list[tuple[int, int]]) -> int:
    tot, cur = 0, None
    for a, b in sorted(iv):
        if cur and a <= cur[1]:
            cur = (cur[0], max(cur[1], b))
        else:
            if cur:
                tot += cur[1] - cur[0] + 1
            cur = (a, b)
    return tot + (cur[1] - cur[0] + 1 if cur else 0)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__); ap.add_argument("--dir", required=True, type=Path)
    ap.add_argument("--local", type=Path, help="a directory of multi-query XML from a local search (blast_<set>.xml)")
    a = ap.parse_args()
    rows = []
    if a.local:
        files = [(f, LOCAL.get(f.stem[len("blast_"):], f.stem), None) for f in sorted(a.local.glob("blast_*.xml"))]
    else:
        files = [(f, SETS.get(f.stem[len("blast_"):].split("__", 1)[0], f.stem), f.stem.split("__", 1)[1]) for f in sorted(a.dir.glob("blast_*__*.xml"))]
    for f, setname, region0 in files:
        with open(f) as fh:
            for rec in NCBIXML.parse(fh):
                region = region0 or rec.query.split("|")[0].split()[0]
                L = rec.query_length
                iv, best = [], None
                for al in rec.alignments:
                    for h in al.hsps:
                        if h.expect > E_MAX:
                            continue
                        iv.append((min(h.query_start, h.query_end), max(h.query_start, h.query_end)))
                        if best is None or h.align_length > best[0]:
                            best = (h.align_length, 100 * h.identities / h.align_length, al.hit_def[:90])
                cov = union(iv) / L if L else 0
                rows.append(dict(search=setname, region=region, kb=round(L / 1000, 1), nt_cov=round(cov, 3),
                                 pass_10pct=cov >= MIN_FRAC, hsps=len(iv), longest_aln=best[0] if best else 0,
                                 longest_pident=round(best[1], 1) if best else None, longest_hit=best[2] if best else ""))
    t = pd.DataFrame(rows)
    tag = "nucleotide_local" if a.local else "nucleotide"
    t.to_csv(a.dir / f"{tag}.regions.tsv", sep="\t", index=False)
    s = t.groupby("search", sort=False).agg(regions=("region", "size"), passed=("pass_10pct", "sum"),
                                            median_cov=("nt_cov", "median"), max_cov=("nt_cov", "max"))
    s.to_csv(a.dir / f"{tag}.summary.tsv", sep="\t")
    with pd.option_context("display.width", 230, "display.max_colwidth", 90):
        print(s.to_string()); print()
        print(t[t.search.str.startswith(("default", "exploratory"))].sort_values("nt_cov", ascending=False).head(15).to_string(index=False))


if __name__ == "__main__":
    main()
