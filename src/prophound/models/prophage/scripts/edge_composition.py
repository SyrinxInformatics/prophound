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
"""Composition change-points for estimated edges (docs/att-boundaries.md "Composition edges").

Prophages usually differ from their host in base composition. As in Alien_Hunter's
change-point step and CheckV's adjacent-window comparison, an estimated edge is moved to
where the sequence turns from host-like to element-like:

  - profiles: tetranucleotide frequencies of the element (the middle 60% of the call)
    and of the local host (FLANK bp beyond the edge's range, on the outer side)
  - windows: WIN bp every STEP bp across the edge's range, each scored by the
    log-likelihood ratio element vs host (multinomial on tetranucleotides, pseudocounts)
  - change-point: the position maximising (sum of window scores on the element side) -
    (sum on the host side), i.e. the best single host -> element switch in the range
  - gate: move the edge only if the profiles differ (Jensen-Shannon divergence >= jsd) and
    the change-point's gain over the reported edge is >= gain

Only `model` edges that are not on a contig end move; att edges are unchanged. Input: the
refine+att calls from boundary_eval.py (<sample>_pipeline.calls.tsv, all modes). Output:
metrics per gate setting, and the edges table.

Example:
  python scripts/edge_composition.py --calls analysis/att/tune42_pipeline.calls.tsv \\
      --genomes_dir ../../../../arc_data/raw/genomes --elements scripts/eval_sizestrat_elements.csv \\
      --out analysis/att/tune42_comp
"""

from __future__ import annotations

import argparse
import itertools
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import att_boundaries as ab
import boundary_eval as be
import numpy as np
import pandas as pd

K = 4
WIN, STEP = 1000, 100
FLANK = 20000
PSEUDO = 0.5


def kmer_codes(seq: str) -> np.ndarray:
    a = ab.encode(seq)
    n = len(a) - K + 1
    if n <= 0:
        return np.zeros(0, np.int64)
    bad = np.concatenate([[0], np.cumsum(a == 255)])
    ok = (bad[K:K + n] - bad[:n]) == 0
    code = np.zeros(n, np.int64)
    a64 = np.where(a == 255, 0, a).astype(np.int64)
    for i in range(K):
        code = code * 4 + a64[i:i + n]
    return np.where(ok, code, -1)


def profile(codes: np.ndarray) -> np.ndarray:
    c = np.bincount(codes[codes >= 0], minlength=4 ** K).astype(float) + PSEUDO
    return c / c.sum()


def jsd(p: np.ndarray, q: np.ndarray) -> float:
    m = (p + q) / 2
    return float(0.5 * (p * np.log2(p / m)).sum() + 0.5 * (q * np.log2(q / m)).sum())


def change_point(contig: str, lo: int, hi: int, elem: np.ndarray, host: np.ndarray, side: str, edge: int):
    """Best host->element switch in [lo, hi]: (position, gain over `edge`)."""
    a, b = max(0, lo - WIN // 2), min(len(contig), hi + WIN // 2)
    codes = kmer_codes(contig[a:b])
    if len(codes) < WIN:
        return edge, 0.0
    w = np.log(elem) - np.log(host)
    per = np.where(codes >= 0, w[np.maximum(codes, 0)], 0.0)
    c = np.concatenate([[0.0], np.cumsum(per)])
    starts = np.arange(0, len(per) - WIN + 1, STEP)
    llr = c[starts + WIN] - c[starts]                   # each window's element-vs-host score
    centres = a + starts + WIN // 2
    keep = (centres >= lo) & (centres <= hi)
    llr, centres = llr[keep], centres[keep]
    if len(llr) < 3:
        return edge, 0.0
    # element on the right of a start edge, on the left of an end edge
    sign = 1.0 if side == "start" else -1.0
    cs = np.concatenate([[0.0], np.cumsum(sign * llr)])
    total = cs[-1]
    split = total - 2 * cs                               # (sum right) - (sum left) for each split
    j = int(np.argmax(split))
    pos = int(centres[j]) if j < len(centres) else int(centres[-1])
    j0 = int(np.searchsorted(centres, edge))
    return pos, float(split[j] - split[min(j0, len(split) - 1)])


def genome_job(job: dict) -> list:
    seqs, _, _ = ab.load_genome(Path(job["fasta"]))
    out = []
    for r in job["calls"]:
        contig = seqs[r["contig_id"]]
        s, e = int(r["start"]), int(r["end"])
        L = e - s
        rec = {"key": r["key"]}
        if r["boundary"] != "model" or L < 3000:
            out.append(rec)
            continue
        elem = profile(kmer_codes(contig[s + int(0.2 * L):e - int(0.2 * L)]))
        for side in ("start", "end"):
            if r["contig_edge"] in (side, "both"):
                continue
            lo, hi = int(r[f"{side}_lo"]), int(r[f"{side}_hi"])
            host_seq = contig[max(0, lo - FLANK):lo] if side == "start" else contig[hi:hi + FLANK]
            if len(host_seq) < 5000:
                continue
            host = profile(kmer_codes(host_seq))
            pos, gain = change_point(contig, lo, hi, elem, host, side, s if side == "start" else e)
            rec.update({f"{side}_cp": pos, f"{side}_gain": gain, f"{side}_jsd": jsd(elem, host)})
        out.append(rec)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--calls", required=True)
    ap.add_argument("--genomes_dir", required=True)
    ap.add_argument("--elements", required=True)
    ap.add_argument("--jsd", type=float, nargs="+", default=[0.0, 0.02, 0.04, 0.06])
    ap.add_argument("--gain", type=float, nargs="+", default=[0.0, 5.0, 20.0, 50.0])
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    calls = pd.read_csv(a.calls, sep="\t")
    calls["key"] = np.arange(len(calls))
    jobs = [{"fasta": str(Path(a.genomes_dir) / f"{g}.fna"), "calls": sub.to_dict("records")}
            for g, sub in calls.groupby("genome_id")]
    rows = []
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        for r in ex.map(genome_job, jobs):
            rows += r
    cp = pd.DataFrame(rows).set_index("key").reindex(calls.key)
    for col in ("start_cp", "end_cp", "start_gain", "end_gain", "start_jsd", "end_jsd"):
        calls[col] = cp[col].to_numpy() if col in cp else np.nan
    calls.to_csv(f"{a.out}.cp.tsv", sep="\t", index=False)
    truth = pd.read_csv(a.elements)
    ev = truth[truth.accession.isin(calls.genome_id.unique())]
    res = []
    n_g = calls.genome_id.nunique()
    for mode, cm in calls.groupby("mode", sort=False):
        res.append({"mode": mode, "variant": "shipped", **be.metrics(ev, cm, n_g), "moved_edges": 0})
        for t_jsd, t_gain in itertools.product(a.jsd, a.gain):
            v = cm.copy()
            moved = 0
            for side in ("start", "end"):
                ok = v[f"{side}_cp"].notna() & (v[f"{side}_jsd"] >= t_jsd) & (v[f"{side}_gain"] >= t_gain)
                v.loc[ok, side] = v.loc[ok, f"{side}_cp"].astype(int)
                moved += int(ok.sum())
            v = v[v.end > v.start]
            res.append({"mode": mode, "variant": f"jsd>={t_jsd:g} gain>={t_gain:g}", **be.metrics(ev, v, n_g),
                        "moved_edges": moved})
    out = pd.DataFrame(res)
    out.to_csv(f"{a.out}.csv", index=False)
    print(out[["mode", "variant", "covered50", "covered", "matched", "median_iou", "Mb_called", "Mb_off",
               "edge_med_bp", "edge_le1k", "moved_edges"]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
