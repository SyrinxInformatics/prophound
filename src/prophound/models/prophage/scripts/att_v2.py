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
"""Att search v2, DEPhT-style: short seeds, ungapped extension with mismatches, decoy-set
threshold, model evidence instead of integrase proximity (docs/att-boundaries.md "Att v2").

The shipped search (src/att.py, v1) accepts only exact repeats >=16 bp strong enough to
stand on their own statistically, and finds pairs at ~15% of known elements. DEPhT finds
more by accepting short, imperfect cores (BLASTN word 5) and judging candidates with other
evidence. Here: 10-mer seeds between the two accepted zones (10 kb in / 5 kb out of the
coarse edges), each diagonal's seeds extended without gaps under +1 / -2 scoring (a
maximal-scoring segment), and the segment scored by E = m * n * exp(-lambda * S), lambda
from the genome's base composition. The threshold is set by the decoys, so v2 is compared
with v1 at equal false-att rates. Optional model evidence: the pair must enclose the
call's high-scoring core.

Input: the items tables `att_boundaries.py run` wrote (calls and decoys, coarse
coordinates, core span). Subcommands:
  run    every item's best v2 candidate -> <out>.v2.tsv
  score  v1 vs v2: att rate on calls and decoys, recall, IoU, bp, and edge errors against
         truth, at decoy-matched thresholds -> <out>.v2score.csv

Example:
  python scripts/att_v2.py run --items analysis/att/tune42.items.tsv \\
      --genomes_dir ../../../../arc_data/raw/genomes --out analysis/att/tune42
  python scripts/att_v2.py score --items analysis/att/tune42.items.tsv --v1 analysis/att/tune42 \\
      --v2 analysis/att/tune42.v2.tsv --elements scripts/eval_sizestrat_elements.csv --out analysis/att/tune42
"""

from __future__ import annotations

import argparse
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import att_boundaries as ab
import evalkit as ek
import numpy as np
import pandas as pd

SEED_K = 10
MATCH, MISMATCH = 1, -2
FLANK = 60            # bp examined beyond a diagonal's outermost seeds
GROUP_GAP = 100       # seeds on one diagonal closer than this are one candidate
MIN_LEN, MAX_LEN = 12, 300
MIN_ELEMENT = 2000
ACCEPT_IN, ACCEPT_OUT = 10000, 5000


def lam(p: float) -> float:
    """Karlin-Altschul lambda for +MATCH/MISMATCH scoring at match probability p."""
    lo, hi = 1e-6, 10.0
    f = lambda x: p * math.exp(MATCH * x) + (1 - p) * math.exp(MISMATCH * x) - 1   # noqa: E731
    for _ in range(100):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if f(mid) < 0 else (lo, mid)
    return (lo + hi) / 2


def best_segment(sc: np.ndarray) -> tuple[int, int, float]:
    """Maximal-scoring contiguous segment, by prefix sums: start, end (exclusive), score."""
    c = np.concatenate([[0.0], np.cumsum(sc)])
    run_min = np.minimum.accumulate(c)
    gain = c[1:] - run_min[:-1]
    j = int(np.argmax(gain))
    i = int(np.argmin(c[:j + 1]))
    return i, j + 1, float(gain[j])


def candidates(left: str, right: str) -> list:
    """Ungapped near-identical segments between two strings: (li, ri, length, score, identity)."""
    a, b = ab.seed_table(ab.encode(left), SEED_K), ab.seed_table(ab.encode(right), SEED_K)
    hits = a.merge(b, on="code", suffixes=("_l", "_r"))
    if hits.empty:
        return []
    hits["diag"] = hits.pos_r - hits.pos_l
    L = np.frombuffer(left.encode(), dtype=np.uint8)
    R = np.frombuffer(right.encode(), dtype=np.uint8)
    out = []
    for diag, g in hits.groupby("diag"):
        pos = np.sort(g.pos_l.to_numpy())
        cuts = np.flatnonzero(np.diff(pos) > GROUP_GAP)
        for grp in np.split(pos, cuts + 1):
            i0 = max(0, int(grp[0]) - FLANK, -int(diag))
            i1 = min(len(L), int(grp[-1]) + SEED_K + FLANK, len(R) - int(diag))
            if i1 - i0 < MIN_LEN:
                continue
            same = L[i0:i1] == R[i0 + diag:i1 + diag]
            s, e, score = best_segment(np.where(same, MATCH, MISMATCH))
            if e - s >= MIN_LEN:
                out.append((i0 + s, i0 + s + diag, e - s, float(score), float(same[s:e].mean())))
    return out


def item_job(job: dict) -> list:
    seqs, text, _ = ab.load_genome(Path(job["fasta"]))
    rows = []
    for it in job["items"]:
        if not it["searchable"]:
            continue
        contig = seqs[it["contig_id"]]
        s, e = int(it["s"]), int(it["e"])
        mid = (s + e) // 2
        l0, l1 = max(0, s - ACCEPT_OUT), min(mid, s + ACCEPT_IN)
        r0, r1 = max(mid, e - ACCEPT_IN), min(len(contig), e + ACCEPT_OUT)
        if l1 - l0 < SEED_K or r1 - r0 < SEED_K:
            continue
        lmb, m, n = lam(float(it["p_match"])), l1 - l0, r1 - r0
        best = None
        for li, ri, length, score, ident in candidates(contig[l0:l1], contig[r0:r1]):
            start, end = l0 + li, r0 + ri + length
            if length > MAX_LEN or end - start < MIN_ELEMENT:
                continue
            ev = m * n * math.exp(-lmb * score)
            if best is not None and ev >= best["v2_evalue"]:
                continue
            core = contig[start:start + length].upper()
            if ab.entropy_bits(core) < 1.5:
                continue
            probe = core[(length - 16) // 2:(length - 16) // 2 + 16] if length >= 16 else core
            rc = probe.translate(ab.COMP)[::-1]
            copies = text.count(probe) + (text.count(rc) if rc != probe else 0)
            if copies > 2:
                continue
            best = {"item_id": it["item_id"], "v2_start": start, "v2_end": end, "v2_len": length, "v2_score": score,
                    "v2_identity": ident, "v2_evalue": ev,
                    "v2_core_ok": bool(start <= int(it["core_lo"]) and end >= int(it["core_hi"]))}
        if best:
            rows.append(best)
    return rows


def cmd_run(a) -> None:
    items = pd.read_csv(a.items, sep="\t")
    jobs = [{"fasta": str(Path(a.genomes_dir) / f"{g}.fna"), "items": sub.to_dict("records")}
            for g, sub in items.groupby("genome_id")]
    rows = []
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        for i, r in enumerate(ex.map(item_job, jobs)):
            rows += r
            print(f"[v2] {i + 1}/{len(jobs)} genomes", end="\r", flush=True)
    pd.DataFrame(rows).to_csv(f"{a.out}.v2.tsv", sep="\t", index=False)
    print(f"\n[v2] {len(rows)} items with a candidate -> {a.out}.v2.tsv")


def evaluate(items, chosen_start, chosen_end, ev, mode):
    c = items[(items.kind == "call") & (items["mode"] == mode)]
    s = np.where(chosen_start.reindex(c.item_id).notna(), chosen_start.reindex(c.item_id), c.call_start).astype(int)
    e = np.where(chosen_end.reindex(c.item_id).notna(), chosen_end.reindex(c.item_id), c.call_end).astype(int)
    reg = pd.DataFrame({"genome_id": c.genome_id.to_numpy(), "contig_id": c.contig_id.to_numpy(), "start": s, "end": e})
    att = chosen_start.reindex(c.item_id).notna().to_numpy()
    err = ab.edge_errors(reg[att], ev)
    v = ek.score_variant(ev, reg)
    return {"covered50": ek.covered50(ev, reg), "matched": v["matched_recall"], "median_iou": v["median_iou"],
            "Mb_called": (reg.end - reg.start).sum() / 1e6, "att_calls": float(att.mean()), "att_n": int(att.sum()),
            "att_edge_n": int(len(err) // 2), "att_edge_med_bp": float(np.median(err)) if len(err) else np.nan,
            "att_edge_le100": float(np.mean(err <= 100)) if len(err) else np.nan}


def cmd_score(a) -> None:
    items = pd.read_csv(a.items, sep="\t")
    cands = pd.read_csv(f"{a.v1}.cands.tsv", sep="\t")
    v2 = pd.read_csv(a.v2, sep="\t").set_index("item_id")
    truth = pd.read_csv(a.elements)
    ev = truth[truth.accession.isin(items.genome_id.unique())]
    v1 = ab.select(items, cands, 10000, 5000, 300, 2, 1.5, None, 1e-6).set_index("item_id")
    rows = []
    for mode in ("exploratory", "default", "conservative"):
        it = items[items["mode"] == mode]
        dec = it[(it.kind == "decoy") & it.searchable]
        v1_dec = float(dec.item_id.isin(v1.index).mean())
        rows.append({"mode": mode, "variant": "v1 (shipped)", "decoy_rate": v1_dec,
                     **evaluate(items, v1.att_start, v1.att_end, ev, mode)})
        for core_rule in (False, True):
            pool = v2[v2.v2_core_ok] if core_rule else v2
            dev = pool.reindex(dec.item_id).v2_evalue.dropna()
            for target in sorted({v1_dec, *a.decoy_rates}):
                # the largest E threshold whose decoy acceptance does not exceed the target
                k = int(math.floor(target * len(dec)))
                thr = np.sort(dev.to_numpy())[k - 1] if 0 < k <= len(dev) else (np.inf if k > len(dev) else -np.inf)
                sel = pool[pool.v2_evalue <= thr]
                rows.append({"mode": mode, "variant": f"v2{' +core' if core_rule else ''} @decoy {target:.3f}",
                             "decoy_rate": float(dec.item_id.isin(sel.index).mean()), "E_threshold": thr,
                             **evaluate(items, sel.v2_start, sel.v2_end, ev, mode)})
    out = pd.DataFrame(rows)
    out.to_csv(f"{a.out}.v2score.csv", index=False)
    print(out.round(4).to_string(index=False))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    r = sp.add_parser("run")
    r.add_argument("--items", required=True)
    r.add_argument("--genomes_dir", required=True)
    r.add_argument("--workers", type=int, default=8)
    r.add_argument("--out", required=True)
    r.set_defaults(fn=cmd_run)
    s = sp.add_parser("score")
    s.add_argument("--items", required=True)
    s.add_argument("--v1", required=True, help="prefix of att_boundaries.py run (for <prefix>.cands.tsv)")
    s.add_argument("--v2", required=True)
    s.add_argument("--elements", required=True)
    s.add_argument("--decoy_rates", type=float, nargs="+", default=[0.03, 0.05])
    s.add_argument("--out", required=True)
    s.set_defaults(fn=cmd_score)
    a = ap.parse_args()
    a.fn(a)


if __name__ == "__main__":
    main()
