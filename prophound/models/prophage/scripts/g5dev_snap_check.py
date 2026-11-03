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
"""Confirms the shipped edge snap (predict_genome.snap_contig_edges at
DEFAULT_EDGE_SNAP_FRAC x the run's window length) against what was measured in
development (a flat 5 kb, g5dev_eval.edge_snap), on the corpus shreds and the
phage-free genomes. Freeze check, 2026-10-05; at 10 kb windows the two must agree.

Run from prophound/models/prophage:
  python scripts/g5dev_snap_check.py --out analysis/g5dev/snap_check.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import evalkit as ek
import numpy as np
import pandas as pd
from evalkit import DEFAULT_EDGE_SNAP_FRAC, cov_frac, snap_contig_edges
from g5dev_compare import best_iou, off_element_bp
from g5dev_eval import edge_snap, load, to_genome

MAGS = ["mag_s0", "mag_s1", "mag_s2"]
CONDS = ["complete", *MAGS, "draft_s0"]
MODES = ["default", "exploratory", "conservative"]


def shipped_snap(base: pd.DataFrame, lengths: dict, window_len: dict) -> pd.DataFrame:
    """predict_genome's snap, genome by genome, at each run's own window length."""
    parts = []
    for g, sub in base.groupby("genome_id"):
        L = {c: n for (gg, c), n in lengths.items() if gg == g}
        r, _ = snap_contig_edges(sub.assign(contig_id=sub.contig_id.astype(str)), L,
                                 int(DEFAULT_EDGE_SNAP_FRAC * window_len[g]))
        parts.append(r)
    return pd.concat(parts, ignore_index=True) if parts else base


def score(T, calls, free_mb, n_genomes):
    h = cov_frac(T, calls) >= 0.5
    on = np.array([((T.assembly == c.assembly) & (T.sequence == c.sequence) & (T.end > c.start)
                    & (T.start < c.end)).any() for c in calls.itertuples()], bool)
    iou = best_iou(T, calls)
    return h, {"calls_per_genome": len(calls) / n_genomes, "false_per_mb": float((~on).sum()) / free_mb,
               "bp_called": int((calls.end - calls.start).sum()), "bp_off_element": int(off_element_bp(calls, T).sum()),
               "iou": iou[h]}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run_root", default="analysis/g5dev/prophound")
    ap.add_argument("--lengths", default="analysis/g5dev/contig_lengths.tsv")
    ap.add_argument("--set", action="append", default=None, metavar="NAME=ELEMENTS_CSV")
    ap.add_argument("--phagefree", default="analysis/phagefree")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    sets = a.set or ["tune=scripts/eval_sizestrat_elements.csv", "check=scripts/eval_confirm18_elements.csv"]
    Lt = pd.read_csv(a.lengths, sep="\t", dtype={"contig_id": str})
    variants = ("shipped", "flat5kb", "frac")
    res = {}
    for sname, path in (x.split("=", 1) for x in sets):
        T = pd.read_csv(path).rename(columns={"accession": "assembly", "contig": "sequence"})
        T = T[["assembly", "sequence", "start", "end"]].reset_index(drop=True)
        G = sorted(T.assembly.unique())
        bp = Lt[(Lt.condition == "complete") & Lt.genome.isin(G)].groupby("genome").length.sum()
        el = (T.end - T.start).groupby(T.assembly).sum()
        free_mb = float((bp - el.reindex(bp.index).fillna(0)).sum()) / 1e6
        cells, hits = {}, {}
        for cond in CONDS:
            ens, seeds = load(Path(a.run_root), cond, G)
            lc = Lt[(Lt.condition == cond) & Lt.genome.isin(G)]
            lengths = {(g, c): int(n) for g, c, n in zip(lc.genome, lc.contig_id, lc.length)}
            wl = {g: v[1] for g, v in ens.items()}
            for mode in MODES:
                base = ek.shipped_calls(ens, seeds, mode)
                regs = {"shipped": base, "flat5kb": edge_snap(base, lengths, 5000),
                        "frac": shipped_snap(base, lengths, wl)}
                for v in variants:
                    hits[(cond, mode, v)], cells[(cond, mode, v)] = score(T, to_genome(regs[v]), free_mb, len(G))
        out = res[sname] = {}
        for mode in MODES:
            out[mode] = {}
            for v in variants:
                mag = np.concatenate([hits[(c, mode, v)] for c in MAGS])
                out[mode][v] = {
                    "complete": float(hits[("complete", mode, v)].mean()), "mag": float(mag.mean()),
                    "draft": float(hits[("draft_s0", mode, v)].mean()),
                    "mag_calls_per_genome": float(np.mean([cells[(c, mode, v)]["calls_per_genome"] for c in MAGS])),
                    "mag_false_per_mb": float(np.mean([cells[(c, mode, v)]["false_per_mb"] for c in MAGS])),
                    "complete_false_per_mb": cells[("complete", mode, v)]["false_per_mb"],
                    "mag_bp_called": int(sum(cells[(c, mode, v)]["bp_called"] for c in MAGS)),
                    "mag_bp_off_element": int(sum(cells[(c, mode, v)]["bp_off_element"] for c in MAGS)),
                    "complete_bp_called": cells[("complete", mode, v)]["bp_called"],
                    "mag_median_iou": float(np.median(np.concatenate([cells[(c, mode, v)]["iou"] for c in MAGS]))),
                    "complete_median_iou": float(np.median(cells[("complete", mode, v)]["iou"]))}
    # phage-free complete genomes: contig length = furthest window end (within one stride)
    root = Path(a.phagefree)
    ens = ek.load_tracks(str(root / "ens"))
    seeds = [ek.load_tracks(str(root / r)) for r in ("r1", "r2", "r3")]
    lengths = {(g, str(c)): int(e) for g, (w, _) in ens.items() for c, e in w.groupby("contig_id").end.max().items()}
    wl = {g: v[1] for g, v in ens.items()}
    pf = res["phagefree"] = {"n_genomes": len(ens)}
    for mode in MODES:
        base = ek.shipped_calls(ens, seeds, mode)
        regs = {"shipped": base, "flat5kb": edge_snap(base, lengths, 5000), "frac": shipped_snap(base, lengths, wl)}
        pf[mode] = {v: {"calls_per_genome": len(r) / len(ens), "bp_called": int((r.end - r.start).sum())}
                    for v, r in regs.items()}
    Path(a.out).write_text(json.dumps(res, indent=1))
    for sname in [s.split("=", 1)[0] for s in sets]:
        print(f"== {sname}")
        for mode in MODES:
            for v in variants:
                r = res[sname][mode][v]
                print(f"  {mode:12s} {v:8s} complete {r['complete']:.3f} MAG {r['mag']:.3f} draft {r['draft']:.3f} "
                      f"calls/g {r['mag_calls_per_genome']:5.1f} false/Mb {r['mag_false_per_mb']:.2f} "
                      f"(complete {r['complete_false_per_mb']:.2f}) MAG bp {r['mag_bp_called']/1e6:6.2f} Mb "
                      f"off {r['mag_bp_off_element']/1e6:5.2f} IoU MAG {r['mag_median_iou']:.3f} "
                      f"complete {r['complete_median_iou']:.3f} complete bp {r['complete_bp_called']/1e6:.2f}")
    print("== phagefree")
    for mode in MODES:
        print("  " + mode + " " + "  ".join(f"{v} {x['calls_per_genome']:.2f}/g {x['bp_called']/1e6:.3f} Mb"
                                         for v, x in pf[mode].items()))
    print(f"[done] wrote {a.out}")


if __name__ == "__main__":
    main()
