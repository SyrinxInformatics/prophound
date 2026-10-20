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
"""Chance-adjusted covered50 on LAMBDA (docs/claims-and-evidence.md "LAMBDA result",
docs/release-text-v1.md section 2). Per tool, every call on a chromosome is rotated
by one uniform random offset (circularly, wrapping calls split at the origin). That
keeps the tool's call lengths, count and spacing, and moves them off the elements.
Chance covered50 is the mean over draws, and recall above chance is observed minus
chance. prophound Default and Conservative are compared with each annotation tool on
the genomes both ran: the difference in recall above chance, with a paired 95% genome
bootstrap (per-genome observed hits minus mean chance hits, 10,000 resamples). Also reports the fraction of each chromosome a tool calls. Lengths are
LAMBDA's `genome_length` (the true length rounded up to 1 kb, plus 1 kb), which
shifts chance by well under 1%.

Run from src/prophound/models/prophage:
  python scripts/lambda_chance.py --data analysis/external_benchmarks/lambda \\
      --prophound_run analysis/lambda/ens --out analysis/external_benchmarks/lambda/chance.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from lambda_score import prophound_calls, union_cov

TOOLS = ["prophound exploratory", "prophound default", "prophound conservative",
         "geNomad", "PHASTER", "PIDE", "VIBRANT"]


def union_len(c: pd.DataFrame) -> int:
    tot, cs, ce = 0, None, None
    for s, e in c.sort_values("start")[["start", "end"]].to_numpy():
        if ce is None or s > ce:
            tot += 0 if ce is None else ce - cs
            cs, ce = s, e
        else:
            ce = max(ce, e)
    return int(tot + (0 if ce is None else ce - cs))


def rotate(c: pd.DataFrame, L: int, off: int) -> pd.DataFrame:
    s, e = (c.start.to_numpy() + off) % L, (c.end.to_numpy() - c.start.to_numpy())
    rows = []
    for a, n in zip(s, e):
        if a + n <= L:
            rows.append((a, a + n))
        else:
            rows += [(a, L), (0, a + n - L)]
    return pd.DataFrame(rows, columns=["start", "end"])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", required=True)
    ap.add_argument("--prophound_run", required=True)
    ap.add_argument("--draws", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=20260929)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    D = Path(a.data)
    G = pd.read_csv(D / "lambda_genomes.tsv", sep="\t")
    T = pd.read_csv(D / "lambda_truth.tsv", sep="\t")
    C = pd.read_csv(D / "lambda_calls.tsv", sep="\t")
    RAN = pd.read_csv(D / "lambda_tools_run.tsv", sep="\t")
    tools = {t: C[C.tool == t] for t in TOOLS if t in set(C.tool)}
    ran = {t: set(RAN[RAN.tool == t].assembly) for t in tools}
    for k, v in prophound_calls(a.prophound_run, G).items():
        tools[k], ran[k] = v, set(G.assembly)
    rng = np.random.default_rng(a.seed)
    res, excess = {}, {}
    for tool in [t for t in TOOLS if t in tools]:
        gs = G[G.assembly.isin(ran[tool])]
        t = T[T.assembly.isin(gs.assembly)].reset_index(drop=True)
        c = tools[tool][tools[tool].assembly.isin(gs.assembly)]
        obs = float((union_cov(t, c) >= 0.5).mean())
        frac = float(sum(union_len(c[c.sequence == s]) for s in gs.sequence) / gs.genome_length.sum())
        chance, per_hits = [], np.zeros(len(t))
        for _ in range(a.draws):
            parts = []
            for r in gs.itertuples():
                cc = c[c.sequence == r.sequence]
                if len(cc):
                    parts.append(rotate(cc, int(r.genome_length), int(rng.integers(0, r.genome_length))).assign(sequence=r.sequence))
            rc = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=["start", "end", "sequence"])
            h = union_cov(t, rc) >= 0.5
            per_hits += h
            chance.append(float(h.mean()))
        ch = float(np.mean(chance))
        res[tool] = {"n_genomes": len(gs), "n_elements": len(t), "covered50": obs, "fraction_called": frac,
                     "chance_covered50": ch, "chance_95": [float(np.percentile(chance, 2.5)), float(np.percentile(chance, 97.5))],
                     "above_chance": obs - ch}
        excess[tool] = pd.DataFrame({"assembly": t.assembly, "excess": (union_cov(t, c) >= 0.5) - per_hits / a.draws})
        print(f"{tool:24s} genomes {len(gs):2d} covered50 {obs:.3f} called {100 * frac:5.1f}% "
              f"chance {ch:.3f} [{res[tool]['chance_95'][0]:.3f}, {res[tool]['chance_95'][1]:.3f}] above chance {obs - ch:.3f}")
    res["paired_above_chance"] = {}
    for ours in ("prophound default", "prophound conservative"):
        for other in [t for t in TOOLS if not t.startswith("prophound") and t in excess]:
            gs = sorted(ran[ours] & ran[other])
            A = excess[ours].groupby("assembly").excess.sum().reindex(gs, fill_value=0).to_numpy()
            B = excess[other].groupby("assembly").excess.sum().reindex(gs, fill_value=0).to_numpy()
            n = T[T.assembly.isin(gs)].groupby("assembly").size().reindex(gs, fill_value=0).to_numpy()
            idx = rng.integers(0, len(gs), size=(10000, len(gs)))
            d = (A[idx].sum(1) - B[idx].sum(1)) / np.maximum(1, n[idx].sum(1))
            est = float((A.sum() - B.sum()) / n.sum())
            lo, hi = np.percentile(d, [2.5, 97.5])
            res["paired_above_chance"][f"{ours} - {other}"] = {"n_genomes": len(gs), "diff": est, "boot95": [float(lo), float(hi)]}
            print(f"{ours} - {other}: above-chance difference {est:+.3f} [{lo:+.3f}, {hi:+.3f}] on {len(gs)} genomes")
    Path(a.out).write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
