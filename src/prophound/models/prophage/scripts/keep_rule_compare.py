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
"""Descriptive comparison of call-keeping rules on fixed candidate pools (user, 2026-10-20).

docs/error-analysis.md "Combined keep score vs the vote". Run from the repo root:
  pixi run python src/prophound/models/prophage/scripts/keep_rule_compare.py \
      src/prophound/models/prophage/analysis/keep_rules/results.csv
Tracks: the 42 (analysis/multiscale/*_3000), the 18 (analysis/confirm18) and the
phage-free 20 (analysis/phagefree). Hallmark counts use analysis/offtarget/annotation.csv,
which covers exactly the 0.95 pool, so they are reported for Default only.

Pool = ensemble calls at a mode's seed threshold (k=1), rebuilt with mode_recheck.calls.
Rules: the shipped seed vote (k), mean_prob cut, a combined logistic score
(support, mean_prob, log length, max_prob) fit on the 42 and applied unchanged to
the 18 and phage-free, and a rescue rule (3 of 3, or mean_prob >= c).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler

M = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(M / "scripts"))
sys.path.insert(0, str(M))
import mode_recheck as mr  # noqa: E402
from refine_sweep import covered50  # noqa: E402
from seed_vote_confirm import calls_at, mr_support  # noqa: E402

OUT = Path(sys.argv[1])
SAMPLES = {
    "tune42": ({r: M / f"analysis/multiscale/{r}_3000" for r in ("r1", "r2", "r3", "ens")},
               M / "scripts/eval_sizestrat_elements.csv"),
    "check18": ({r: M / f"analysis/confirm18/{r}" for r in ("r1", "r2", "r3", "ens")},
                M / "scripts/eval_confirm18_elements.csv"),
    "phagefree": ({r: M / f"analysis/phagefree/{r}" for r in ("r1", "r2", "r3", "ens")}, None),
}
POOLS = {"exploratory": 0.90, "default": 0.95, "conservative": 0.995}
FEATS = ["support", "mean_prob", "max_prob", "log_len"]
COMBO_CUTS = sorted(set(np.round(np.r_[np.arange(0.005, 0.1, 0.005), np.arange(0.1, 0.96, 0.05)], 3)))


def pool(tr: dict, t: float, ev: pd.DataFrame | None) -> pd.DataFrame:
    ens = calls_at(tr["ens"], t)
    seeds = [calls_at(tr[r], t) for r in ("r1", "r2", "r3")]
    sup = {}
    for g, sub in ens.groupby("genome_id", sort=False):
        sup.update(zip(sub.index, mr_support(sub, [s[s.genome_id == g] for s in seeds])))
    ens = ens.assign(support=pd.Series(sup), log_len=np.log10((ens.end - ens.start).clip(lower=1)))
    if ev is not None:
        ov = np.array([mr.element_overlap_bp(ev, r) for r in ens.itertuples()])
        ens = ens.assign(on_target=ov > 0, off_bp=(ens.end - ens.start) - ov)
    else:
        ens = ens.assign(on_target=False, off_bp=ens.end - ens.start)
    return ens.reset_index(drop=True)


def metrics(kept: pd.DataFrame, ev, n_gen: int, hallmark: pd.DataFrame | None) -> dict:
    m = {"calls": len(kept) / n_gen, "bp_kb": (kept.end - kept.start).sum() / n_gen / 1e3,
         "off_bp_kb": kept.off_bp.sum() / n_gen / 1e3}
    if ev is not None:
        m["covered50"] = covered50(ev, kept)
        m["off_calls"] = (~kept.on_target).sum() / n_gen
    if hallmark is not None:
        k = kept.merge(hallmark, on=["genome_id", "contig_id", "start", "end"], how="inner")
        m["off_hallmark_kept"] = int(((~k.on_target) & k.hallmark).sum())
    return m


def main() -> None:
    data = {}
    for name, (dirs, elf) in SAMPLES.items():
        tr = {r: mr.load_tracks(str(d)) for r, d in dirs.items()}
        ev = pd.read_csv(elf) if elf else None
        data[name] = (tr, ev, len(tr["ens"]))
    ann = pd.read_csv(M / "analysis/offtarget/annotation.csv")
    ann = ann[(ann["mode"] == "exploratory") & ~ann.group.str.startswith("control")]
    hall = ann.assign(hallmark=ann.category != "no_mge_hallmark")[["genome_id", "contig_id", "start", "end",
                                                                   "hallmark"]]
    rows = []
    for mode, t in POOLS.items():
        pools = {n: pool(tr, t, ev) for n, (tr, ev, _) in data.items()}
        tune = pools["tune42"]
        # combined score: fit on the 42; out-of-fold scores for the 42 itself (grouped by genome)
        X, y, grp = tune[FEATS].to_numpy(float), tune.on_target.to_numpy(), tune.genome_id.to_numpy()
        oof = np.zeros(len(tune))
        for tr_i, te_i in GroupKFold(5).split(X, y, grp):
            sc = StandardScaler().fit(X[tr_i])
            fit = LogisticRegression().fit(sc.transform(X[tr_i]), y[tr_i])
            oof[te_i] = fit.predict_proba(sc.transform(X[te_i]))[:, 1]
        sc = StandardScaler().fit(X)
        lr = LogisticRegression().fit(sc.transform(X), y)
        print(f"[{mode}] coef {dict(zip(FEATS, np.round(lr.coef_[0], 2)))}", flush=True)
        for n, p in pools.items():
            p["combo"] = oof if n == "tune42" else lr.predict_proba(sc.transform(p[FEATS].to_numpy(float)))[:, 1]
        for n, p in pools.items():
            _, ev, ng = data[n]
            rules = [("vote", k, p.support >= k) for k in (1, 2, 3)]
            for c in np.round(np.arange(0.30, 0.86, 0.05), 2):
                rules += [("mean_prob", c, p.mean_prob >= c),
                          ("rescue3_or_mean", c, (p.support >= 3) | (p.mean_prob >= c))]
            rules += [("combo", c, p.combo >= c) for c in COMBO_CUTS]
            for rule, param, keep in rules:
                rows.append({"mode": mode, "t": t, "sample": n, "rule": rule, "param": param,
                             **metrics(p[keep], ev, ng, hall if (t == 0.95 and n == "tune42") else None)})
        print(f"[{mode}] done", flush=True)
    pd.DataFrame(rows).to_csv(OUT, index=False)
    print(f"[done] {OUT}")


if __name__ == "__main__":
    main()
