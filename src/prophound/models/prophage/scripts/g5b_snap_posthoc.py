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
"""Post-hoc, not evidence for G5b: the edge snap's effect on the G5b tracks.
Rebuilds prophound calls with the snap off and at the shipped 0.5 x window, then
reports per condition and mode the elements whose coverage changed, covered50
off -> on, and (Default) the sequence added, split on/off element. 2026-10-16.

Run from src/prophound/models/prophage:
  python scripts/g5b_snap_posthoc.py
"""

from __future__ import annotations

import g5_score as g
import numpy as np
import pandas as pd

CONDS = ["complete", "mag_s0", "mag_s1", "mag_s2", "draft_s0"]
MODES = ["exploratory", "default", "conservative"]


def union_bp(c: pd.DataFrame) -> int:
    tot = 0
    for _, x in c.groupby(["assembly", "sequence"]):
        cs = ce = None
        for s, e in x.sort_values("start")[["start", "end"]].to_numpy():
            if ce is None or s > ce:
                tot += 0 if ce is None else ce - cs
                cs, ce = s, e
            else:
                ce = max(ce, e)
        tot += 0 if ce is None else ce - cs
    return int(tot)


def main() -> None:
    T = pd.read_csv("analysis/g5b/g5b_truth.tsv", sep="\t")
    lengths = pd.read_csv("analysis/g5b/contig_lengths.tsv", sep="\t")
    L = (T.end - T.start).to_numpy()
    for cond in CONDS:
        run = f"analysis/g5b/prophound/{cond}"
        off, on = g.prophound_calls(cond, run), g.prophound_calls(cond, run, lengths, 0.5)
        for mode in MODES:
            a, b = off[off.tool == f"prophound {mode}"], on[on.tool == f"prophound {mode}"]
            ca, cb = g.cov_frac(T, a), g.cov_frac(T, b)
            d = cb - ca
            line = (f"{cond:9s} {mode:12s} calls {len(b):4d} | elements with coverage changed "
                    f"{int((np.abs(d) > 1e-9).sum()):2d}, max gain {d.max():.3f}, "
                    f"covered50 {int((ca >= .5).sum())} -> {int((cb >= .5).sum())}")
            if mode == "default":
                ua, ub = union_bp(a), union_bp(b)
                on_el = ((cb - ca) * L).sum()
                line += (f" | called {ua / 1e3:.0f} -> {ub / 1e3:.0f} kb (+{100 * (ub - ua) / ua:.1f}%), "
                         f"on-element {on_el / 1e3:.0f} of {(ub - ua) / 1e3:.0f} kb added")
            print(line)


if __name__ == "__main__":
    main()
