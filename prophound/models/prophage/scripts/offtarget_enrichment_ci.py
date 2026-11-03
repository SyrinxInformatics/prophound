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
"""Interval and test for the off-target hallmark enrichment (analysis/offtarget/,
docs/research-integrity.md "Off-target enrichment: controls", docs/release-text-0.1.0.md section 8). Per call group
(Default, tuning sample), the hallmark rate against its length-matched controls
drawn clear of geNomad elements: ratio with a 95% genome bootstrap (10,000
resamples; calls and controls resampled together by genome), plus a one-sided
Fisher exact test. Fisher treats controls as independent, so the bootstrap is the
interval to quote.

Run from prophound/models/prophage:
  python scripts/offtarget_enrichment_ci.py \\
      --annotation analysis/offtarget/annotation_ctrl_elements_only.csv \\
      --out analysis/offtarget/enrichment_ci.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

GROUPS = [("off_3seed", "control_off_3seed"), ("off_1_2seed", "control_off_1_2seed"),
          ("on_target", "control_on_target")]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--annotation", required=True)
    ap.add_argument("--mode", default="default")
    ap.add_argument("--seed", type=int, default=20260929)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    d = pd.read_csv(a.annotation)
    d = d[d["mode"] == a.mode].assign(hm=lambda x: x.category != "no_mge_hallmark")
    rng = np.random.default_rng(a.seed)
    res = {}
    for g, c in GROUPS:
        x, y = d[d.group == g], d[d.group == c]
        kx, nx, ky, ny = int(x.hm.sum()), len(x), int(y.hm.sum()), len(y)
        p = float(fisher_exact([[kx, nx - kx], [ky, ny - ky]], alternative="greater")[1])
        gen = sorted(set(x.genome_id) | set(y.genome_id))
        X = x.groupby("genome_id").hm.agg(["sum", "count"]).reindex(gen, fill_value=0).to_numpy()
        Y = y.groupby("genome_id").hm.agg(["sum", "count"]).reindex(gen, fill_value=0).to_numpy()
        idx = rng.integers(0, len(gen), size=(10000, len(gen)))
        rx = X[idx, 0].sum(1) / np.maximum(1, X[idx, 1].sum(1))
        ry = Y[idx, 0].sum(1) / np.maximum(1, Y[idx, 1].sum(1))
        ok = ry > 0
        lo, hi = np.percentile(rx[ok] / ry[ok], [2.5, 97.5])
        res[g] = {"calls": [kx, nx], "controls": [ky, ny], "rate": kx / nx, "control_rate": ky / ny,
                  "ratio": (kx / nx) / (ky / ny), "ratio_genome_boot95": [float(lo), float(hi)],
                  "fisher_one_sided_p": p, "n_genomes": int(x.genome_id.nunique())}
        print(f"{g}: {kx}/{nx} = {kx / nx:.3f} vs {ky}/{ny} = {ky / ny:.3f}; ratio {res[g]['ratio']:.2f} "
              f"[{lo:.2f}, {hi:.2f}]; Fisher p = {p:.1e}; {res[g]['n_genomes']} genomes")
    Path(a.out).write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
