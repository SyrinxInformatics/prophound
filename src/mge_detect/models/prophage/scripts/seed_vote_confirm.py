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
"""Call-level seed vote, re-measured offline (docs/error-analysis.md).

For each seed threshold t, calls are rebuilt from the ensemble track and from
each single-seed track with the shipped pipeline (mode_recheck.calls logic at
seed t, extend min(t, 0.99)). An ensemble call is kept at vote k if >= k seeds'
own calls at the same t overlap it (predict_genome.seed_support). k=1 is `max`.

Short-contig windows use the thresholds of the mode whose seed is the largest
not exceeding t (exploratory below 0.995).

Example:
  python scripts/seed_vote_confirm.py --elements scripts/eval_confirm18_elements.csv \\
      --root analysis/confirm18 --out analysis/confirm18/seed_vote.csv
  python scripts/seed_vote_confirm.py --root analysis/phagefree --out analysis/phagefree/seed_vote.csv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mode_recheck as mr  # noqa: E402
from refine_sweep import covered50  # noqa: E402

THRESHOLDS = [0.9, 0.95, 0.98, 0.99, 0.995, 0.999]


def calls_at(tracks: dict, t: float) -> pd.DataFrame:
    mode = max((m for m, s in mr.MODE_THRESHOLDS.items() if s <= t), key=mr.MODE_THRESHOLDS.get,
               default="exploratory")
    saved = mr.MODE_THRESHOLDS[mode]
    mr.MODE_THRESHOLDS[mode] = t   # calls() reads the seed from here; restore after
    try:
        return mr.calls(tracks, mode)
    finally:
        mr.MODE_THRESHOLDS[mode] = saved


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", required=True, help="dir holding r1/ r2/ r3/ ens/ run directories")
    ap.add_argument("--elements", default=None, help="omit for a geNomad-negative sample")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    ev = pd.read_csv(args.elements) if args.elements else None
    tr = {r: mr.load_tracks(str(Path(args.root) / r)) for r in ("r1", "r2", "r3", "ens")}
    n_gen = len(tr["ens"])
    rows = []
    for t in THRESHOLDS:
        ens = calls_at(tr["ens"], t)
        seeds = [calls_at(tr[r], t) for r in ("r1", "r2", "r3")]
        sup = []
        for g, sub in ens.groupby("genome_id", sort=False):
            per = [s[s.genome_id == g] for s in seeds]
            sup += list(zip(sub.index, mr_support(sub, per)))
        ens = ens.assign(support=pd.Series(dict(sup)))
        for k in (1, 2, 3):
            reg = ens[ens.support >= k]
            row = {"threshold": t, "k": k, "calls_per_genome": len(reg) / n_gen}
            if ev is not None:
                off = sum(mr.element_overlap_bp(ev, r) == 0 for r in reg.itertuples())
                row.update({"covered50": covered50(ev, reg), "off_target_per_genome": off / n_gen})
            rows.append(row)
            print(" ".join(f"{a} {b:.3f}" if isinstance(b, float) else f"{a} {b}" for a, b in row.items()),
                  flush=True)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.out, index=False)
    print(f"[done] wrote {args.out}")


def mr_support(sub: pd.DataFrame, per: list) -> list:
    try:
        from mge_detect.models.prophage.src.predict_genome import seed_support
    except ImportError:
        from src.predict_genome import seed_support
    return seed_support(sub, per)


if __name__ == "__main__":
    main()
