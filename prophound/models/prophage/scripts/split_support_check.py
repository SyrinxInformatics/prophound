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
"""seed_support before and after the split-piece fix (docs/lockbox2.md, amendment 2026-10-26).

Rebuilds every mode with the shipped pipeline on the four boundary samples, once with
calls.vote_extent replaced by the identity (the old count, on the piece itself) and once
as shipped; asserts the coordinates are identical and prints, per sample and mode, how
many calls change support and how many read below the mode's vote. Writes the changed
calls to the given TSV. Run from prophound/models/prophage:
  python scripts/split_support_check.py analysis/att/split_support_fix.changed.tsv
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, "scripts")
import evalkit as ek  # noqa: E402
import prophound.models.prophage.src.calls as calls_mod  # noqa: E402

RUNS = Path("../../../arc_data/att_runs")
SAMPLES = {"tune42": "../../../arc_data/raw/genomes", "confirm18": "../../../arc_data/ncbi_confirm18",
           "lit7": "../../../arc_data/raw/genomes", "lambda": "analysis/g5/fasta"}
VOTE = {"exploratory": 1, "default": 3, "conservative": 3}
fixed_fn = calls_mod.vote_extent
rows, detail = [], []
for s, gdir in SAMPLES.items():
    r = RUNS / s
    ens = ek.load_tracks(str(r))
    seeds = [ek.load_tracks(str(r), f"prob_ckpt{k}") for k in (1, 2, 3)]
    fine = {g: pd.read_csv(r / g / "fine_windows.tsv", sep="\t", float_precision="round_trip") for g in ens}
    for mode in VOTE:
        out = {}
        for label, fn in (("old", lambda d: d), ("new", fixed_fn)):
            calls_mod.vote_extent = fn
            out[label] = ek.shipped_calls(ens, seeds, mode, gdir, fine_tracks=fine).reset_index(drop=True)
        calls_mod.vote_extent = fixed_fn
        o, n = out["old"], out["new"]
        key = ["genome_id", "contig_id", "start", "end"]
        assert o[key].equals(n[key]), (s, mode)
        st = n.get("refine_status", pd.Series("?", index=n.index))
        changed = o.seed_support != n.seed_support
        rows.append({"sample": s, "mode": mode, "calls": len(n), "split": int((st == "split").sum()),
                     "changed": int(changed.sum()), "changed_split": int((changed & (st == "split")).sum()),
                     "changed_refined": int((changed & (st == "refined")).sum()),
                     "below_vote_old": int((o.seed_support < VOTE[mode]).sum()),
                     "below_vote_new": int((n.seed_support < VOTE[mode]).sum()),
                     "lowered": int((n.seed_support < o.seed_support).sum())})
        d = n[changed].assign(sample=s, mode=mode, old=o.seed_support[changed])
        detail.append(d[["sample", "mode", *key, "refine_status", "old", "seed_support"]])
print(pd.DataFrame(rows).to_string(index=False))
pd.concat(detail).to_csv(sys.argv[1], sep="\t", index=False)
