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
"""predict = offline: the shipped CLI reproduces the offline rebuild of every boundary stage
(docs/att-boundaries.md; docs/project-principles.md method-change procedure, step 4).

  reference  calls rebuilt by evalkit.shipped_calls (refinement from stored fine tracks,
             att, composition edges, ranges) from Exploratory refinement-on runs with
             --write-seed-scores -> one TSV
  compare    each <outdir>/<mode>/<sample>/<genome>/regions.detailed.tsv written by
             `prophound predict` against the reference: coordinates, boundary and both
             edge sources must be identical

Run from prophound/models/prophage. The 2026-10-21 check (`arc_data/att_verify2/run.sh`):
  python scripts/boundary_verify.py reference --out analysis/att/verify2_reference.tsv
  python scripts/boundary_verify.py compare --reference analysis/att/verify2_reference.tsv \\
      --predict_dir ../../../arc_data/att_verify2
The 2026-10-26 Exploratory check (`arc_data/att_verify3/run.sh`):
  python scripts/boundary_verify.py reference --plan exploratory --out analysis/att/verify3_reference.tsv
  python scripts/boundary_verify.py compare --support --reference analysis/att/verify3_reference.tsv \\
      --predict_dir ../../../arc_data/att_verify3
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import evalkit as ek  # noqa: E402

RUNS = Path("../../../arc_data/att_runs")
GENOMES = {"tune42": "../../../arc_data/raw/genomes", "confirm18": "../../../arc_data/ncbi_confirm18",
           "lit7": "../../../arc_data/raw/genomes"}
PLANS = {  # sample -> modes
    "verify2": (("confirm18", ("default", "conservative")), ("lit7", ("default", "conservative")),
                ("tune42", ("default",))),
    # Exploratory composition edges and split-piece support (2026-10-26)
    "exploratory": (("confirm18", ("exploratory",)), ("lit7", ("exploratory",))),
}
KEY = ["genome_id", "contig_id", "start", "end", "boundary", "start_source", "end_source"]


def cmd_reference(a) -> None:
    out = []
    for sample, modes in PLANS[a.plan]:
        r = Path(a.runs) / sample
        ens = ek.load_tracks(str(r))
        seeds = [ek.load_tracks(str(r), f"prob_ckpt{k}") for k in (1, 2, 3)]
        fine = {d.name: pd.read_csv(d / "fine_windows.tsv", sep="\t", float_precision="round_trip")
                for d in r.iterdir() if d.is_dir()}
        for mode in modes:
            out.append(ek.shipped_calls(ens, seeds, mode, GENOMES[sample], fine_tracks=fine)
                       .assign(sample=sample, mode=mode))
    ref = pd.concat(out, ignore_index=True)
    ref.to_csv(a.out, sep="\t", index=False)
    print(f"[reference] {len(ref)} calls -> {a.out}")


def cmd_compare(a) -> None:
    ref = pd.read_csv(a.reference, sep="\t")
    total, bad = 0, 0
    for (sample, mode), r in ref.groupby(["sample", "mode"]):
        d = Path(a.predict_dir) / mode / sample
        if not (d / "batch_summary.tsv").exists():
            print(f"{mode:12} {sample:9} not finished")
            bad += 1
            continue
        parts = [pd.read_csv(g / "regions.detailed.tsv", sep="\t") for g in sorted(d.iterdir())
                 if g.is_dir() and (g / "regions.detailed.tsv").stat().st_size > 1]
        key = KEY + (["seed_support"] if a.support else [])
        cli = pd.concat(parts, ignore_index=True)[key].sort_values(key).reset_index(drop=True)
        off = r[key].sort_values(key).reset_index(drop=True)
        same = cli.equals(off)
        total += len(cli)
        bad += not same
        print(f"{mode:12} {sample:9} predict {len(cli):5} calls, offline {len(off):5}: identical = {same}")
        if not same:
            diff = cli.merge(off, how="outer", indicator=True)
            print(diff[diff._merge != "both"].head(20).to_string(index=False))
    verdict = "all identical" if not bad else f"{bad} sample(s) differ or unfinished"
    print(f"[compare] {total} calls compared; {verdict}")
    raise SystemExit(1 if bad else 0)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    r = sp.add_parser("reference")
    r.add_argument("--runs", default=str(RUNS))
    r.add_argument("--out", required=True)
    r.add_argument("--plan", choices=sorted(PLANS), default="verify2")
    r.set_defaults(fn=cmd_reference)
    c = sp.add_parser("compare")
    c.add_argument("--reference", required=True)
    c.add_argument("--predict_dir", required=True)
    c.add_argument("--support", action="store_true", help="also compare seed_support")
    c.set_defaults(fn=cmd_compare)
    a = ap.parse_args()
    a.fn(a)


if __name__ == "__main__":
    main()
