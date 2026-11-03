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
"""Operating modes under the shipped pipeline, rebuilt offline from score tracks.

For each run directory (`<genome>/windows.tsv` + `run_metadata.json`) and each
mode, calls are rebuilt the way predict_genome builds them (with runs labelled
r1, r2, r3 and ens, an extra `ens_shipped` row applies each mode's seed vote,
MODE_MIN_SEED_SUPPORT; the other rows are unvoted): full-length windows
with the mode's seed/extend, coverage builder and fragment merging; short-contig
windows (shorter than the run's window_len) per scale at
SHORT_CONTIG_THRESHOLDS[mode], no merging. Refinement is not applied (off by
default).

With --elements: covered and covered50 recall (overall and per stratum), calls
per genome, the on-target fraction of calls (a cost proxy, not precision), and
false calls per Mb of element-free sequence. Without it (a geNomad-negative
sample): calls per genome and per Mb, the fraction of bp called, and the calls'
max_prob -- which cannot distinguish false positives from elements geNomad
missed.

Example:
  python scripts/mode_recheck.py --elements scripts/eval_confirm18_elements.csv \\
      --genomes_dir ../raw/genomes --run ens=analysis/confirm18/ens --out analysis/confirm18/modes.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
from evalkit import (
    MODE_THRESHOLDS,
    calls,
    covered50,
    element_overlap_bp,
    iter_fasta_lengths,
    load_tracks,
    score_variant,
    shipped_calls,
)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--elements", default=None, help="omit for a geNomad-negative sample")
    ap.add_argument("--genomes_dir", default=None, help="<genome>.fna, for sequence lengths")
    ap.add_argument("--lengths_csv", default=None, help="alternative: accession,total_bp (e.g. summary.csv)")
    ap.add_argument("--run", action="append", required=True, metavar="LABEL=DIR")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    if not (args.genomes_dir or args.lengths_csv):
        raise SystemExit("need --genomes_dir or --lengths_csv")
    lengths = pd.read_csv(args.lengths_csv).set_index("accession").total_bp if args.lengths_csv else None
    ev = pd.read_csv(args.elements) if args.elements else None
    res, flat = {}, []
    specs = [s.partition("=")[::2] for s in args.run]
    loaded = {label: load_tracks(d) for label, d in specs}
    seeds = [loaded[s] for s in ("r1", "r2", "r3") if s in loaded]
    jobs = [(label, loaded[label], None) for label, _ in specs]
    if "ens" in loaded and len(seeds) == 3:   # the shipped pipeline: ensemble + the mode's vote
        jobs.append(("ens_shipped", loaded["ens"], seeds))
    for label, tracks, vote_seeds in jobs:
        if lengths is not None:
            total_bp = float(lengths.loc[list(tracks)].sum())
        else:
            total_bp = float(sum(n for g in tracks for _, n in iter_fasta_lengths(str(Path(args.genomes_dir) / f"{g}.fna"))))
        el_bp = float((ev.end - ev.start).sum()) if ev is not None else 0.0
        free_bp = total_bp - el_bp
        res[label] = {"n_genomes": len(tracks), "total_bp": total_bp, "modes": {}}
        for mode in MODE_THRESHOLDS:
            reg = (shipped_calls(tracks, vote_seeds, mode, args.genomes_dir) if vote_seeds
                   else calls(tracks, mode))
            n = len(reg)
            m = {"calls_per_genome": n / len(tracks), "short_contig_calls": int(reg["short"].sum()) if n else 0,
                 "calls_per_mb": n / (total_bp / 1e6)}
            if ev is not None:
                ov = [element_overlap_bp(ev, r) for r in reg.itertuples()]
                n_false = sum(o == 0 for o in ov)
                m.update({"covered_recall": score_variant(ev, reg)["covered_recall"] if n else 0.0,
                          "covered50_recall": covered50(ev, reg),
                          "on_target": 1 - n_false / max(1, n),
                          "false_calls_per_mb": n_false / (free_bp / 1e6),
                          "by_stratum": {str(s): {"n": int(len(sub)), "covered50_recall": covered50(sub, reg)}
                                         for s, sub in ev.groupby("stratum")}})
            else:
                bp = float((reg.end - reg.start).sum()) if n else 0.0
                m.update({"frac_bp_called": bp / total_bp,
                          "genomes_with_calls": int(reg.genome_id.nunique()) if n else 0,
                          "max_prob_median": float(reg.max_prob.median()) if n else None,
                          "calls_by_genome": reg.genome_id.value_counts().to_dict() if n else {}})
            res[label]["modes"][mode] = m
            flat.append({"run": label, "mode": mode, **{k: v for k, v in m.items() if not isinstance(v, dict)}})
            print(f"[{label} {mode}] " + " ".join(f"{k} {v:.3f}" if isinstance(v, float) else f"{k} {v}"
                                                  for k, v in m.items() if not isinstance(v, dict)), flush=True)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(res, indent=2) + "\n")
    pd.DataFrame(flat).to_csv(Path(args.out).with_suffix(".csv"), index=False)
    print(f"[done] wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
