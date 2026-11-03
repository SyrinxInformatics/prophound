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
"""Where composition edges reach, and what extending them to Exploratory would cost.

Composition edges in every mode (analysis/att/). For each sample, rebuilds every mode's calls
with the shipped pipeline (evalkit.shipped_calls on a refinement-on run's tracks) twice:
composition as shipped (`auto`: Default and Conservative) and composition `on` (all modes).
Per mode and variant it reports the method-change costs (docs/project-principles.md): covered50, covered,
matched, median IoU, absorbed, Mb called and off-element, calls per genome, edge error,
edge-range coverage (calls on exactly one element, estimated edges not on a contig end),
and each edge's source. For the composition stage's input (the `off` calls) it also tags
why each estimated edge is or is not moved: contig end, call under MIN_CALL_BP, host flank
under MIN_HOST_BP, profiles closer than MIN_JSD, change-point on the edge, moved, or moved
to the inward ceiling.

Example (from prophound/models/prophage):
  python scripts/composition_coverage.py --sample tune42 ../../../arc_data/att_runs/tune42 \\
      ../../../arc_data/raw/genomes scripts/eval_sizestrat_elements.csv \\
      --out analysis/att/composition_coverage
"""

from __future__ import annotations

import argparse
from pathlib import Path

import att_boundaries as ab
import evalkit as ek
import numpy as np
import pandas as pd

try:
    from prophound.models.prophage.src import composition as comp
except ImportError:
    from src import composition as comp

MODES = ["exploratory", "default", "conservative"]


def range_coverage(calls: pd.DataFrame, truth: pd.DataFrame) -> tuple[int, int]:
    """(true edges inside their range, edges checked): estimated edges, not on a contig end,
    of calls overlapping exactly one element."""
    hit = n = 0
    for c in calls.itertuples():
        t = truth[(truth.accession == c.genome_id) & (truth.contig == c.contig_id)
                  & (truth.end > c.start) & (truth.start < c.end)]
        if len(t) != 1 or c.boundary == "att":
            continue
        for side, true in (("start", int(t.start.iloc[0])), ("end", int(t.end.iloc[0]))):
            if c.contig_edge in (side, "both"):
                continue
            lo, hi = int(getattr(c, f"{side}_lo")), int(getattr(c, f"{side}_hi"))
            n += 1
            hit += lo <= true <= hi
    return hit, n


def skip_reasons(calls: pd.DataFrame, seqs_of) -> pd.Series:
    """Why the composition stage moves, or leaves, each edge of its input calls."""
    out = []
    for r in calls.itertuples(index=False):
        s, e = int(r.start), int(r.end)
        L = e - s
        for side in ("start", "end"):
            if r.boundary == "att":
                out.append("att")
                continue
            if r.contig_edge in (side, "both"):
                out.append("contig_end")
                continue
            if L < comp.MIN_CALL_BP:
                out.append("short_call")
                continue
            contig = seqs_of(r.genome_id)[r.contig_id]
            lo, hi = int(getattr(r, f"{side}_lo")), int(getattr(r, f"{side}_hi"))
            host_seq = contig[max(0, lo - comp.HOST_FLANK):lo] if side == "start" else contig[hi:hi + comp.HOST_FLANK]
            if len(host_seq) < comp.MIN_HOST_BP:
                out.append("host_flank")
                continue
            elem = comp.profile(comp.kmer_codes(contig[s + int(0.2 * L):e - int(0.2 * L)]))
            host = comp.profile(comp.kmer_codes(host_seq))
            if comp.jsd(elem, host) < comp.MIN_JSD:
                out.append("low_jsd")
                continue
            old = s if side == "start" else e
            new = comp.change_point(contig, lo, hi, elem, host, side, old)
            inward = (new - old) if side == "start" else (old - new)
            if inward > min(comp.MAX_INWARD_FRAC * L, comp.MAX_INWARD_BP):
                out.append("moved_to_ceiling")
            else:
                out.append("moved" if int(np.round(new)) != old else "no_change_point")
    return pd.Series(out, dtype=object)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", nargs=4, action="append", required=True,
                    metavar=("NAME", "RUN_DIR", "GENOMES_DIR", "ELEMENTS"),
                    help="a refinement-on run, its <genome>.fna directory and its truth CSV; repeatable")
    ap.add_argument("--strip_contig_version", action="append", default=[],
                    help="sample whose truth contigs carry no '.N' version (LAMBDA)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rows, reasons = [], []
    for name, run, gdir, elements in a.sample:
        ens = ek.load_tracks(run)
        seeds = [ek.load_tracks(run, f"prob_ckpt{k}") for k in (1, 2, 3)]
        fine = {g: pd.read_csv(Path(run) / g / "fine_windows.tsv", sep="\t", float_precision="round_trip")
                for g in ens}
        ev = pd.read_csv(elements)
        ev = ev[ev.accession.isin(set(ens))]
        cache: dict = {}

        def seqs_of(g, gdir=gdir, cache=cache):
            if g not in cache:
                cache.clear()
                cache[g] = ab.load_genome(Path(gdir) / f"{g}.fna")[0]
            return cache[g]

        for mode in MODES:
            variants = ("off", "on") if mode == "exploratory" else ("auto",)  # on == auto in the others
            built = {v: ek.shipped_calls(ens, seeds, mode, gdir, fine_tracks=fine, composition=v)
                     for v in variants + (() if mode == "exploratory" else ("off",))}
            for v in variants:
                c = built[v]
                if name in a.strip_contig_version:
                    c = c.assign(contig_id=c.contig_id.astype(str).str.split(".").str[0])
                hit, n = range_coverage(c, ev)
                src = pd.concat([c.start_source, c.end_source]).value_counts()
                rows.append({"sample": name, "mode": mode, "composition": "on" if v != "off" else "off",
                             **ab.metrics(ev, c, ek), "calls_per_genome": len(c) / len(ens),
                             "range_cov": hit / max(1, n), "range_n": n,
                             **{f"src_{k}": int(src.get(k, 0)) for k in
                                ("att", "composition", "fine", "coarse", "contig_end")}})
                print(f"[{name} {mode} {v}] " + ", ".join(f"{k} {x:.3f}" if isinstance(x, float) else f"{k} {x}"
                                                          for k, x in rows[-1].items()), flush=True)
            r = skip_reasons(built["off"].sort_values("genome_id", kind="stable").reset_index(drop=True), seqs_of).value_counts()
            reasons.append({"sample": name, "mode": mode, "edges": int(r.sum()), **r.to_dict()})
            print(f"[{name} {mode}] reasons {reasons[-1]}", flush=True)
    pd.DataFrame(rows).to_csv(f"{a.out}.csv", index=False)
    pd.DataFrame(reasons).fillna(0).to_csv(f"{a.out}.reasons.csv", index=False)
    print(f"[done] wrote {a.out}.csv and {a.out}.reasons.csv")


if __name__ == "__main__":
    main()
