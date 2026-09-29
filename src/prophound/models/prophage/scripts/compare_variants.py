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


"""Compare region-level metrics across inference variants (docs/scientific-validity.md).

Consumes directories of `<genome>/regions.tsv` produced by `prophound predict`
under different flags and reports, per variant, covered recall against geNomad
truth plus calls/genome, median IoU and fragmentation. It is deliberately
variant-agnostic: the ambiguity-policy job and the grid-phase job both use it.

`--group_by_phase` additionally reports *stability*: for each true element, in
how many variants it is covered, and for each call, in how many variants an
overlapping call appears. That is the number the grid-phase experiment needs --
a marginal call present in 1 of 5 phases is the failure mode, and a mean
recall would hide it.

Example:
  python scripts/compare_variants.py \
      --elements scripts/eval_sizestrat_elements.csv \
      --variant nchannel_r1=analysis/ambiguity_impact/nchannel_r1 \
      --variant uniform_r1=analysis/ambiguity_impact/uniform_r1 \
      --out analysis/ambiguity_impact/results.json
"""

from __future__ import annotations

import argparse
import glob
import json
import os
from collections import defaultdict

import numpy as np
import pandas as pd


def load_regions(d: str) -> pd.DataFrame:
    dfs = []
    for f in glob.glob(os.path.join(d, "*", "regions.tsv")):
        try:
            dfs.append(pd.read_csv(f, sep="\t"))
        except pd.errors.EmptyDataError:  # header-only: a genome with no calls
            pass
    if not dfs:
        return pd.DataFrame(columns=["genome_id", "contig_id", "start", "end", "max_prob"])
    return pd.concat(dfs, ignore_index=True)


def score_variant(ev: pd.DataFrame, reg: pd.DataFrame) -> dict:
    covered, ious, n_overlaps = [], [], []
    per_call: dict = {}
    for r in ev.itertuples():
        g = reg[(reg.genome_id == r.accession) & (reg.contig_id == r.contig)]
        o = g[(g.end > r.start) & (g.start < r.end)]
        covered.append(len(o) > 0)
        n_overlaps.append(len(o))
        for c in o.itertuples():
            per_call.setdefault((c.genome_id, c.contig_id, c.start, c.end), []).append(float(r.kb) if hasattr(r, "kb") else 0.0)
        if len(o):
            inter = np.minimum(o.end, r.end) - np.maximum(o.start, r.start)
            union = np.maximum(o.end, r.end) - np.minimum(o.start, r.start)
            ious.append(float((inter / union).max()))
        else:
            ious.append(0.0)
    n_genomes = ev.accession.nunique()
    # Conflation: one call claiming several distinct elements. Covered recall is
    # blind to this by construction (3 elements under 1 call score 3/3), so it is
    # reported explicitly -- see docs/downstream-postprocessing.md "Conflation".
    conflated = {k: v for k, v in per_call.items() if len(v) >= 2}
    absorbed = sum(len(v) for v in conflated.values())
    n_covered = int(np.sum(covered))
    return {
        "covered_recall": float(np.mean(covered)),
        "n_covered": n_covered,
        # matched recall: covered, minus elements that merely share a call with
        # another. The 1-1 achievable ceiling; the metric that *does* see conflation.
        "matched_recall": float((n_covered - (absorbed - len(conflated))) / max(1, len(ev))),
        "n_elements": int(len(ev)),
        "n_calls": int(len(reg)),
        "calls_per_genome": float(len(reg) / max(1, n_genomes)),
        "median_iou": float(np.median(ious)),
        "iou_ge_50": float(np.mean(np.asarray(ious) >= 0.5)),
        "fragmented_elements": int(np.sum(np.asarray(n_overlaps) > 1)),
        "conflated_calls": len(conflated),
        "absorbed_elements": absorbed,
        "_covered_vec": [bool(c) for c in covered],
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--elements", required=True, help="element manifest: accession,contig,start,end[,stratum]")
    ap.add_argument("--variant", action="append", required=True, metavar="LABEL=DIR",
                    help="repeatable; DIR holds <genome>/regions.tsv")
    ap.add_argument("--out", required=True)
    ap.add_argument("--group_by_phase", action="store_true",
                    help="also report per-element coverage stability across variants")
    ap.add_argument("--strata", action="store_true", help="break covered recall down by the stratum column")
    args = ap.parse_args()

    ev = pd.read_csv(args.elements)
    res: dict = {"n_elements": int(len(ev)), "variants": {}}
    cov_by_variant = {}
    for spec in args.variant:
        label, _, d = spec.partition("=")
        if not d:
            raise SystemExit(f"--variant needs LABEL=DIR, got {spec!r}")
        reg = load_regions(d)
        m = score_variant(ev, reg)
        cov_by_variant[label] = m.pop("_covered_vec")
        if args.strata and "stratum" in ev:
            byst = {}
            for st, sub in ev.groupby("stratum"):
                sm = score_variant(sub, reg)
                sm.pop("_covered_vec")
                byst[str(st)] = {k: sm[k] for k in ("covered_recall", "n_elements", "median_iou")}
            m["by_stratum"] = byst
        res["variants"][label] = m
        print(f"[{label}] covered {m['covered_recall']:.3f} / matched {m['matched_recall']:.3f} "
              f"({m['n_covered']}/{m['n_elements']}), {m['calls_per_genome']:.1f} calls/genome, "
              f"median IoU {m['median_iou']:.3f}, conflated {m['conflated_calls']} calls "
              f"absorbing {m['absorbed_elements']}", flush=True)

    if args.group_by_phase and len(cov_by_variant) > 1:
        labels = list(cov_by_variant)
        mat = np.array([cov_by_variant[l] for l in labels])  # variants x elements
        per_element = mat.sum(axis=0)
        hist = {int(k): int(v) for k, v in zip(*np.unique(per_element, return_counts=True))}
        res["stability"] = {
            "n_variants": len(labels),
            "elements_covered_in_all": int((per_element == len(labels)).sum()),
            "elements_covered_in_none": int((per_element == 0).sum()),
            "elements_unstable": int(((per_element > 0) & (per_element < len(labels))).sum()),
            "coverage_count_histogram": hist,
            "note": ("An element covered in some variants but not others is phase- or "
                     "policy-dependent. For the grid-phase experiment, elements_unstable "
                     "is the headline: it is invisible in any single-variant recall figure."),
        }
        print(f"\n[stability] covered in all {len(labels)}: {res['stability']['elements_covered_in_all']}; "
              f"unstable: {res['stability']['elements_unstable']}; "
              f"never covered: {res['stability']['elements_covered_in_none']}", flush=True)

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(res, f, indent=2)
    print(f"\n[done] wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
