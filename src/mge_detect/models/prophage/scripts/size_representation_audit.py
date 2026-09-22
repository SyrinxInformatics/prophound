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


"""Are small elements represented as positives in training, and at what dilution?

`docs/v1-scope.md` gating audit A.1. This is the central *model-side* diversity
question and it gates the retraining data design, so it runs before any training
code is written.

The problem it settles: `window_len` is 10,000 bp, so no window can fit inside an
element smaller than that. A 5 kb element's best window is 50% host sequence, and
measured on real data the best achievable is 66.8% element content at 5-8 kb. What
is unknown -- because the BED->window labelling step is not in this repo -- is
whether such windows were labelled positive at all, and if so how diluted the
positive class is at small sizes.

Reports, per element-size stratum:
  * number of positive training windows whose overlapping element is that size
  * the distribution of *element content* in those windows (what fraction of the
    window is actually element)
  * positives per element, i.e. how many training examples each element generated

How to read it:
  * If small elements produce few or no positive windows, the model never learned
    them and no downstream adaptation can fully recover it -- the retraining needs
    variable-length windows so a 3 kb window can sit inside a 5 kb element.
  * If they produce plenty of positive windows but at low element content, the
    model learned them from mostly-host examples, which predicts exactly the weak,
    smeared signal observed (0.857 recall, IoU 0.334 at 5-8 kb) -- again arguing
    for variable-length training, but for a different reason.
  * If they are well represented at high content, the small-element gap is
    architectural rather than data-driven, and downstream multi-scale is the right
    fix instead.

Reads FASTA-free: needs only the manifest coordinates and the BED positives.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

STRATA = [(0, 5000, "<5kb"), (5000, 8000, "5-8kb"), (8000, 10000, "8-10kb"),
          (10000, 15000, "10-15kb"), (15000, 30000, "15-30kb"),
          (30000, 60000, "30-60kb"), (60000, 10**9, "60kb+")]


def stratum_of(n: int) -> str:
    for lo, hi, lab in STRATA:
        if lo <= n < hi:
            return lab
    return "?"


def read_bed_positives(bed_path: Path):
    out = []
    for line in open(bed_path):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 4 or parts[3] != "positive":
            continue
        out.append((parts[0], int(parts[1]), int(parts[2])))
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--beds_dir", required=True, help="beds_pos/")
    ap.add_argument("--splits_dir", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--split", default="train")
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    print("[data] reading manifest", flush=True)
    m = pd.read_csv(args.manifest, usecols=["genome_id", "contig_id", "start", "end", "label", "neg_type"])
    ids = {l.strip() for l in open(Path(args.splits_dir) / f"{args.split}.txt") if l.strip()}
    m = m[m.genome_id.isin(ids)]
    pos = m[m.label == 1]
    print(f"[data] {len(m)} {args.split} windows, {len(pos)} positive "
          f"({len(pos) / max(1, len(m)) * 100:.2f}%)", flush=True)

    # element intervals per genome, from the BEDs
    elements: dict = {}
    missing = 0
    for g in sorted(ids):
        bed = Path(args.beds_dir) / f"{g}.bed"
        if not bed.exists():
            missing += 1
            continue
        for contig, s, e in read_bed_positives(bed):
            elements.setdefault((g, contig), []).append((s, e))
    n_el = sum(len(v) for v in elements.values())
    print(f"[data] {n_el} elements from BEDs ({missing} genomes had no BED)", flush=True)

    per_stratum = defaultdict(lambda: {"n_windows": 0, "content": []})
    per_element = defaultdict(int)
    unmatched = 0
    for r in pos.itertuples():
        ivs = elements.get((r.genome_id, r.contig_id))
        if not ivs:
            unmatched += 1
            continue
        ws, we = int(r.start), int(r.end)
        wlen = max(1, we - ws)
        best = None
        for (s, e) in ivs:
            ov = min(we, e) - max(ws, s)
            if ov > 0 and (best is None or ov > best[0]):
                best = (ov, s, e)
        if best is None:
            unmatched += 1
            continue
        ov, s, e = best
        lab = stratum_of(e - s)
        per_stratum[lab]["n_windows"] += 1
        per_stratum[lab]["content"].append(ov / wlen)
        per_element[(r.genome_id, r.contig_id, s, e)] += 1

    # elements that generated NO positive window at all
    all_el = [(g, c, s, e) for (g, c), ivs in elements.items() for (s, e) in ivs]
    zero = [k for k in all_el if per_element.get(k, 0) == 0]
    zero_by_stratum = defaultdict(int)
    total_by_stratum = defaultdict(int)
    for (g, c, s, e) in all_el:
        total_by_stratum[stratum_of(e - s)] += 1
    for (g, c, s, e) in zero:
        zero_by_stratum[stratum_of(e - s)] += 1

    res = {"split": args.split, "n_windows": int(len(m)), "n_positive_windows": int(len(pos)),
           "n_elements": n_el, "unmatched_positive_windows": unmatched,
           "elements_with_zero_positive_windows": len(zero), "by_stratum": {}}
    order = [lab for _, _, lab in STRATA]
    print(f"\n{'stratum':<10} {'elements':>9} {'zero-win':>9} {'pos windows':>12} "
          f"{'win/elem':>9} {'median content':>15} {'content>=0.9':>13}", flush=True)
    for lab in order:
        d = per_stratum.get(lab, {"n_windows": 0, "content": []})
        ne = total_by_stratum.get(lab, 0)
        c = np.asarray(d["content"]) if d["content"] else np.asarray([0.0])
        row = {"n_elements": ne, "elements_with_zero_windows": zero_by_stratum.get(lab, 0),
               "n_positive_windows": d["n_windows"],
               "windows_per_element": d["n_windows"] / ne if ne else 0.0,
               "median_element_content": float(np.median(c)) if d["content"] else 0.0,
               "frac_windows_content_ge_0.9": float(np.mean(c >= 0.9)) if d["content"] else 0.0}
        res["by_stratum"][lab] = row
        print(f"{lab:<10} {ne:>9} {row['elements_with_zero_windows']:>9} {d['n_windows']:>12} "
              f"{row['windows_per_element']:>9.1f} {row['median_element_content']:>15.3f} "
              f"{row['frac_windows_content_ge_0.9']:>13.3f}", flush=True)

    res["interpretation"] = (
        "median_element_content is the fraction of a positive training window that is "
        "actually element. Where it is low, the model learned that size class from "
        "mostly-host examples, which predicts weak smeared signal -- the observed "
        "0.857 recall / 0.334 IoU at 5-8 kb. elements_with_zero_windows counts elements "
        "that generated no positive training example at all. Either finding argues for "
        "variable-length training windows so a short window can sit inside a short "
        "element; a high-content, well-populated small stratum would instead point at "
        "an architectural cause and make downstream multi-scale the right fix."
    )
    (out / "results.json").write_text(json.dumps(res, indent=2) + "\n")
    print(f"\n[done] wrote {out / 'results.json'}", flush=True)


if __name__ == "__main__":
    main()
