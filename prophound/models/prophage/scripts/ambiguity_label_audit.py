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


"""Why does the model treat ambiguity codes as a hard negative? (docs/scientific-validity.md item 1a)

Local work on 2026-09-22 established that ambiguity codes (N and IUPAC codes,
all of which route to the encoder's N channel) collapse the model's output:
in EDL933, **no** window free of ambiguity codes is saturated, while a window
with >=5 of them is saturated 77% of the time, and cleaning 2-7 such bases
turns p=0.000 into p=1.000. See `scripts/rc_n_resilience.py`.

This job asks the training-side question that decides whether an inference-time
fix is sufficient or a retrain is needed: **what did the model see?** If
positive training windows are essentially ambiguity-free while negatives carry
them, the model had every reason to learn "ambiguity code => not a prophage",
and the behaviour is a training-data artifact rather than a modelling choice.

Outputs `analysis/ambiguity_labels/results.json`:
  * fraction of positive vs negative training windows containing >=1 ambiguity
    code, and the distribution of counts per window
  * the same split by `neg_type` (near vs far negatives)
  * how many positive windows would have been the model's only exposure to an
    ambiguity code inside a real prophage

Reads the FASTAs directly (the manifest stores coordinates, not sequence), so
it is I/O bound; it opens each genome once and closes it (see the EMFILE
incident in docs/environment.md).
"""

from __future__ import annotations

import argparse
import json
import time
from collections import Counter
from pathlib import Path

import pandas as pd
from pyfaidx import Fasta, IndexNotFoundError

AMBIG = None  # set in main from the complement of ACGT


def open_fasta(path: str, retries: int = 5, wait_s: float = 2.0) -> Fasta:
    """pyfaidx open tolerant of a concurrent .fai rewrite (docs/environment.md)."""
    for attempt in range(retries):
        try:
            return Fasta(path, as_raw=True, sequence_always_upper=True)
        except IndexNotFoundError as e:
            if attempt == retries - 1:
                raise
            print(f"[fasta] {path}: {e}; retry {attempt + 1}/{retries - 1}", flush=True)
            time.sleep(wait_s * (attempt + 1))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--splits_dir", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--n_sample", type=int, default=400000,
                    help="windows to sample per class; positives are ~1.4%% of the manifest so they are taken whole")
    ap.add_argument("--seed", type=int, default=13)
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    print("[data] reading manifest", flush=True)
    m = pd.read_csv(args.manifest,
                    usecols=["genome_id", "contig_id", "start", "end", "label", "neg_type", "fasta_path"])
    train = {l.strip() for l in open(Path(args.splits_dir) / "train.txt") if l.strip()}
    m = m[m.genome_id.isin(train)]
    pos = m[m.label == 1]
    neg = m[m.label == 0]
    print(f"[data] train windows: {len(m)} ({len(pos)} positive, {len(neg)} negative)", flush=True)
    if len(neg) > args.n_sample:
        neg = neg.sample(args.n_sample, random_state=args.seed)
    if len(pos) > args.n_sample:
        pos = pos.sample(args.n_sample, random_state=args.seed)
    sub = pd.concat([pos, neg]).sort_values("fasta_path", kind="stable")

    counts: dict = {}
    fa = None
    cur = None
    skipped = 0
    try:
        for i, row in enumerate(sub.itertuples()):
            if row.fasta_path != cur:
                if fa is not None:
                    fa.close()
                fa = open_fasta(row.fasta_path)
                cur = row.fasta_path
            try:
                seq = fa[row.contig_id][int(row.start) : int(row.end)]
            except Exception:
                skipped += 1
                continue
            n_amb = sum(1 for b in seq if b not in "ACGT")
            key = (int(row.label), str(row.neg_type))
            c = counts.setdefault(key, Counter())
            c["n"] += 1
            c["with_ambig"] += 1 if n_amb else 0
            c["total_ambig"] += n_amb
            c[f"bucket_{min(n_amb, 10) if n_amb < 10 else '10plus'}"] += 1
            if (i + 1) % 50000 == 0:
                print(f"[scan] {i + 1}/{len(sub)}", flush=True)
    finally:
        if fa is not None:
            fa.close()

    res = {"n_train_windows": int(len(m)), "n_scanned": int(len(sub)), "n_skipped": skipped, "by_class": {}}
    for (label, neg_type), c in sorted(counts.items()):
        key = f"label{label}_{neg_type}"
        res["by_class"][key] = {
            "n": int(c["n"]),
            "frac_with_ambiguity": c["with_ambig"] / max(1, c["n"]),
            "mean_ambiguity_per_window": c["total_ambig"] / max(1, c["n"]),
            "buckets": {k: int(v) for k, v in c.items() if k.startswith("bucket_")},
        }
    pos_keys = [k for k in res["by_class"] if k.startswith("label1")]
    neg_keys = [k for k in res["by_class"] if k.startswith("label0")]
    pn = sum(res["by_class"][k]["n"] for k in pos_keys)
    pa = sum(res["by_class"][k]["frac_with_ambiguity"] * res["by_class"][k]["n"] for k in pos_keys)
    nn = sum(res["by_class"][k]["n"] for k in neg_keys)
    na = sum(res["by_class"][k]["frac_with_ambiguity"] * res["by_class"][k]["n"] for k in neg_keys)
    res["summary"] = {
        "positive_frac_with_ambiguity": pa / max(1, pn),
        "negative_frac_with_ambiguity": na / max(1, nn),
        "positive_windows_with_ambiguity": int(round(pa)),
        "interpretation": (
            "If positives are far cleaner than negatives, the model learned "
            "'ambiguity code => not a prophage' from the data and an "
            "inference-time fix (impute, or encode ambiguity as 0.25 across "
            "ACGT) addresses the symptom while N-augmented retraining is the "
            "durable fix. If both classes carry ambiguity codes at similar "
            "rates, the collapse is a modelling artifact and retraining is "
            "the only route."
        ),
    }
    (out / "results.json").write_text(json.dumps(res, indent=2) + "\n")
    print(json.dumps(res["summary"], indent=2), flush=True)


if __name__ == "__main__":
    main()
