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


"""Composition-only baseline: can GC/k-mer statistics alone match the CNN?

docs/scientific-validity.md item 1c -- the cheapest falsification available.
The CNN reaches window-level val/test AUPRC 0.63-0.65. If a logistic regression
on plain composition features lands near that, the CNN's advantage is not
phage-specific signal and the project's framing has to change. If it lands well
below, that is the first positive evidence that the model uses something beyond
nucleotide composition.

Features are deliberately dumb and database-free: GC fraction, GC skew, AT skew,
Shannon entropy, and canonical k-mer frequencies for k=1..K (canonical because
the model is reverse-complement invariant, so the baseline should be too).
No ORF density -- that needs gene calling; noted as a limitation rather than
faked.

Runs on the training manifest with the project's own splits, so the comparison
is like-for-like. Windows are subsampled (stratified by label) because the full
manifest is ~10.5M windows; the default 150k/50k is ample for separating "near
0.65" from "near 0.2".

Usage (see composition_baseline.slurm):
    python composition_baseline.py --manifest data/manifest_windows.filtered.csv \
        --splits_dir splits --out_dir analysis/composition_baseline \
        --n_train 150000 --n_test 50000 --k 4
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

# Imported at module scope on purpose. The HPCC login node cannot load the
# miniconda module at all ("available on the job nodes" -- docs/environment.md),
# so whether sklearn exists in the env is only discoverable inside a job. Import
# it here and the job dies in seconds with a clear error instead of after ~an
# hour of featurization.
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.preprocessing import StandardScaler

BASES = np.full(256, 255, dtype=np.uint8)
for i, b in enumerate(b"ACGT"):
    BASES[b] = i
    BASES[b + 32] = i  # lowercase


def canonical_kmer_map(k: int) -> np.ndarray:
    """Map every k-mer index to its canonical (min of itself and its revcomp)."""
    n = 4 ** k
    idx = np.arange(n, dtype=np.int64)
    rc = np.zeros(n, dtype=np.int64)
    tmp = idx.copy()
    for _ in range(k):
        rc = (rc << 2) | (3 - (tmp & 3))
        tmp >>= 2
    return np.minimum(idx, rc)


def window_features(seq: str, k: int, canon: dict) -> np.ndarray:
    """GC/skew/entropy + canonical k-mer frequencies for k=1..k."""
    arr = BASES[np.frombuffer(seq.encode("ascii", "ignore"), dtype=np.uint8)]
    arr = arr[arr < 4]
    if arr.size < k:
        return None
    counts = np.bincount(arr, minlength=4).astype(np.float64)
    total = counts.sum()
    a, c, g, t = counts / total
    gc = g + c
    gc_skew = (g - c) / (g + c) if (g + c) > 0 else 0.0
    at_skew = (a - t) / (a + t) if (a + t) > 0 else 0.0
    p = counts[counts > 0] / total
    entropy = float(-(p * np.log2(p)).sum())
    feats = [gc, gc_skew, at_skew, entropy]

    for kk in range(1, k + 1):
        codes = arr[: arr.size - kk + 1].astype(np.int64)
        for shift in range(1, kk):
            codes = (codes << 2) | arr[shift : arr.size - kk + 1 + shift]
        cmap = canon[kk]
        binned = np.bincount(cmap[codes], minlength=4 ** kk).astype(np.float64)
        keep = binned[np.unique(cmap)]
        feats.extend(keep / max(1.0, keep.sum()))
    return np.asarray(feats, dtype=np.float32)


def load_split_ids(splits_dir: Path, name: str) -> set:
    return {l.strip() for l in open(splits_dir / f"{name}.txt") if l.strip()}


def build_matrix(df: pd.DataFrame, k: int, canon: dict, label: str):
    from pyfaidx import Fasta

    cache: dict = {}
    X, y, skipped = [], [], 0
    for i, row in enumerate(df.itertuples()):
        fa = cache.get(row.fasta_path)
        if fa is None:
            if len(cache) > 8:
                cache.pop(next(iter(cache)))
            fa = Fasta(row.fasta_path, as_raw=True, sequence_always_upper=True)
            cache[row.fasta_path] = fa
        try:
            seq = fa[row.contig_id][int(row.start) : int(row.end)]
        except Exception:
            skipped += 1
            continue
        f = window_features(seq, k, canon)
        if f is None:
            skipped += 1
            continue
        X.append(f)
        y.append(int(row.label))
        if (i + 1) % 25000 == 0:
            print(f"[{label}] {i + 1}/{len(df)} windows featurized", flush=True)
    print(f"[{label}] done: {len(X)} windows, {skipped} skipped", flush=True)
    return np.vstack(X), np.asarray(y)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--splits_dir", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--n_train", type=int, default=150000)
    ap.add_argument("--n_test", type=int, default=50000)
    ap.add_argument("--k", type=int, default=4)
    ap.add_argument("--seed", type=int, default=13)
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    canon = {kk: canonical_kmer_map(kk) for kk in range(1, args.k + 1)}

    print("[data] reading manifest", flush=True)
    m = pd.read_csv(args.manifest, usecols=["genome_id", "contig_id", "start", "end", "label", "neg_type", "fasta_path"])
    splits = Path(args.splits_dir)
    train_ids, test_ids = load_split_ids(splits, "train"), load_split_ids(splits, "test")

    def subsample(ids: set, n: int, name: str) -> pd.DataFrame:
        sub = m[m.genome_id.isin(ids)]
        pos, neg = sub[sub.label == 1], sub[sub.label == 0]
        # keep the manifest's own class ratio rather than balancing: AUPRC is
        # prevalence-dependent, and the CNN's 0.63-0.65 was measured at the
        # natural ~1.4% positive rate.
        n_pos = max(1, int(round(n * len(pos) / max(1, len(sub)))))
        n_neg = n - n_pos
        take = pd.concat([
            pos.sample(min(n_pos, len(pos)), random_state=args.seed),
            neg.sample(min(n_neg, len(neg)), random_state=args.seed),
        ]).sample(frac=1.0, random_state=args.seed)
        print(f"[{name}] {len(take)} windows ({take.label.mean() * 100:.2f}% positive) from {len(ids)} genomes", flush=True)
        return take

    tr = subsample(train_ids, args.n_train, "train")
    te = subsample(test_ids, args.n_test, "test")

    Xtr, ytr = build_matrix(tr, args.k, canon, "train")
    Xte, yte = build_matrix(te, args.k, canon, "test")

    results = {"n_train": int(len(ytr)), "n_test": int(len(yte)),
               "test_positive_rate": float(yte.mean()), "k": args.k,
               "n_features": int(Xtr.shape[1])}

    scaler = StandardScaler().fit(Xtr)
    lr = LogisticRegression(max_iter=2000, class_weight="balanced")
    lr.fit(scaler.transform(Xtr), ytr)
    p = lr.predict_proba(scaler.transform(Xte))[:, 1]
    results["logreg"] = {"auprc": float(average_precision_score(yte, p)),
                         "auroc": float(roc_auc_score(yte, p))}

    gb = HistGradientBoostingClassifier(max_iter=300, random_state=args.seed)
    gb.fit(Xtr, ytr)
    p = gb.predict_proba(Xte)[:, 1]
    results["hist_gbm"] = {"auprc": float(average_precision_score(yte, p)),
                           "auroc": float(roc_auc_score(yte, p))}

    # the number this whole script exists to be compared against
    results["cnn_reference_auprc"] = "0.63-0.65 (docs/model-reengineering.md)"
    results["interpretation"] = (
        "If either baseline AUPRC approaches 0.63-0.65, composition alone explains "
        "the CNN and docs/scientific-validity.md claim 2 fails. If both land well "
        "below, the CNN uses signal beyond nucleotide composition."
    )
    (out / "results.json").write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results, indent=2), flush=True)


if __name__ == "__main__":
    main()
