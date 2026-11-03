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

from __future__ import annotations

import argparse
import os
import random
from typing import List, Tuple


def stem(filename: str) -> str:
    base = os.path.basename(filename)
    # remove common extensions
    for ext in [".bed", ".tsv", ".txt", ".fasta", ".fa", ".fna", ".fas"]:
        if base.endswith(ext):
            base = base[: -len(ext)]
    return base


def bed_has_positive(bed_path: str) -> bool:
    with open(bed_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            seg_type = parts[3]
            if seg_type == "positive":
                return True
    return False


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--beds_dir",
        required=True,
        help="Directory containing BED maps with far_neg/near_neg/positive",
    )
    ap.add_argument(
        "--out_dir", default="splits", help="Output directory for train/val/test.txt"
    )
    ap.add_argument("--seed", type=int, default=13)
    ap.add_argument("--train_frac", type=float, default=0.70)
    ap.add_argument("--val_frac", type=float, default=0.15)
    ap.add_argument("--test_frac", type=float, default=0.15)
    args = ap.parse_args()

    if abs((args.train_frac + args.val_frac + args.test_frac) - 1.0) > 1e-6:
        raise SystemExit("Fractions must sum to 1.0")

    bed_files = [
        os.path.join(args.beds_dir, f)
        for f in os.listdir(args.beds_dir)
        if f.endswith((".bed", ".tsv", ".txt"))
    ]
    if not bed_files:
        raise SystemExit(f"No BED files found in {args.beds_dir}")

    pos_genomes: List[str] = []
    negonly_genomes: List[str] = []

    for bed in sorted(bed_files):
        gid = stem(bed)
        if bed_has_positive(bed):
            pos_genomes.append(gid)
        else:
            negonly_genomes.append(gid)

    rng = random.Random(args.seed)
    rng.shuffle(pos_genomes)
    rng.shuffle(negonly_genomes)

    def split_list(xs: List[str]) -> Tuple[List[str], List[str], List[str]]:
        n = len(xs)
        n_train = int(round(n * args.train_frac))
        n_val = int(round(n * args.val_frac))
        # ensure sum <= n; remainder goes to test
        n_train = min(n_train, n)
        n_val = min(n_val, max(0, n - n_train))
        train = xs[:n_train]
        val = xs[n_train : n_train + n_val]
        test = xs[n_train + n_val :]
        return train, val, test

    pos_tr, pos_va, pos_te = split_list(pos_genomes)
    neg_tr, neg_va, neg_te = split_list(negonly_genomes)

    train = pos_tr + neg_tr
    val = pos_va + neg_va
    test = pos_te + neg_te

    # Shuffle within each split to avoid grouping
    rng.shuffle(train)
    rng.shuffle(val)
    rng.shuffle(test)

    os.makedirs(args.out_dir, exist_ok=True)

    def write(path: str, ids: List[str]) -> None:
        with open(path, "w") as f:
            for gid in ids:
                f.write(gid + "\n")

    write(os.path.join(args.out_dir, "train.txt"), train)
    write(os.path.join(args.out_dir, "val.txt"), val)
    write(os.path.join(args.out_dir, "test.txt"), test)

    def pct(a: int, b: int) -> str:
        return "0.0%" if b == 0 else f"{100.0 * a / b:.1f}%"

    print("Split summary:")
    print(f"  total genomes: {len(pos_genomes) + len(negonly_genomes)}")
    print(f"  genomes w/ positives: {len(pos_genomes)}")
    print(f"  negative-only genomes: {len(negonly_genomes)}")
    print("  train:", len(train), f"(pos {len(pos_tr)} / neg-only {len(neg_tr)})")
    print("  val:  ", len(val), f"(pos {len(pos_va)} / neg-only {len(neg_va)})")
    print("  test: ", len(test), f"(pos {len(pos_te)} / neg-only {len(neg_te)})")

    # Leakage check
    dup = set(train) & set(val) | set(train) & set(test) | set(val) & set(test)
    if dup:
        raise SystemExit(
            f"ERROR: genome IDs overlap between splits: {sorted(list(dup))[:10]}"
        )
    print("  leakage check: OK (no overlapping genome IDs)")


if __name__ == "__main__":
    main()

