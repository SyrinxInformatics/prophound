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
"""Compare re-run figure outputs (scripts/reproduce_figures.sh) with the recorded ones.

Each file is reported as exact, within 1e-12, or DIFFERS with its first differences.
  python scripts/compare_figures.py <model dir> <reproduce_figures.sh output dir>
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import pandas as pd

TOLERANCE = 1e-12
# (recorded file, relative to the model dir; re-run file, relative to the output dir)
PAIRS = [
    ("analysis/lockbox/results.json", "lockbox/results.json"),
    ("analysis/lockbox/results.elements.csv", "lockbox/results.elements.csv"),
    ("analysis/g5/results.json", "g5/results.json"),
    ("analysis/g5/results.csv", "g5/results.csv"),
    ("analysis/g5b/results.json", "g5b/results.json"),
    ("analysis/g5b/results.csv", "g5b/results.csv"),
    ("analysis/external_benchmarks/lambda/scores.json", "lambda/scores.json"),
    ("analysis/external_benchmarks/lambda/scores.csv", "lambda/scores.csv"),
    ("analysis/external_benchmarks/lambda/chance.json", "lambda/chance.json"),
    ("analysis/offtarget/enrichment_ci.json", "offtarget/enrichment_ci.json"),
    ("analysis/integrity/recall42.json", "integrity/recall42.json"),
    ("analysis/integrity/recall42.csv", "integrity/recall42.csv"),
    ("analysis/integrity/recall18.json", "integrity/recall18.json"),
    ("analysis/integrity/recall18.csv", "integrity/recall18.csv"),
    ("analysis/integrity/recall_phagefree.json", "integrity/recall_phagefree.json"),
    ("analysis/confirm18/seed_vote.csv", "seedvote/confirm18.csv"),
    ("analysis/phagefree/seed_vote.csv", "seedvote/phagefree.csv"),
    ("analysis/keep_rules/results.csv", "keep_rules/results.csv"),
]


def close(a, b) -> bool:
    if not (isinstance(a, (int, float)) and isinstance(b, (int, float))):
        return False
    if a == b or (math.isnan(a) and math.isnan(b)):
        return True
    return abs(a - b) <= TOLERANCE * max(1, abs(a))


def json_diffs(a, b, path: str = "") -> list[str]:
    if isinstance(a, dict) and isinstance(b, dict):
        out = []
        for key in sorted(set(a) | set(b), key=str):
            if key not in a or key not in b:
                out.append(f"{path}.{key}: only in {'re-run' if key not in a else 'recorded'}")
            else:
                out += json_diffs(a[key], b[key], f"{path}.{key}")
        return out
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            return [f"{path}: length {len(a)} vs {len(b)}"]
        return [d for i, (x, y) in enumerate(zip(a, b)) for d in json_diffs(x, y, f"{path}[{i}]")]
    if isinstance(a, float) or isinstance(b, float):
        return [] if close(a, b) else [f"{path}: {a!r} vs {b!r}"]
    return [] if a == b else [f"{path}: {a!r} vs {b!r}"]


def csv_diffs(a: pd.DataFrame, b: pd.DataFrame) -> list[str]:
    if a.shape != b.shape or list(a.columns) != list(b.columns):
        return [f"shape or columns: {a.shape} vs {b.shape}"]
    out = []
    for col in a.columns:
        if a[col].dtype.kind == "f":
            bad = ~(((a[col] - b[col]).abs() <= TOLERANCE) | (a[col].isna() & b[col].isna()))
        else:
            bad = a[col].astype(str) != b[col].astype(str)
        if bad.any():
            out.append(f"column {col}: {int(bad.sum())} rows")
    return out


def compare(recorded: Path, rerun: Path) -> str:
    if recorded.suffix == ".json":
        a, b = json.loads(recorded.read_text()), json.loads(rerun.read_text())
        exact, diffs = a == b, json_diffs(a, b)
    else:
        a, b = (pd.read_csv(p, float_precision="round_trip") for p in (recorded, rerun))
        exact, diffs = a.equals(b), csv_diffs(a, b)
    if exact:
        return "exact"
    if not diffs:
        return f"within {TOLERANCE:g}"
    return f"DIFFERS ({len(diffs)}): " + "; ".join(diffs[:4])


def main() -> None:
    model_dir, out_dir = Path(sys.argv[1]), Path(sys.argv[2])
    for rec, new in PAIRS:
        recorded, rerun = model_dir / rec, out_dir / new
        if not recorded.exists():
            status = "no recorded file"
        elif not rerun.exists():
            status = "MISSING re-run output"
        else:
            status = compare(recorded, rerun)
        print(f"{rec:55s} {status}")


if __name__ == "__main__":
    main()
