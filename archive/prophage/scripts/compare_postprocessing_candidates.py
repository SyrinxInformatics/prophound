#!/usr/bin/env python3
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

import argparse
import csv
import json
from pathlib import Path
from statistics import mean
from typing import Dict, List, Tuple


KEY_METRICS = [
    "n_pred_regions",
    "well_matched_true_region",
    "fragmented_true_region",
    "truncated_true_region",
    "missed_true_region",
    "isolated_spurious_pred_region",
]


def extract_run_name(path: Path) -> str:
    parts = list(path.parts)
    if "runs" in parts:
        idx = parts.index("runs")
        if idx + 1 < len(parts):
            return parts[idx + 1]
    return path.parent.parent.name


def extract_variant(path: Path) -> str:
    name = path.stem
    prefix = "region_boundary_summary_"
    if name.startswith(prefix):
        return name[len(prefix) :]
    return name


def summary_to_row(path: Path) -> Dict[str, object]:
    data = json.loads(path.read_text())
    row: Dict[str, object] = {
        "run_name": extract_run_name(path),
        "variant": extract_variant(path),
        "summary_json": str(path),
        "threshold": data.get("threshold"),
        "max_gap_windows": data.get("max_gap_windows", 0),
        "min_run_windows": data.get("min_run_windows", 1),
        "min_pred_region_windows": data.get("min_pred_region_windows", 1),
        "min_pred_region_bp": data.get("min_pred_region_bp", 0),
        "n_true_regions": data.get("n_true_regions", 0),
        "n_pred_regions": data.get("n_pred_regions", 0),
    }

    for item in data.get("true_region_summary", []):
        row[item["match_type"]] = item["count"]
    for item in data.get("pred_region_summary", []):
        row[item["match_type"]] = item["count"]
    return row


def build_delta_rows(rows: List[Dict[str, object]], baseline_variant: str) -> List[Dict[str, object]]:
    by_run_variant: Dict[Tuple[str, str], Dict[str, object]] = {}
    for row in rows:
        by_run_variant[(str(row["run_name"]), str(row["variant"]))] = row

    delta_rows: List[Dict[str, object]] = []
    run_names = sorted({str(row["run_name"]) for row in rows})
    variants = sorted({str(row["variant"]) for row in rows if str(row["variant"]) != baseline_variant})

    for run_name in run_names:
        baseline = by_run_variant.get((run_name, baseline_variant))
        if baseline is None:
            continue
        for variant in variants:
            candidate = by_run_variant.get((run_name, variant))
            if candidate is None:
                continue
            delta_row: Dict[str, object] = {"run_name": run_name, "variant": variant}
            for metric in KEY_METRICS:
                base_val = int(baseline.get(metric, 0))
                cand_val = int(candidate.get(metric, 0))
                delta_row[f"{metric}_delta"] = cand_val - base_val
            delta_rows.append(delta_row)
    return delta_rows


def build_rollup(delta_rows: List[Dict[str, object]]) -> Dict[str, Dict[str, float]]:
    by_variant: Dict[str, Dict[str, List[float]]] = {}
    for row in delta_rows:
        variant = str(row["variant"])
        stats = by_variant.setdefault(variant, {f"{metric}_delta": [] for metric in KEY_METRICS})
        for metric in KEY_METRICS:
            stats[f"{metric}_delta"].append(float(row[f"{metric}_delta"]))

    out: Dict[str, Dict[str, float]] = {}
    for variant, stats in by_variant.items():
        out[variant] = {}
        for metric_name, values in stats.items():
            out[variant][metric_name] = mean(values) if values else 0.0
    return out


def write_csv(path: Path, rows: List[Dict[str, object]]) -> None:
    if not rows:
        return
    fieldnames = sorted({key for row in rows for key in row.keys()})
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Compare region-boundary postprocessing summaries across runs."
    )
    ap.add_argument(
        "--summary_json",
        action="append",
        required=True,
        help="Path to a region_boundary_summary_*.json file. Repeat for each input.",
    )
    ap.add_argument(
        "--baseline_variant",
        default="baseline",
        help="Variant name to treat as baseline. Default: baseline.",
    )
    ap.add_argument("--out_csv", required=True, help="Flat comparison CSV output.")
    ap.add_argument(
        "--out_delta_csv",
        required=True,
        help="Per-run delta CSV relative to the baseline variant.",
    )
    ap.add_argument(
        "--out_rollup_json",
        required=True,
        help="Aggregate mean delta summary across runs.",
    )
    args = ap.parse_args()

    summary_paths = [Path(p) for p in args.summary_json]
    rows = [summary_to_row(path) for path in summary_paths]
    delta_rows = build_delta_rows(rows, args.baseline_variant)
    rollup = build_rollup(delta_rows)

    out_csv = Path(args.out_csv)
    out_delta_csv = Path(args.out_delta_csv)
    out_rollup_json = Path(args.out_rollup_json)
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    write_csv(out_csv, rows)
    write_csv(out_delta_csv, delta_rows)
    out_rollup_json.write_text(json.dumps(rollup, indent=2) + "\n")

    print(f"[done] wrote {out_csv}")
    print(f"[done] wrote {out_delta_csv}")
    print(f"[done] wrote {out_rollup_json}")


if __name__ == "__main__":
    main()
