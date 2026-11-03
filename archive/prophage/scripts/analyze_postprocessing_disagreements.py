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
import json
from pathlib import Path
from typing import Dict, List

import pandas as pd


REGION_KEYS = ["genome_id", "contig_id", "start", "end"]


def build_transition_table(
    baseline_df: pd.DataFrame, candidate_df: pd.DataFrame
) -> pd.DataFrame:
    merged_df = baseline_df.merge(
        candidate_df,
        on=REGION_KEYS,
        how="outer",
        suffixes=("_baseline", "_candidate"),
        indicator=True,
    )
    merged_df["match_type_baseline"] = merged_df["match_type_baseline"].fillna("missing")
    merged_df["match_type_candidate"] = merged_df["match_type_candidate"].fillna("missing")
    transition_df = (
        merged_df.groupby(["match_type_baseline", "match_type_candidate"], dropna=False)
        .size()
        .reset_index(name="count")
        .sort_values("count", ascending=False)
        .reset_index(drop=True)
    )
    return merged_df, transition_df


def transition_examples(
    merged_df: pd.DataFrame,
    baseline_type: str,
    candidate_type: str,
    limit: int,
) -> List[Dict[str, object]]:
    subset = merged_df[
        (merged_df["match_type_baseline"] == baseline_type)
        & (merged_df["match_type_candidate"] == candidate_type)
    ].copy()
    if subset.empty:
        return []

    sort_columns = [
        col
        for col in [
            "best_iou_candidate",
            "best_iou_baseline",
            "left_truncation_bp_candidate",
            "right_truncation_bp_candidate",
        ]
        if col in subset.columns
    ]
    if sort_columns:
        ascending = [False, False, True, True][: len(sort_columns)]
        subset = subset.sort_values(sort_columns, ascending=ascending)

    example_columns = [
        "genome_id",
        "contig_id",
        "start",
        "end",
        "match_type_baseline",
        "match_type_candidate",
        "n_overlapping_pred_regions_baseline",
        "n_overlapping_pred_regions_candidate",
        "best_iou_baseline",
        "best_iou_candidate",
        "left_truncation_bp_baseline",
        "left_truncation_bp_candidate",
        "right_truncation_bp_baseline",
        "right_truncation_bp_candidate",
    ]
    keep_columns = [col for col in example_columns if col in subset.columns]
    return subset[keep_columns].head(limit).to_dict(orient="records")


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Compare baseline and candidate true-region classifications."
    )
    ap.add_argument("--baseline_true_csv", required=True)
    ap.add_argument("--candidate_true_csv", required=True)
    ap.add_argument("--out_transition_csv", required=True)
    ap.add_argument("--summary_json", required=True)
    ap.add_argument("--example_limit", type=int, default=10)
    args = ap.parse_args()

    baseline_df = pd.read_csv(args.baseline_true_csv)
    candidate_df = pd.read_csv(args.candidate_true_csv)

    merged_df, transition_df = build_transition_table(baseline_df, candidate_df)

    out_transition_csv = Path(args.out_transition_csv)
    out_transition_csv.parent.mkdir(parents=True, exist_ok=True)
    transition_df.to_csv(out_transition_csv, index=False)

    summary = {
        "baseline_true_csv": args.baseline_true_csv,
        "candidate_true_csv": args.candidate_true_csv,
        "n_baseline_rows": int(len(baseline_df)),
        "n_candidate_rows": int(len(candidate_df)),
        "transition_counts": transition_df.to_dict(orient="records"),
        "examples": {
            "fragmented_to_well_matched": transition_examples(
                merged_df,
                "fragmented_true_region",
                "well_matched_true_region",
                args.example_limit,
            ),
            "fragmented_to_truncated": transition_examples(
                merged_df,
                "fragmented_true_region",
                "truncated_true_region",
                args.example_limit,
            ),
            "truncated_to_well_matched": transition_examples(
                merged_df,
                "truncated_true_region",
                "well_matched_true_region",
                args.example_limit,
            ),
            "well_matched_to_truncated": transition_examples(
                merged_df,
                "well_matched_true_region",
                "truncated_true_region",
                args.example_limit,
            ),
        },
    }

    summary_path = Path(args.summary_json)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")

    print(f"[done] wrote transition csv -> {out_transition_csv}")
    print(f"[done] wrote summary json -> {summary_path}")


if __name__ == "__main__":
    main()
