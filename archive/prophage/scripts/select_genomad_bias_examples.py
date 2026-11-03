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
from typing import Callable, Dict, List

import pandas as pd


def filter_cohort(
    data_frame: pd.DataFrame,
    cohort_name: str,
    predicate: Callable[[pd.DataFrame], pd.Series],
    sort_by: List[str],
    ascending: List[bool],
    top_n: int,
) -> pd.DataFrame:
    cohort_df = data_frame[predicate(data_frame)].copy()
    if cohort_df.empty:
        return cohort_df
    cohort_df = cohort_df.sort_values(sort_by, ascending=ascending).head(top_n).copy()
    cohort_df.insert(0, "cohort", cohort_name)
    return cohort_df


def build_cohorts(data_frame: pd.DataFrame, top_n: int) -> Dict[str, pd.DataFrame]:
    adjacent_categories = {"neg_adjacent_1", "neg_adjacent_2_3", "neg_adjacent_4_10"}
    cohorts: Dict[str, pd.DataFrame] = {}

    cohorts["tp_unsupported_highconf"] = filter_cohort(
        data_frame,
        "tp_unsupported_highconf",
        lambda df: (
            (df["label"].astype(int) == 1)
            & (df["pred_label"].astype(int) == 1)
            & (df["protein_any_virus_hallmark_gene"].astype(int) == 0)
        ),
        sort_by=["prob", "protein_virus_hallmark_gene_count"],
        ascending=[False, False],
        top_n=top_n,
    )

    cohorts["fn_supported_nearmiss"] = filter_cohort(
        data_frame,
        "fn_supported_nearmiss",
        lambda df: (
            (df["label"].astype(int) == 1)
            & (df["pred_label"].astype(int) == 0)
            & (df["protein_any_virus_hallmark_gene"].astype(int) == 1)
        ),
        sort_by=["prob", "protein_virus_hallmark_gene_count"],
        ascending=[False, False],
        top_n=top_n,
    )

    cohorts["fp_supported_adjacent"] = filter_cohort(
        data_frame,
        "fp_supported_adjacent",
        lambda df: (
            (df["label"].astype(int) == 0)
            & (df["pred_label"].astype(int) == 1)
            & (df["protein_any_virus_hallmark_gene"].astype(int) == 1)
            & (df["boundary_category"].isin(adjacent_categories))
        ),
        sort_by=["prob", "protein_virus_hallmark_gene_count"],
        ascending=[False, False],
        top_n=top_n,
    )

    cohorts["fp_supported_far"] = filter_cohort(
        data_frame,
        "fp_supported_far",
        lambda df: (
            (df["label"].astype(int) == 0)
            & (df["pred_label"].astype(int) == 1)
            & (df["protein_any_virus_hallmark_gene"].astype(int) == 1)
            & (df["boundary_category"] == "neg_far_gt10")
        ),
        sort_by=["prob", "protein_virus_hallmark_gene_count"],
        ascending=[False, False],
        top_n=top_n,
    )

    cohorts["tp_supported_edge"] = filter_cohort(
        data_frame,
        "tp_supported_edge",
        lambda df: (
            (df["label"].astype(int) == 1)
            & (df["pred_label"].astype(int) == 1)
            & (df["protein_any_virus_hallmark_gene"].astype(int) == 1)
            & (df["boundary_category"] == "positive_edge")
        ),
        sort_by=["prob", "protein_virus_hallmark_gene_count"],
        ascending=[False, False],
        top_n=top_n,
    )

    cohorts["tp_unsupported_edge"] = filter_cohort(
        data_frame,
        "tp_unsupported_edge",
        lambda df: (
            (df["label"].astype(int) == 1)
            & (df["pred_label"].astype(int) == 1)
            & (df["protein_any_virus_hallmark_gene"].astype(int) == 0)
            & (df["boundary_category"] == "positive_edge")
        ),
        sort_by=["prob", "protein_virus_hallmark_gene_count"],
        ascending=[False, False],
        top_n=top_n,
    )

    return cohorts


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Select representative cases for geNomad-bias inspection from joined prediction tables."
    )
    parser.add_argument("--input_csv", required=True)
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--top_n", type=int, default=25)
    args = parser.parse_args()

    input_df = pd.read_csv(args.input_csv)
    cohorts = build_cohorts(input_df, top_n=args.top_n)

    output_dir = Path(args.out_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    summary_rows = []
    for cohort_name, cohort_df in cohorts.items():
        out_csv = output_dir / f"{cohort_name}.csv"
        cohort_df.to_csv(out_csv, index=False)
        summary_rows.append({"cohort": cohort_name, "n_rows": int(len(cohort_df)), "path": str(out_csv)})

    summary_path = output_dir / "cohort_summary.json"
    summary_path.write_text(json.dumps(summary_rows, indent=2) + "\n")
    print(f"[done] wrote {summary_path}")
    for row in summary_rows:
        print(f"[done] {row['cohort']} -> {row['n_rows']} rows")


if __name__ == "__main__":
    main()
