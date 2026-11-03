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
from typing import Dict, List, Tuple

import pandas as pd


def build_hysteresis_windows(
    prob_values: List[float],
    seed_threshold: float,
    extend_threshold: float,
) -> List[int]:
    if not prob_values:
        return []

    positive_values = [0] * len(prob_values)
    for seed_index, seed_prob in enumerate(prob_values):
        if seed_prob < seed_threshold or positive_values[seed_index] == 1:
            continue

        left_index = seed_index
        while left_index - 1 >= 0 and prob_values[left_index - 1] >= extend_threshold:
            left_index -= 1

        right_index = seed_index
        while (
            right_index + 1 < len(prob_values)
            and prob_values[right_index + 1] >= extend_threshold
        ):
            right_index += 1

        for region_index in range(left_index, right_index + 1):
            positive_values[region_index] = 1

    return positive_values


def overlap_len(
    left_start: int,
    left_end: int,
    right_start: int,
    right_end: int,
) -> int:
    return max(0, min(left_end, right_end) - max(left_start, right_start))


def iou(
    left_start: int,
    left_end: int,
    right_start: int,
    right_end: int,
) -> float:
    overlap = overlap_len(left_start, left_end, right_start, right_end)
    if overlap <= 0:
        return 0.0
    union = max(left_end, right_end) - min(left_start, right_start)
    return overlap / union if union > 0 else 0.0


def overlapping_regions(
    region_df: pd.DataFrame,
    genome_id: str,
    contig_id: str,
    start: int,
    end: int,
) -> pd.DataFrame:
    same_contig_df = region_df[
        (region_df["genome_id"] == genome_id)
        & (region_df["contig_id"] == contig_id)
    ].copy()
    if same_contig_df.empty:
        return same_contig_df
    return same_contig_df[
        (same_contig_df["end"] > start) & (same_contig_df["start"] < end)
    ].copy()


def classify_predicted_region(
    predicted_region: pd.Series,
    overlapping_true_regions: pd.DataFrame,
) -> Dict[str, object]:
    if overlapping_true_regions.empty:
        return {
            "match_type": "isolated_spurious_pred_region",
            "n_overlapping_true_regions": 0,
            "best_iou": 0.0,
        }

    n_overlaps = len(overlapping_true_regions)
    best_iou = -1.0
    best_true_region = None
    for _, true_region in overlapping_true_regions.iterrows():
        score = iou(
            predicted_region.start,
            predicted_region.end,
            true_region.start,
            true_region.end,
        )
        if score > best_iou:
            best_iou = score
            best_true_region = true_region

    if n_overlaps > 1:
        match_type = "merged_pred_region"
    else:
        left_truncation = max(0, int(predicted_region.start - best_true_region.start))
        right_truncation = max(0, int(best_true_region.end - predicted_region.end))
        left_overextension = max(0, int(best_true_region.start - predicted_region.start))
        right_overextension = max(0, int(predicted_region.end - best_true_region.end))
        if left_truncation > 0 or right_truncation > 0:
            match_type = "truncated_pred_region"
        elif left_overextension > 0 or right_overextension > 0:
            match_type = "overextended_pred_region"
        else:
            match_type = "well_matched_pred_region"

    return {
        "match_type": match_type,
        "n_overlapping_true_regions": int(n_overlaps),
        "best_iou": float(best_iou),
    }


def true_regions_from_windows(predictions_df: pd.DataFrame) -> pd.DataFrame:
    region_rows: List[Dict[str, object]] = []
    for (genome_id, contig_id), contig_df in predictions_df.groupby(
        ["genome_id", "contig_id"],
        sort=False,
    ):
        contig_df = contig_df.sort_values(["start", "end"]).copy()
        true_labels = contig_df["label"].astype(int).to_list()
        window_starts = contig_df["start"].astype(int).to_list()
        window_ends = contig_df["end"].astype(int).to_list()

        region_start = None
        region_end = None
        region_open = False
        for window_index, is_positive in enumerate(true_labels):
            if is_positive and not region_open:
                region_open = True
                region_start = window_starts[window_index]
                region_end = window_ends[window_index]
            elif is_positive and region_open:
                region_end = max(region_end, window_ends[window_index])
            elif (not is_positive) and region_open:
                region_rows.append(
                    {
                        "genome_id": genome_id,
                        "contig_id": contig_id,
                        "start": int(region_start),
                        "end": int(region_end),
                    }
                )
                region_open = False
                region_start = None
                region_end = None

        if region_open:
            region_rows.append(
                {
                    "genome_id": genome_id,
                    "contig_id": contig_id,
                    "start": int(region_start),
                    "end": int(region_end),
                }
            )

    return pd.DataFrame(region_rows)


def region_record(
    region_windows_df: pd.DataFrame,
    genome_id: str,
    contig_id: str,
    seed_threshold: float,
) -> Dict[str, object]:
    prob_values = region_windows_df["prob"].astype(float)
    first_prob = float(prob_values.iloc[0])
    last_prob = float(prob_values.iloc[-1])
    core_mean_prob = (
        float(prob_values.iloc[1:-1].mean()) if len(prob_values) > 2 else float(prob_values.mean())
    )

    return {
        "genome_id": genome_id,
        "contig_id": contig_id,
        "start": int(region_windows_df["start"].min()),
        "end": int(region_windows_df["end"].max()),
        "n_windows": int(len(region_windows_df)),
        "region_len_bp": int(
            region_windows_df["end"].max() - region_windows_df["start"].min()
        ),
        "mean_prob": float(prob_values.mean()),
        "max_prob": float(prob_values.max()),
        "min_prob": float(prob_values.min()),
        "std_prob": float(prob_values.std(ddof=0)),
        "first_prob": first_prob,
        "last_prob": last_prob,
        "edge_mean_prob": float((first_prob + last_prob) / 2.0),
        "core_mean_prob": core_mean_prob,
        "frac_windows_ge_seed": float((prob_values >= seed_threshold).mean()),
        "any_virus_support": int(
            region_windows_df["protein_any_virus_hallmark_gene"].astype(int).max()
        ),
        "virus_support_window_count": int(
            region_windows_df["protein_any_virus_hallmark_gene"].astype(int).sum()
        ),
        "mean_virus_hallmark_gene_count": float(
            region_windows_df["protein_virus_hallmark_gene_count"].astype(float).mean()
        ),
        "max_virus_hallmark_gene_count": int(
            region_windows_df["protein_virus_hallmark_gene_count"].astype(int).max()
        ),
        "boundary_categories": ",".join(
            sorted(set(region_windows_df["boundary_category"].astype(str).tolist()))
        ),
    }


def predicted_regions_from_windows(
    predictions_df: pd.DataFrame,
    seed_threshold: float,
    extend_threshold: float,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    annotated_contigs: List[pd.DataFrame] = []
    predicted_region_rows: List[Dict[str, object]] = []

    for (genome_id, contig_id), contig_df in predictions_df.groupby(
        ["genome_id", "contig_id"],
        sort=False,
    ):
        contig_df = contig_df.sort_values(["start", "end"]).reset_index(drop=True).copy()
        prob_values = contig_df["prob"].astype(float).to_list()
        contig_df["pred_positive_window"] = build_hysteresis_windows(
            prob_values,
            seed_threshold=seed_threshold,
            extend_threshold=extend_threshold,
        )
        annotated_contigs.append(contig_df)

        region_open = False
        start_window_index = None
        end_window_index = None
        for window_index, is_positive in enumerate(
            contig_df["pred_positive_window"].astype(int).to_list()
        ):
            if is_positive and not region_open:
                region_open = True
                start_window_index = window_index
                end_window_index = window_index
            elif is_positive and region_open:
                end_window_index = window_index
            elif (not is_positive) and region_open:
                region_windows_df = contig_df.iloc[
                    start_window_index : end_window_index + 1
                ].copy()
                predicted_region_rows.append(
                    region_record(
                        region_windows_df,
                        genome_id,
                        contig_id,
                        seed_threshold,
                    )
                )
                region_open = False
                start_window_index = None
                end_window_index = None

        if region_open:
            region_windows_df = contig_df.iloc[
                start_window_index : end_window_index + 1
            ].copy()
            predicted_region_rows.append(
                region_record(
                    region_windows_df,
                    genome_id,
                    contig_id,
                    seed_threshold,
                )
            )

    return pd.concat(annotated_contigs, ignore_index=True), pd.DataFrame(predicted_region_rows)


def summarize_cohort(region_df: pd.DataFrame, cohort_name: str) -> Dict[str, object]:
    if region_df.empty:
        return {"cohort": cohort_name, "n_regions": 0}

    return {
        "cohort": cohort_name,
        "n_regions": int(len(region_df)),
        "mean_n_windows": float(region_df["n_windows"].mean()),
        "mean_prob": float(region_df["mean_prob"].mean()),
        "mean_max_prob": float(region_df["max_prob"].mean()),
        "mean_min_prob": float(region_df["min_prob"].mean()),
        "mean_std_prob": float(region_df["std_prob"].mean()),
        "mean_edge_prob": float(region_df["edge_mean_prob"].mean()),
        "mean_core_prob": float(region_df["core_mean_prob"].mean()),
        "mean_frac_windows_ge_seed": float(region_df["frac_windows_ge_seed"].mean()),
        "mean_best_iou": float(region_df["best_iou"].mean())
        if "best_iou" in region_df.columns and not region_df["best_iou"].empty
        else float("nan"),
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Analyze score coherence of predicted regions under the strict hysteresis builder."
    )
    parser.add_argument("--input_csv", required=True)
    parser.add_argument("--out_region_csv", required=True)
    parser.add_argument("--out_summary_json", required=True)
    parser.add_argument("--seed_threshold", type=float, default=0.995)
    parser.add_argument("--extend_threshold", type=float, default=0.99)
    args = parser.parse_args()

    predictions_df = pd.read_csv(args.input_csv)
    required_columns = {
        "genome_id",
        "contig_id",
        "start",
        "end",
        "label",
        "prob",
        "pred_label",
        "protein_any_virus_hallmark_gene",
        "protein_virus_hallmark_gene_count",
        "boundary_category",
    }
    missing_columns = sorted(required_columns - set(predictions_df.columns))
    if missing_columns:
        raise ValueError(f"input missing required columns: {missing_columns}")

    predictions_df, predicted_regions_df = predicted_regions_from_windows(
        predictions_df,
        seed_threshold=args.seed_threshold,
        extend_threshold=args.extend_threshold,
    )
    true_regions_df = true_regions_from_windows(predictions_df)

    annotated_region_rows: List[Dict[str, object]] = []
    for predicted_region in predicted_regions_df.itertuples(index=False):
        overlapping_true_regions = overlapping_regions(
            true_regions_df,
            predicted_region.genome_id,
            predicted_region.contig_id,
            predicted_region.start,
            predicted_region.end,
        )
        predicted_region_dict = predicted_region._asdict()
        predicted_region_dict.update(
            classify_predicted_region(
                pd.Series(predicted_region_dict),
                overlapping_true_regions,
            )
        )
        predicted_region_dict["any_true_overlap"] = int(
            not overlapping_true_regions.empty
        )
        predicted_region_dict["cohort_fp_supported"] = int(
            predicted_region_dict["match_type"] == "isolated_spurious_pred_region"
            and predicted_region_dict["any_virus_support"] == 1
        )
        predicted_region_dict["cohort_tp_unsupported"] = int(
            predicted_region_dict["any_true_overlap"] == 1
            and predicted_region_dict["any_virus_support"] == 0
        )
        annotated_region_rows.append(predicted_region_dict)

    regions_df = pd.DataFrame(annotated_region_rows)
    cohort_frames = {
        "fp_supported_regions": regions_df[
            (regions_df["match_type"] == "isolated_spurious_pred_region")
            & (regions_df["any_virus_support"] == 1)
        ].copy(),
        "tp_unsupported_regions": regions_df[
            (regions_df["any_true_overlap"] == 1)
            & (regions_df["any_virus_support"] == 0)
        ].copy(),
        "tp_supported_regions": regions_df[
            (regions_df["any_true_overlap"] == 1)
            & (regions_df["any_virus_support"] == 1)
        ].copy(),
        "fp_unsupported_regions": regions_df[
            (regions_df["match_type"] == "isolated_spurious_pred_region")
            & (regions_df["any_virus_support"] == 0)
        ].copy(),
    }

    summary = {
        "seed_threshold": args.seed_threshold,
        "extend_threshold": args.extend_threshold,
        "n_pred_regions": int(len(regions_df)),
        "cohorts": [
            summarize_cohort(region_df, cohort_name)
            for cohort_name, region_df in cohort_frames.items()
        ],
    }

    out_region_path = Path(args.out_region_csv)
    out_summary_path = Path(args.out_summary_json)
    out_region_path.parent.mkdir(parents=True, exist_ok=True)
    regions_df.to_csv(out_region_path, index=False)
    out_summary_path.write_text(json.dumps(summary, indent=2) + "\n")

    print(f"[done] wrote {out_region_path}")
    print(f"[done] wrote {out_summary_path}")


if __name__ == "__main__":
    main()
