import argparse
import json
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd


GROUP_KEYS = ["genome_id", "contig_id"]


def nearest_positive_window_distance(labels: np.ndarray) -> np.ndarray:
    positive_indices = np.flatnonzero(labels == 1)
    boundary_distance = np.full(len(labels), -1, dtype=int)
    if len(positive_indices) == 0:
        return boundary_distance

    window_indices = np.arange(len(labels))
    insertion_points = np.searchsorted(positive_indices, window_indices)

    has_left_positive = insertion_points > 0
    has_right_positive = insertion_points < len(positive_indices)

    left_dist = np.full(len(labels), np.inf)
    right_dist = np.full(len(labels), np.inf)

    left_dist[has_left_positive] = (
        window_indices[has_left_positive]
        - positive_indices[insertion_points[has_left_positive] - 1]
    )
    right_dist[has_right_positive] = (
        positive_indices[insertion_points[has_right_positive]]
        - window_indices[has_right_positive]
    )

    boundary_distance = np.minimum(left_dist, right_dist).astype(int)
    return boundary_distance


def positive_block_positions(labels: np.ndarray) -> Dict[str, np.ndarray]:
    n = len(labels)
    edge = np.zeros(n, dtype=bool)
    core = np.zeros(n, dtype=bool)
    if n == 0:
        return {"positive_edge": edge, "positive_core": core}

    pos = labels == 1
    for i in range(n):
        if not pos[i]:
            continue
        prev_pos = i > 0 and pos[i - 1]
        next_pos = i < n - 1 and pos[i + 1]
        if prev_pos and next_pos:
            core[i] = True
        else:
            edge[i] = True
    return {"positive_edge": edge, "positive_core": core}


def assign_boundary_category(predictions_df: pd.DataFrame) -> pd.DataFrame:
    annotated_parts: List[pd.DataFrame] = []

    for _, group_df in predictions_df.groupby(GROUP_KEYS, sort=False):
        group_df = group_df.sort_values(["start", "end"]).copy()
        true_labels = group_df["label"].astype(int).to_numpy()
        predicted_labels = group_df["pred_label"].astype(int).to_numpy()
        boundary_distance = nearest_positive_window_distance(true_labels)
        positive_positions = positive_block_positions(true_labels)

        boundary_category = np.full(len(group_df), "", dtype=object)
        boundary_category[positive_positions["positive_edge"]] = "positive_edge"
        boundary_category[positive_positions["positive_core"]] = "positive_core"

        negative_mask = true_labels == 0
        boundary_category[negative_mask & (boundary_distance == 1)] = "neg_adjacent_1"
        boundary_category[
            negative_mask & np.isin(boundary_distance, [2, 3])
        ] = "neg_adjacent_2_3"
        boundary_category[
            negative_mask & (boundary_distance >= 4) & (boundary_distance <= 10)
        ] = "neg_adjacent_4_10"
        boundary_category[negative_mask & (boundary_distance > 10)] = "neg_far_gt10"
        boundary_category[
            negative_mask & (boundary_distance < 0)
        ] = "neg_no_positive_on_contig"

        group_df["boundary_window_distance"] = boundary_distance
        group_df["boundary_category"] = boundary_category
        group_df["is_false_positive"] = (
            (true_labels == 0) & (predicted_labels == 1)
        ).astype(int)
        group_df["is_false_negative"] = (
            (true_labels == 1) & (predicted_labels == 0)
        ).astype(int)
        group_df["is_true_positive"] = (
            (true_labels == 1) & (predicted_labels == 1)
        ).astype(int)
        annotated_parts.append(group_df)

    return pd.concat(annotated_parts, ignore_index=True)


def maybe_metric(group_df: pd.DataFrame, column_name: str) -> Optional[float]:
    if column_name not in group_df.columns or len(group_df) == 0:
        return None
    return float(group_df[column_name].mean())


def summarize_by_category(annotated_df: pd.DataFrame) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for category, category_df in annotated_df.groupby("boundary_category", sort=False):
        row: Dict[str, object] = {
            "boundary_category": category,
            "rows": int(len(category_df)),
            "positive_rate": float(category_df["label"].mean()),
            "pred_positive_rate": float(category_df["pred_label"].mean()),
            "prob_mean": float(category_df["prob"].mean()),
            "false_positive_rate": float(category_df["is_false_positive"].mean()),
            "false_negative_rate": float(category_df["is_false_negative"].mean()),
            "true_positive_rate": float(category_df["is_true_positive"].mean()),
        }
        virus_rate = maybe_metric(category_df, "protein_any_virus_hallmark_gene")
        if virus_rate is not None:
            row["virus_hallmark_rate"] = virus_rate
        rows.append(row)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_csv", required=True)
    parser.add_argument("--out_csv", required=True)
    parser.add_argument("--summary_json", required=True)
    args = parser.parse_args()

    input_df = pd.read_csv(args.input_csv)
    required_columns = {
        "genome_id",
        "contig_id",
        "start",
        "end",
        "label",
        "pred_label",
        "prob",
    }
    missing = sorted(required_columns - set(input_df.columns))
    if missing:
        raise ValueError(f"input missing required columns: {missing}")

    annotated = assign_boundary_category(input_df)

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    annotated.to_csv(out_csv, index=False)

    summary = {
        "rows": int(len(annotated)),
        "categories": summarize_by_category(annotated),
    }

    summary_json = Path(args.summary_json)
    summary_json.parent.mkdir(parents=True, exist_ok=True)
    with open(summary_json, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"[done] wrote boundary csv -> {out_csv}")
    print(f"[done] wrote boundary summary -> {summary_json}")


if __name__ == "__main__":
    main()
