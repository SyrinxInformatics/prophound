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


MODE_THRESHOLDS = {
    "exploratory": 0.95,
    "default": 0.995,
    "conservative": 0.999,
}


def compute_counts(predictions_df: pd.DataFrame, prediction_col: str) -> Dict[str, float]:
    true_labels = predictions_df["label"].astype(int)
    predicted_labels = predictions_df[prediction_col].astype(int)
    tp = int(((true_labels == 1) & (predicted_labels == 1)).sum())
    fp = int(((true_labels == 0) & (predicted_labels == 1)).sum())
    fn = int(((true_labels == 1) & (predicted_labels == 0)).sum())
    tn = int(((true_labels == 0) & (predicted_labels == 0)).sum())
    precision = tp / (tp + fp) if (tp + fp) else 0.0
    recall = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0
    return {
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "tn": tn,
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "predicted_positive_rows": int((predicted_labels == 1).sum()),
    }


def by_boundary(predictions_df: pd.DataFrame, prediction_col: str) -> List[Dict[str, float]]:
    rows: List[Dict[str, float]] = []
    predicted_labels = predictions_df[prediction_col].astype(int)
    for category, category_df in predictions_df.groupby("boundary_category", sort=False):
        true_labels = category_df["label"].astype(int)
        category_predictions = predicted_labels.loc[category_df.index]
        row: Dict[str, float] = {
            "boundary_category": str(category),
            "rows": int(len(category_df)),
            "pred_positive_rate": float(category_predictions.mean()),
        }
        if int(true_labels.iloc[0]) == 1:
            row["recall_within_category"] = float(category_predictions.mean())
            row["false_negative_rate_within_category"] = float(
                1.0 - category_predictions.mean()
            )
        else:
            row["false_positive_rate_within_category"] = float(
                category_predictions.mean()
            )
        rows.append(row)
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input_csv", required=True)
    ap.add_argument("--out_json", required=True)
    args = ap.parse_args()

    predictions_df = pd.read_csv(args.input_csv)
    required = {"label", "prob", "boundary_category"}
    missing = sorted(required - set(predictions_df.columns))
    if missing:
        raise ValueError(f"input missing required columns: {missing}")

    summary: Dict[str, object] = {
        "input_csv": args.input_csv,
        "mode_thresholds": MODE_THRESHOLDS,
        "modes": {},
    }

    for mode, threshold in MODE_THRESHOLDS.items():
        pred_col = f"pred_{mode}"
        predictions_df[pred_col] = (
            predictions_df["prob"].astype(float) >= threshold
        ).astype(int)
        summary["modes"][mode] = {
            "threshold": threshold,
            "overall": compute_counts(predictions_df, pred_col),
            "boundary": by_boundary(predictions_df, pred_col),
        }

    out = Path(args.out_json)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"[done] wrote operating mode summary -> {out}")


if __name__ == "__main__":
    main()
