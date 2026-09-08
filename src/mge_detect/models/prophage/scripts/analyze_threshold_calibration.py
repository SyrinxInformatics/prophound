import argparse
import json
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd


def compute_metrics(true_labels: np.ndarray, predicted_labels: np.ndarray) -> Dict[str, float]:
    tp = int(((predicted_labels == 1) & (true_labels == 1)).sum())
    fp = int(((predicted_labels == 1) & (true_labels == 0)).sum())
    fn = int(((predicted_labels == 0) & (true_labels == 1)).sum())
    tn = int(((predicted_labels == 0) & (true_labels == 0)).sum())

    precision = tp / (tp + fp) if (tp + fp) else 0.0
    recall = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = (
        2.0 * precision * recall / (precision + recall)
        if (precision + recall)
        else 0.0
    )
    pred_pos_rate = float((predicted_labels == 1).mean())
    return {
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "tn": tn,
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "predicted_positive_rate": pred_pos_rate,
    }


def rows_for_strategy(
    predictions_df: pd.DataFrame, thresholds: np.ndarray, high_conf_threshold: float
) -> pd.DataFrame:
    true_labels = predictions_df["label"].astype(int).to_numpy()
    probabilities = predictions_df["prob"].astype(float).to_numpy()
    support_labels = predictions_df["protein_any_virus_hallmark_gene"].astype(int).to_numpy()

    rows: List[Dict[str, float]] = []
    for threshold in thresholds:
        baseline_predictions = (probabilities >= threshold).astype(int)
        supported_predictions = (
            (probabilities >= threshold) & (support_labels == 1)
        ).astype(int)
        supported_or_highconf = (
            ((probabilities >= threshold) & (support_labels == 1))
            | (probabilities >= high_conf_threshold)
        ).astype(int)

        for strategy_name, predicted_labels in [
            ("baseline_prob", baseline_predictions),
            ("supported_prob", supported_predictions),
            ("supported_or_highconf", supported_or_highconf),
        ]:
            row = {"strategy": strategy_name, "threshold": float(threshold)}
            row.update(compute_metrics(true_labels, predicted_labels))
            rows.append(row)

    return pd.DataFrame(rows)


def pick_best_rows(curve_df: pd.DataFrame) -> Dict[str, object]:
    summary: Dict[str, object] = {}
    for strategy, strategy_df in curve_df.groupby("strategy", sort=False):
        strategy_df = strategy_df.sort_values(
            ["f1", "precision", "recall"], ascending=[False, False, False]
        )
        best_f1 = strategy_df.iloc[0].to_dict()

        best_precision_at_recall_08 = None
        recall_08_candidates = strategy_df[strategy_df["recall"] >= 0.80].sort_values(
            ["precision", "recall", "f1"], ascending=[False, False, False]
        )
        if not recall_08_candidates.empty:
            best_precision_at_recall_08 = recall_08_candidates.iloc[0].to_dict()

        best_precision_at_recall_06 = None
        recall_06_candidates = strategy_df[strategy_df["recall"] >= 0.60].sort_values(
            ["precision", "recall", "f1"], ascending=[False, False, False]
        )
        if not recall_06_candidates.empty:
            best_precision_at_recall_06 = recall_06_candidates.iloc[0].to_dict()

        summary[strategy] = {
            "best_f1": best_f1,
            "best_precision_at_recall_ge_0_80": best_precision_at_recall_08,
            "best_precision_at_recall_ge_0_60": best_precision_at_recall_06,
        }
    return summary


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input_csv", required=True)
    ap.add_argument("--out_csv", required=True)
    ap.add_argument("--summary_json", required=True)
    ap.add_argument("--threshold_min", type=float, default=0.05)
    ap.add_argument("--threshold_max", type=float, default=0.99)
    ap.add_argument("--threshold_step", type=float, default=0.01)
    ap.add_argument("--high_conf_threshold", type=float, default=0.99)
    args = ap.parse_args()

    predictions_df = pd.read_csv(args.input_csv)
    required = {"label", "prob", "protein_any_virus_hallmark_gene"}
    missing = sorted(required - set(predictions_df.columns))
    if missing:
        raise ValueError(f"input missing required columns: {missing}")

    thresholds = np.arange(
        args.threshold_min, args.threshold_max + (args.threshold_step / 2.0), args.threshold_step
    )
    curve = rows_for_strategy(predictions_df, thresholds, args.high_conf_threshold)

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    curve.to_csv(out_csv, index=False)

    summary = {
        "input_csv": args.input_csv,
        "high_conf_threshold": args.high_conf_threshold,
        "threshold_min": args.threshold_min,
        "threshold_max": args.threshold_max,
        "threshold_step": args.threshold_step,
        "best_operating_points": pick_best_rows(curve),
    }

    summary_json = Path(args.summary_json)
    summary_json.parent.mkdir(parents=True, exist_ok=True)
    with open(summary_json, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"[done] wrote calibration csv -> {out_csv}")
    print(f"[done] wrote calibration summary -> {summary_json}")


if __name__ == "__main__":
    main()
