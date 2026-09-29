#!/usr/bin/env python3
import argparse
import json
from pathlib import Path
from typing import Dict, List

import pandas as pd


def rate(numerator: int, denominator: int) -> float:
    if denominator <= 0:
        return float("nan")
    return float(numerator / denominator)


def support_stratum_rows(
    data_frame: pd.DataFrame, support_column: str, true_label: int
) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    subset = data_frame[data_frame["label"].astype(int) == true_label].copy()
    if subset.empty:
        return rows

    for support_value, support_df in subset.groupby(support_column, dropna=False):
        support_int = int(support_value)
        total = int(len(support_df))
        predicted_positive = int((support_df["pred_label"].astype(int) == 1).sum())
        predicted_negative = int((support_df["pred_label"].astype(int) == 0).sum())
        rows.append(
            {
                "label": int(true_label),
                "support_column": support_column,
                "support_value": support_int,
                "n_windows": total,
                "predicted_positive": predicted_positive,
                "predicted_negative": predicted_negative,
                "positive_rate": rate(predicted_positive, total),
                "negative_rate": rate(predicted_negative, total),
            }
        )
    return rows


def boundary_rows(
    data_frame: pd.DataFrame, error_type: str, support_column: str
) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    subset = data_frame[data_frame["error_type"] == error_type].copy()
    if subset.empty:
        return rows

    for (boundary_category, support_value), boundary_df in subset.groupby(
        ["boundary_category", support_column], dropna=False
    ):
        rows.append(
            {
                "error_type": error_type,
                "support_column": support_column,
                "boundary_category": str(boundary_category),
                "support_value": int(support_value),
                "n_windows": int(len(boundary_df)),
                "mean_prob": float(boundary_df["prob"].astype(float).mean()),
            }
        )
    return rows


def boundary_stratum_rows(
    data_frame: pd.DataFrame, support_column: str, true_label: int
) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    subset = data_frame[data_frame["label"].astype(int) == true_label].copy()
    if subset.empty:
        return rows

    for (boundary_category, support_value), boundary_df in subset.groupby(
        ["boundary_category", support_column], dropna=False
    ):
        total = int(len(boundary_df))
        predicted_positive = int((boundary_df["pred_label"].astype(int) == 1).sum())
        predicted_negative = int((boundary_df["pred_label"].astype(int) == 0).sum())
        rows.append(
            {
                "label": int(true_label),
                "support_column": support_column,
                "boundary_category": str(boundary_category),
                "support_value": int(support_value),
                "n_windows": total,
                "predicted_positive": predicted_positive,
                "predicted_negative": predicted_negative,
                "positive_rate": rate(predicted_positive, total),
                "negative_rate": rate(predicted_negative, total),
                "mean_prob": float(boundary_df["prob"].astype(float).mean()),
            }
        )
    return rows


def build_summary(data_frame: pd.DataFrame) -> Dict[str, object]:
    required_columns = {
        "label",
        "pred_label",
        "prob",
        "error_type",
        "boundary_category",
        "protein_any_virus_hallmark_gene",
        "protein_any_marker_gene",
    }
    missing = sorted(required_columns - set(data_frame.columns))
    if missing:
        raise ValueError(f"input missing required columns: {missing}")

    summary: Dict[str, object] = {
        "n_windows": int(len(data_frame)),
        "n_positive_windows": int((data_frame["label"].astype(int) == 1).sum()),
        "n_negative_windows": int((data_frame["label"].astype(int) == 0).sum()),
    }

    support_rows: List[Dict[str, object]] = []
    for support_column in [
        "protein_any_virus_hallmark_gene",
        "protein_any_marker_gene",
    ]:
        support_rows.extend(support_stratum_rows(data_frame, support_column, true_label=1))
        support_rows.extend(support_stratum_rows(data_frame, support_column, true_label=0))
    summary["support_strata"] = support_rows

    positive_subset = data_frame[data_frame["label"].astype(int) == 1].copy()
    negative_subset = data_frame[data_frame["label"].astype(int) == 0].copy()

    summary["positive_support_recall"] = {
        "virus_hallmark_supported": rate(
            int(
                (
                    (positive_subset["protein_any_virus_hallmark_gene"].astype(int) == 1)
                    & (positive_subset["pred_label"].astype(int) == 1)
                ).sum()
            ),
            int((positive_subset["protein_any_virus_hallmark_gene"].astype(int) == 1).sum()),
        ),
        "virus_hallmark_unsupported": rate(
            int(
                (
                    (positive_subset["protein_any_virus_hallmark_gene"].astype(int) == 0)
                    & (positive_subset["pred_label"].astype(int) == 1)
                ).sum()
            ),
            int((positive_subset["protein_any_virus_hallmark_gene"].astype(int) == 0).sum()),
        ),
    }

    summary["negative_support_false_positive_rate"] = {
        "virus_hallmark_supported": rate(
            int(
                (
                    (negative_subset["protein_any_virus_hallmark_gene"].astype(int) == 1)
                    & (negative_subset["pred_label"].astype(int) == 1)
                ).sum()
            ),
            int((negative_subset["protein_any_virus_hallmark_gene"].astype(int) == 1).sum()),
        ),
        "virus_hallmark_unsupported": rate(
            int(
                (
                    (negative_subset["protein_any_virus_hallmark_gene"].astype(int) == 0)
                    & (negative_subset["pred_label"].astype(int) == 1)
                ).sum()
            ),
            int((negative_subset["protein_any_virus_hallmark_gene"].astype(int) == 0).sum()),
        ),
    }

    summary["error_support_breakdown"] = {
        "false_negatives_with_virus_support": int(
            (
                (data_frame["error_type"] == "fn")
                & (data_frame["protein_any_virus_hallmark_gene"].astype(int) == 1)
            ).sum()
        ),
        "false_negatives_without_virus_support": int(
            (
                (data_frame["error_type"] == "fn")
                & (data_frame["protein_any_virus_hallmark_gene"].astype(int) == 0)
            ).sum()
        ),
        "false_positives_with_virus_support": int(
            (
                (data_frame["error_type"] == "fp")
                & (data_frame["protein_any_virus_hallmark_gene"].astype(int) == 1)
            ).sum()
        ),
        "false_positives_without_virus_support": int(
            (
                (data_frame["error_type"] == "fp")
                & (data_frame["protein_any_virus_hallmark_gene"].astype(int) == 0)
            ).sum()
        ),
    }

    boundary_summary_rows: List[Dict[str, object]] = []
    for error_type in ["fp", "fn"]:
        boundary_summary_rows.extend(
            boundary_rows(data_frame, error_type, "protein_any_virus_hallmark_gene")
        )
    summary["boundary_error_support_rows"] = boundary_summary_rows

    boundary_strata_rows: List[Dict[str, object]] = []
    boundary_strata_rows.extend(
        boundary_stratum_rows(data_frame, "protein_any_virus_hallmark_gene", true_label=1)
    )
    boundary_strata_rows.extend(
        boundary_stratum_rows(data_frame, "protein_any_virus_hallmark_gene", true_label=0)
    )
    summary["boundary_support_strata"] = boundary_strata_rows
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Summarize geNomad-label-source sensitivity and support dependence from joined prediction tables."
    )
    parser.add_argument("--input_csv", required=True)
    parser.add_argument("--summary_json", required=True)
    parser.add_argument("--support_csv", required=True)
    parser.add_argument("--boundary_csv", required=True)
    args = parser.parse_args()

    input_df = pd.read_csv(args.input_csv)
    summary = build_summary(input_df)

    out_summary = Path(args.summary_json)
    out_support = Path(args.support_csv)
    out_boundary = Path(args.boundary_csv)
    out_summary.parent.mkdir(parents=True, exist_ok=True)

    out_summary.write_text(json.dumps(summary, indent=2) + "\n")
    pd.DataFrame(summary["support_strata"]).to_csv(out_support, index=False)
    boundary_rows_df = pd.concat(
        [
            pd.DataFrame(summary["boundary_support_strata"]).assign(
                row_type="boundary_support_strata"
            ),
            pd.DataFrame(summary["boundary_error_support_rows"]).assign(
                row_type="boundary_error_support_rows"
            ),
        ],
        ignore_index=True,
    )
    boundary_rows_df.to_csv(out_boundary, index=False)

    print(f"[done] wrote {out_summary}")
    print(f"[done] wrote {out_support}")
    print(f"[done] wrote {out_boundary}")


if __name__ == "__main__":
    main()
