import argparse
import json
from pathlib import Path
from typing import Dict, List

import pandas as pd


JOIN_KEYS = ["genome_id", "contig_id", "start", "end", "label", "neg_type"]


def derive_error_type(joined_df: pd.DataFrame) -> pd.Series:
    predicted_labels = joined_df["pred_label"].astype(int)
    true_labels = joined_df["label"].astype(int)
    error_type = pd.Series(index=joined_df.index, dtype="object")
    error_type[(true_labels == 1) & (predicted_labels == 1)] = "tp"
    error_type[(true_labels == 1) & (predicted_labels == 0)] = "fn"
    error_type[(true_labels == 0) & (predicted_labels == 1)] = "fp"
    error_type[(true_labels == 0) & (predicted_labels == 0)] = "tn"
    return error_type


def maybe_float(value: object) -> float:
    return float(value) if pd.notna(value) else float("nan")


def binary_metrics(joined_df: pd.DataFrame, prediction_col: str) -> Dict[str, object]:
    predicted_labels = joined_df[prediction_col].astype(int)
    true_labels = joined_df["label"].astype(int)
    tp = int(((predicted_labels == 1) & (true_labels == 1)).sum())
    fp = int(((predicted_labels == 1) & (true_labels == 0)).sum())
    fn = int(((predicted_labels == 0) & (true_labels == 1)).sum())
    precision = tp / (tp + fp) if (tp + fp) else float("nan")
    recall = tp / (tp + fn) if (tp + fn) else float("nan")
    return {
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "predicted_positive_rows": int((predicted_labels == 1).sum()),
        "precision": precision,
        "recall": recall,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--predictions_csv", required=True)
    parser.add_argument("--protein_csv", required=True)
    parser.add_argument("--out_csv", required=True)
    parser.add_argument("--summary_json", required=True)
    args = parser.parse_args()

    predictions_df = pd.read_csv(args.predictions_csv)
    protein_df = pd.read_csv(args.protein_csv)

    joined_df = predictions_df.merge(
        protein_df,
        on=JOIN_KEYS,
        how="inner",
        validate="one_to_one",
    )
    joined_df["error_type"] = derive_error_type(joined_df)

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    joined_df.to_csv(out_csv, index=False)

    summary: Dict[str, object] = {
        "rows_predictions": int(len(predictions_df)),
        "rows_protein": int(len(protein_df)),
        "rows_joined": int(len(joined_df)),
        "join_keys": JOIN_KEYS,
        "checkpoint_values": sorted(
            str(x) for x in joined_df["checkpoint"].dropna().unique().tolist()
        ),
        "split_values": sorted(
            str(x) for x in joined_df["split"].dropna().unique().tolist()
        ),
    }

    error_rows: List[Dict[str, object]] = []
    for error_type, error_group in joined_df.groupby("error_type", sort=False):
        error_rows.append(
            {
                "error_type": error_type,
                "count": int(len(error_group)),
                "prob_mean": maybe_float(error_group["prob"].mean()),
                "any_virus_hallmark_rate": maybe_float(
                    error_group["protein_any_virus_hallmark_gene"].mean()
                ),
                "virus_hallmark_gene_count_mean": maybe_float(
                    error_group["protein_virus_hallmark_gene_count"].mean()
                ),
                "any_marker_rate": maybe_float(
                    error_group["protein_any_marker_gene"].mean()
                ),
            }
        )
    summary["by_error_type"] = error_rows

    false_positive_rows = joined_df[joined_df["error_type"] == "fp"]
    false_negative_rows = joined_df[joined_df["error_type"] == "fn"]
    summary["headline_findings"] = {
        "fp_without_virus_hallmark_rate": maybe_float(
            1.0 - false_positive_rows["protein_any_virus_hallmark_gene"].mean()
        )
        if len(false_positive_rows)
        else float("nan"),
        "fn_with_virus_hallmark_rate": maybe_float(
            false_negative_rows["protein_any_virus_hallmark_gene"].mean()
        )
        if len(false_negative_rows)
        else float("nan"),
    }

    fp_by_neg_type: List[Dict[str, object]] = []
    for neg_type, neg_type_group in false_positive_rows.groupby("neg_type", sort=False):
        fp_by_neg_type.append(
            {
                "neg_type": str(neg_type),
                "count": int(len(neg_type_group)),
                "no_virus_rate": maybe_float(
                    1.0 - neg_type_group["protein_any_virus_hallmark_gene"].mean()
                ),
                "prob_mean": maybe_float(neg_type_group["prob"].mean()),
            }
        )
    summary["fp_by_neg_type"] = fp_by_neg_type

    fn_by_neg_type: List[Dict[str, object]] = []
    for neg_type, neg_type_group in false_negative_rows.groupby("neg_type", sort=False):
        fn_by_neg_type.append(
            {
                "neg_type": str(neg_type),
                "count": int(len(neg_type_group)),
                "virus_rate": maybe_float(
                    neg_type_group["protein_any_virus_hallmark_gene"].mean()
                ),
                "prob_mean": maybe_float(neg_type_group["prob"].mean()),
            }
        )
    summary["fn_by_neg_type"] = fn_by_neg_type

    joined_df["rule_supported_only"] = (
        (joined_df["pred_label"].astype(int) == 1)
        & (joined_df["protein_any_virus_hallmark_gene"].astype(int) == 1)
    ).astype(int)
    joined_df["rule_highprob_or_supported"] = (
        (joined_df["pred_label"].astype(int) == 1)
        & (
            (joined_df["protein_any_virus_hallmark_gene"].astype(int) == 1)
            | (joined_df["prob"] >= 0.99)
        )
    ).astype(int)

    summary["optional_refinement_rules"] = {
        "baseline_pred_label": binary_metrics(joined_df, "pred_label"),
        "supported_only": binary_metrics(joined_df, "rule_supported_only"),
        "highprob_or_supported": binary_metrics(
            joined_df, "rule_highprob_or_supported"
        ),
    }

    summary_json = Path(args.summary_json)
    summary_json.parent.mkdir(parents=True, exist_ok=True)
    with open(summary_json, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"[done] wrote joined csv -> {out_csv}")
    print(f"[done] wrote summary json -> {summary_json}")


if __name__ == "__main__":
    main()
