#!/usr/bin/env python3
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


def best_overlap_stats(region_row: pd.Series, overlap_df: pd.DataFrame) -> Dict[str, object]:
    if overlap_df.empty:
        return {
            "n_overlaps": 0,
            "best_iou": 0.0,
            "detected": 0,
        }

    best_iou = 0.0
    for _, overlap_row in overlap_df.iterrows():
        best_iou = max(
            best_iou,
            iou(region_row.start, region_row.end, overlap_row.start, overlap_row.end),
        )
    return {
        "n_overlaps": int(len(overlap_df)),
        "best_iou": float(best_iou),
        "detected": 1,
    }


def true_regions_from_windows(predictions_df: pd.DataFrame) -> pd.DataFrame:
    region_rows: List[Dict[str, object]] = []
    for (genome_id, contig_id), contig_df in predictions_df.groupby(
        ["genome_id", "contig_id"],
        sort=False,
    ):
        contig_df = contig_df.sort_values(["start", "end"]).reset_index(drop=True)
        true_labels = contig_df["label"].astype(int).to_list()

        region_open = False
        start_index = None
        end_index = None
        for window_index, is_positive in enumerate(true_labels):
            if is_positive and not region_open:
                region_open = True
                start_index = window_index
                end_index = window_index
            elif is_positive and region_open:
                end_index = window_index
            elif (not is_positive) and region_open:
                region_slice = contig_df.iloc[start_index : end_index + 1].copy()
                region_rows.append(
                    {
                        "genome_id": genome_id,
                        "contig_id": contig_id,
                        "start": int(region_slice["start"].min()),
                        "end": int(region_slice["end"].max()),
                        "n_windows": int(len(region_slice)),
                        "any_virus_support": int(
                            region_slice["protein_any_virus_hallmark_gene"].astype(int).max()
                        ),
                        "mean_prob": float(region_slice["prob"].astype(float).mean()),
                        "boundary_categories": ",".join(
                            sorted(set(region_slice["boundary_category"].astype(str).tolist()))
                        ),
                    }
                )
                region_open = False
                start_index = None
                end_index = None

        if region_open:
            region_slice = contig_df.iloc[start_index : end_index + 1].copy()
            region_rows.append(
                {
                    "genome_id": genome_id,
                    "contig_id": contig_id,
                    "start": int(region_slice["start"].min()),
                    "end": int(region_slice["end"].max()),
                    "n_windows": int(len(region_slice)),
                    "any_virus_support": int(
                        region_slice["protein_any_virus_hallmark_gene"].astype(int).max()
                    ),
                    "mean_prob": float(region_slice["prob"].astype(float).mean()),
                    "boundary_categories": ",".join(
                        sorted(set(region_slice["boundary_category"].astype(str).tolist()))
                    ),
                }
            )

    return pd.DataFrame(region_rows)


def predicted_regions_from_windows(
    predictions_df: pd.DataFrame,
    seed_threshold: float,
    extend_threshold: float,
) -> pd.DataFrame:
    region_rows: List[Dict[str, object]] = []
    for (genome_id, contig_id), contig_df in predictions_df.groupby(
        ["genome_id", "contig_id"],
        sort=False,
    ):
        contig_df = contig_df.sort_values(["start", "end"]).reset_index(drop=True).copy()
        contig_df["pred_positive_window"] = build_hysteresis_windows(
            contig_df["prob"].astype(float).to_list(),
            seed_threshold=seed_threshold,
            extend_threshold=extend_threshold,
        )

        region_open = False
        start_index = None
        end_index = None
        for window_index, is_positive in enumerate(
            contig_df["pred_positive_window"].astype(int).to_list()
        ):
            if is_positive and not region_open:
                region_open = True
                start_index = window_index
                end_index = window_index
            elif is_positive and region_open:
                end_index = window_index
            elif (not is_positive) and region_open:
                region_slice = contig_df.iloc[start_index : end_index + 1].copy()
                region_rows.append(
                    {
                        "genome_id": genome_id,
                        "contig_id": contig_id,
                        "start": int(region_slice["start"].min()),
                        "end": int(region_slice["end"].max()),
                        "n_windows": int(len(region_slice)),
                        "any_virus_support": int(
                            region_slice["protein_any_virus_hallmark_gene"].astype(int).max()
                        ),
                        "mean_prob": float(region_slice["prob"].astype(float).mean()),
                        "max_prob": float(region_slice["prob"].astype(float).max()),
                        "boundary_categories": ",".join(
                            sorted(set(region_slice["boundary_category"].astype(str).tolist()))
                        ),
                    }
                )
                region_open = False
                start_index = None
                end_index = None

        if region_open:
            region_slice = contig_df.iloc[start_index : end_index + 1].copy()
            region_rows.append(
                {
                    "genome_id": genome_id,
                    "contig_id": contig_id,
                    "start": int(region_slice["start"].min()),
                    "end": int(region_slice["end"].max()),
                    "n_windows": int(len(region_slice)),
                    "any_virus_support": int(
                        region_slice["protein_any_virus_hallmark_gene"].astype(int).max()
                    ),
                    "mean_prob": float(region_slice["prob"].astype(float).mean()),
                    "max_prob": float(region_slice["prob"].astype(float).max()),
                    "boundary_categories": ",".join(
                        sorted(set(region_slice["boundary_category"].astype(str).tolist()))
                    ),
                }
            )

    return pd.DataFrame(region_rows)


def summarize_region_strata(region_df: pd.DataFrame, region_type: str) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for support_value, support_df in region_df.groupby("any_virus_support", dropna=False):
        rows.append(
            {
                "region_type": region_type,
                "support_value": int(support_value),
                "n_regions": int(len(support_df)),
                "mean_n_windows": float(support_df["n_windows"].mean()),
                "mean_prob": float(support_df["mean_prob"].mean()),
                "mean_best_iou": float(support_df["best_iou"].mean()) if "best_iou" in support_df else float("nan"),
                "detection_rate": float(support_df["detected"].mean()) if "detected" in support_df else float("nan"),
                "isolated_spurious_rate": float((support_df["match_type"] == "isolated_spurious_pred_region").mean()) if "match_type" in support_df else float("nan"),
            }
        )
    return rows


def build_summary(true_regions_df: pd.DataFrame, predicted_regions_df: pd.DataFrame) -> Dict[str, object]:
    true_region_rows: List[Dict[str, object]] = []
    for true_region in true_regions_df.itertuples(index=False):
        overlapping_predicted_regions = overlapping_regions(
            predicted_regions_df,
            true_region.genome_id,
            true_region.contig_id,
            true_region.start,
            true_region.end,
        )
        true_region_dict = true_region._asdict()
        true_region_dict.update(
            best_overlap_stats(pd.Series(true_region_dict), overlapping_predicted_regions)
        )
        true_region_rows.append(true_region_dict)
    true_region_eval_df = pd.DataFrame(true_region_rows)

    predicted_region_rows: List[Dict[str, object]] = []
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
            best_overlap_stats(pd.Series(predicted_region_dict), overlapping_true_regions)
        )
        predicted_region_dict["match_type"] = (
            "isolated_spurious_pred_region" if overlapping_true_regions.empty else "true_overlapping_pred_region"
        )
        predicted_region_rows.append(predicted_region_dict)
    predicted_region_eval_df = pd.DataFrame(predicted_region_rows)

    return {
        "n_true_regions": int(len(true_region_eval_df)),
        "n_pred_regions": int(len(predicted_region_eval_df)),
        "true_region_support_strata": summarize_region_strata(true_region_eval_df, "true_region"),
        "pred_region_support_strata": summarize_region_strata(predicted_region_eval_df, "pred_region"),
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run a region-level geNomad dependence check under the frozen strict builder."
    )
    parser.add_argument("--input_csv", required=True)
    parser.add_argument("--summary_json", required=True)
    parser.add_argument("--true_regions_csv", required=True)
    parser.add_argument("--pred_regions_csv", required=True)
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
        "boundary_category",
        "protein_any_virus_hallmark_gene",
    }
    missing_columns = sorted(required_columns - set(predictions_df.columns))
    if missing_columns:
        raise ValueError(f"input missing required columns: {missing_columns}")

    true_regions_df = true_regions_from_windows(predictions_df)
    predicted_regions_df = predicted_regions_from_windows(
        predictions_df,
        seed_threshold=args.seed_threshold,
        extend_threshold=args.extend_threshold,
    )
    summary = build_summary(true_regions_df, predicted_regions_df)
    summary["seed_threshold"] = args.seed_threshold
    summary["extend_threshold"] = args.extend_threshold

    out_summary = Path(args.summary_json)
    out_true = Path(args.true_regions_csv)
    out_pred = Path(args.pred_regions_csv)
    out_summary.parent.mkdir(parents=True, exist_ok=True)

    true_regions_df.to_csv(out_true, index=False)
    predicted_regions_df.to_csv(out_pred, index=False)
    out_summary.write_text(json.dumps(summary, indent=2) + "\n")

    print(f"[done] wrote {out_true}")
    print(f"[done] wrote {out_pred}")
    print(f"[done] wrote {out_summary}")


if __name__ == "__main__":
    main()
