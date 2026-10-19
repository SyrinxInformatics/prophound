import argparse
import json
from pathlib import Path
from typing import Dict, List

import pandas as pd


def build_hysteresis_windows(
    prob_values: List[float], seed_threshold: float, extend_threshold: float
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


def close_small_gaps(values: List[int], max_gap_windows: int) -> List[int]:
    if max_gap_windows <= 0 or not values:
        return values[:]
    filled_values = values[:]
    n_values = len(filled_values)
    index = 0
    while index < n_values:
        if filled_values[index] != 1:
            index += 1
            continue
        scan_index = index + 1
        while scan_index < n_values and filled_values[scan_index] == 1:
            scan_index += 1
        gap_start = scan_index
        while scan_index < n_values and filled_values[scan_index] == 0:
            scan_index += 1
        gap_length = scan_index - gap_start
        if (
            gap_length > 0
            and gap_length <= max_gap_windows
            and scan_index < n_values
            and filled_values[scan_index] == 1
        ):
            for gap_index in range(gap_start, scan_index):
                filled_values[gap_index] = 1
        index = scan_index
    return filled_values


def remove_short_runs(values: List[int], min_run_windows: int) -> List[int]:
    if min_run_windows <= 1 or not values:
        return values[:]
    filtered_values = values[:]
    n_values = len(filtered_values)
    index = 0
    while index < n_values:
        if filtered_values[index] != 1:
            index += 1
            continue
        scan_index = index + 1
        while scan_index < n_values and filtered_values[scan_index] == 1:
            scan_index += 1
        run_length = scan_index - index
        if run_length < min_run_windows:
            for run_index in range(index, scan_index):
                filtered_values[run_index] = 0
        index = scan_index
    return filtered_values


def postprocess_positive_windows(
    values: List[int],
    max_gap_windows: int,
    min_run_windows: int,
    trim_terminal_windows: int = 0,
) -> List[int]:
    processed_values = close_small_gaps(values, max_gap_windows)
    processed_values = remove_short_runs(processed_values, min_run_windows)
    if trim_terminal_windows > 0:
        processed_values = trim_run_edges(processed_values, trim_terminal_windows)
    return processed_values


def trim_run_edges(values: List[int], trim_terminal_windows: int) -> List[int]:
    if trim_terminal_windows <= 0 or not values:
        return values[:]
    trimmed_values = values[:]
    n_values = len(trimmed_values)
    index = 0
    while index < n_values:
        if trimmed_values[index] != 1:
            index += 1
            continue
        scan_index = index + 1
        while scan_index < n_values and trimmed_values[scan_index] == 1:
            scan_index += 1
        run_start = index
        run_end = scan_index
        run_length = run_end - run_start
        edge_trim = min(trim_terminal_windows, run_length // 2)
        for trim_index in range(run_start, run_start + edge_trim):
            trimmed_values[trim_index] = 0
        for trim_index in range(run_end - edge_trim, run_end):
            trimmed_values[trim_index] = 0
        index = scan_index
    return trimmed_values


def windows_to_regions(windows_df: pd.DataFrame, positive_column: str) -> pd.DataFrame:
    regions: List[Dict[str, object]] = []
    for (genome_id, contig_id), group_df in windows_df.groupby(
        ["genome_id", "contig_id"], sort=False
    ):
        group_df = group_df.sort_values(["start", "end"]).copy()
        positive_values = group_df[positive_column].astype(int).to_list()
        window_starts = group_df["start"].astype(int).to_list()
        window_ends = group_df["end"].astype(int).to_list()

        open_region = False
        reg_start = None
        reg_end = None
        start_idx = None
        end_idx = None
        n_windows = 0

        for index, is_positive in enumerate(positive_values):
            if is_positive and not open_region:
                open_region = True
                reg_start = window_starts[index]
                reg_end = window_ends[index]
                start_idx = index
                end_idx = index
                n_windows = 1
            elif is_positive and open_region:
                reg_end = max(reg_end, window_ends[index])
                end_idx = index
                n_windows += 1
            elif (not is_positive) and open_region:
                regions.append(
                    {
                        "genome_id": genome_id,
                        "contig_id": contig_id,
                        "start": int(reg_start),
                        "end": int(reg_end),
                        "start_window_idx": int(start_idx),
                        "end_window_idx": int(end_idx),
                        "n_windows": int(n_windows),
                    }
                )
                open_region = False
                reg_start = reg_end = start_idx = end_idx = None
                n_windows = 0

        if open_region:
            regions.append(
                {
                    "genome_id": genome_id,
                    "contig_id": contig_id,
                    "start": int(reg_start),
                    "end": int(reg_end),
                    "start_window_idx": int(start_idx),
                    "end_window_idx": int(end_idx),
                    "n_windows": int(n_windows),
                }
            )

    region_df = pd.DataFrame(regions)
    if region_df.empty:
        return region_df
    region_df["region_len_bp"] = region_df["end"] - region_df["start"]
    return region_df


def filter_regions(
    region_df: pd.DataFrame, min_region_windows: int, min_region_bp: int
) -> pd.DataFrame:
    if region_df.empty:
        return region_df
    filtered_regions = region_df.copy()
    if min_region_windows > 1:
        filtered_regions = filtered_regions[
            filtered_regions["n_windows"] >= min_region_windows
        ].copy()
    if min_region_bp > 0:
        filtered_regions = filtered_regions[
            filtered_regions["region_len_bp"] >= min_region_bp
        ].copy()
    return filtered_regions.reset_index(drop=True)


def overlap_len(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def iou(a_start: int, a_end: int, b_start: int, b_end: int) -> float:
    ov = overlap_len(a_start, a_end, b_start, b_end)
    if ov <= 0:
        return 0.0
    union = max(a_end, b_end) - min(a_start, b_start)
    return ov / union if union > 0 else 0.0


def classify_true_region(true_row: pd.Series, pred_grp: pd.DataFrame) -> Dict[str, object]:
    if pred_grp.empty:
        return {
            "match_type": "missed_true_region",
            "n_overlapping_pred_regions": 0,
            "best_iou": 0.0,
            "best_pred_start_error_bp": None,
            "best_pred_end_error_bp": None,
            "left_truncation_bp": None,
            "right_truncation_bp": None,
            "left_overextension_bp": None,
            "right_overextension_bp": None,
        }

    n_overlap = len(pred_grp)
    best = None
    best_iou = -1.0
    for _, pred in pred_grp.iterrows():
        score = iou(true_row.start, true_row.end, pred.start, pred.end)
        if score > best_iou:
            best_iou = score
            best = pred

    start_err = int(best.start - true_row.start)
    end_err = int(best.end - true_row.end)

    left_trunc = max(0, int(best.start - true_row.start))
    right_trunc = max(0, int(true_row.end - best.end))
    left_over = max(0, int(true_row.start - best.start))
    right_over = max(0, int(best.end - true_row.end))

    if n_overlap > 1:
        match_type = "fragmented_true_region"
    elif left_trunc > 0 or right_trunc > 0:
        match_type = "truncated_true_region"
    elif left_over > 0 or right_over > 0:
        match_type = "overextended_true_region"
    else:
        match_type = "well_matched_true_region"

    return {
        "match_type": match_type,
        "n_overlapping_pred_regions": int(n_overlap),
        "best_iou": float(best_iou),
        "best_pred_start_error_bp": start_err,
        "best_pred_end_error_bp": end_err,
        "left_truncation_bp": left_trunc,
        "right_truncation_bp": right_trunc,
        "left_overextension_bp": left_over,
        "right_overextension_bp": right_over,
    }


def classify_pred_region(pred_row: pd.Series, true_grp: pd.DataFrame) -> Dict[str, object]:
    if true_grp.empty:
        return {
            "match_type": "isolated_spurious_pred_region",
            "n_overlapping_true_regions": 0,
            "best_iou": 0.0,
            "best_true_start_error_bp": None,
            "best_true_end_error_bp": None,
        }

    n_overlap = len(true_grp)
    best = None
    best_iou = -1.0
    for _, true_row in true_grp.iterrows():
        score = iou(pred_row.start, pred_row.end, true_row.start, true_row.end)
        if score > best_iou:
            best_iou = score
            best = true_row

    start_err = int(pred_row.start - best.start)
    end_err = int(pred_row.end - best.end)

    if n_overlap > 1:
        match_type = "merged_pred_region"
    else:
        left_trunc = max(0, int(pred_row.start - best.start))
        right_trunc = max(0, int(best.end - pred_row.end))
        left_over = max(0, int(best.start - pred_row.start))
        right_over = max(0, int(pred_row.end - best.end))
        if left_trunc > 0 or right_trunc > 0:
            match_type = "truncated_pred_region"
        elif left_over > 0 or right_over > 0:
            match_type = "overextended_pred_region"
        else:
            match_type = "well_matched_pred_region"

    return {
        "match_type": match_type,
        "n_overlapping_true_regions": int(n_overlap),
        "best_iou": float(best_iou),
        "best_true_start_error_bp": start_err,
        "best_true_end_error_bp": end_err,
    }


def subset_overlap(
    region_df: pd.DataFrame, genome_id: str, contig_id: str, start: int, end: int
) -> pd.DataFrame:
    matching_regions = region_df[
        (region_df["genome_id"] == genome_id) & (region_df["contig_id"] == contig_id)
    ].copy()
    if matching_regions.empty:
        return matching_regions
    return matching_regions[
        (matching_regions["end"] > start) & (matching_regions["start"] < end)
    ].copy()


def summarize_regions(region_df: pd.DataFrame, type_column: str) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for match_type, match_df in region_df.groupby(type_column, sort=False):
        row: Dict[str, object] = {
            type_column: str(match_type),
            "count": int(len(match_df)),
        }
        for col in [
            "best_iou",
            "best_pred_start_error_bp",
            "best_pred_end_error_bp",
            "best_true_start_error_bp",
            "best_true_end_error_bp",
            "left_truncation_bp",
            "right_truncation_bp",
            "left_overextension_bp",
            "right_overextension_bp",
        ]:
            if col in match_df.columns:
                values = match_df[col].dropna()
                if not values.empty:
                    row[f"{col}_mean"] = float(values.mean())
        rows.append(row)
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input_csv", required=True)
    ap.add_argument("--threshold", type=float, required=True)
    ap.add_argument("--out_true_csv", required=True)
    ap.add_argument("--out_pred_csv", required=True)
    ap.add_argument("--summary_json", required=True)
    ap.add_argument(
        "--max_gap_windows",
        type=int,
        default=0,
        help="Fill predicted negative gaps of this many windows or fewer when flanked by positives.",
    )
    ap.add_argument(
        "--min_run_windows",
        type=int,
        default=1,
        help="Drop predicted positive runs shorter than this many consecutive windows after gap closing.",
    )
    ap.add_argument(
        "--min_pred_region_windows",
        type=int,
        default=1,
        help="Drop predicted regions shorter than this many windows after region construction.",
    )
    ap.add_argument(
        "--min_pred_region_bp",
        type=int,
        default=0,
        help="Drop predicted regions shorter than this many base pairs after region construction.",
    )
    ap.add_argument(
        "--trim_terminal_windows",
        type=int,
        default=0,
        help="Trim this many positive windows from each end of a predicted run after gap closing.",
    )
    ap.add_argument(
        "--builder",
        choices=["threshold", "hysteresis"],
        default="threshold",
        help="Window-to-region builder policy. 'threshold' reproduces the current threshold-first flow.",
    )
    ap.add_argument(
        "--seed_threshold",
        type=float,
        default=None,
        help="Seed threshold for hysteresis mode. Defaults to --threshold.",
    )
    ap.add_argument(
        "--extend_threshold",
        type=float,
        default=None,
        help="Extension threshold for hysteresis mode. Required when --builder hysteresis.",
    )
    args = ap.parse_args()

    predictions_df = pd.read_csv(args.input_csv)
    required = {"genome_id", "contig_id", "start", "end", "label", "prob"}
    missing = sorted(required - set(predictions_df.columns))
    if missing:
        raise ValueError(f"input missing required columns: {missing}")

    predictions_df = predictions_df.copy()
    predictions_df["true_positive_window"] = predictions_df["label"].astype(int)

    per_contig_parts: List[pd.DataFrame] = []
    seed_threshold = args.seed_threshold if args.seed_threshold is not None else args.threshold
    if args.builder == "hysteresis" and args.extend_threshold is None:
        raise ValueError("--extend_threshold is required when --builder hysteresis")

    for _, group_df in predictions_df.groupby(["genome_id", "contig_id"], sort=False):
        group_df = group_df.sort_values(["start", "end"]).copy()
        prob_values = group_df["prob"].astype(float).to_list()
        group_df["pred_positive_window_raw"] = (
            group_df["prob"].astype(float) >= args.threshold
        ).astype(int)

        if args.builder == "hysteresis":
            group_df["pred_positive_window"] = build_hysteresis_windows(
                prob_values,
                seed_threshold=seed_threshold,
                extend_threshold=float(args.extend_threshold),
            )
        else:
            predicted_values = group_df["pred_positive_window_raw"].astype(int).to_list()
            group_df["pred_positive_window"] = postprocess_positive_windows(
                predicted_values,
                max_gap_windows=args.max_gap_windows,
                min_run_windows=args.min_run_windows,
                trim_terminal_windows=args.trim_terminal_windows,
            )
        per_contig_parts.append(group_df)
    predictions_df = pd.concat(per_contig_parts, ignore_index=True)

    true_regions = windows_to_regions(predictions_df, "true_positive_window")
    predicted_regions = windows_to_regions(predictions_df, "pred_positive_window")
    predicted_regions = filter_regions(
        predicted_regions,
        min_region_windows=args.min_pred_region_windows,
        min_region_bp=args.min_pred_region_bp,
    )

    true_rows: List[Dict[str, object]] = []
    for row in true_regions.itertuples(index=False):
        overlapping_predictions = subset_overlap(
            predicted_regions, row.genome_id, row.contig_id, row.start, row.end
        )
        true_region_record = row._asdict()
        true_region_record.update(
            classify_true_region(pd.Series(true_region_record), overlapping_predictions)
        )
        true_rows.append(true_region_record)

    pred_rows: List[Dict[str, object]] = []
    for row in predicted_regions.itertuples(index=False):
        overlapping_truth = subset_overlap(
            true_regions, row.genome_id, row.contig_id, row.start, row.end
        )
        predicted_region_record = row._asdict()
        predicted_region_record.update(
            classify_pred_region(pd.Series(predicted_region_record), overlapping_truth)
        )
        pred_rows.append(predicted_region_record)

    true_df = pd.DataFrame(true_rows)
    pred_df = pd.DataFrame(pred_rows)

    out_true = Path(args.out_true_csv)
    out_pred = Path(args.out_pred_csv)
    out_true.parent.mkdir(parents=True, exist_ok=True)
    out_pred.parent.mkdir(parents=True, exist_ok=True)
    true_df.to_csv(out_true, index=False)
    pred_df.to_csv(out_pred, index=False)

    summary = {
        "input_csv": args.input_csv,
        "threshold": args.threshold,
        "builder": args.builder,
        "seed_threshold": seed_threshold,
        "extend_threshold": args.extend_threshold,
        "max_gap_windows": args.max_gap_windows,
        "min_run_windows": args.min_run_windows,
        "trim_terminal_windows": args.trim_terminal_windows,
        "min_pred_region_windows": args.min_pred_region_windows,
        "min_pred_region_bp": args.min_pred_region_bp,
        "n_true_regions": int(len(true_df)),
        "n_pred_regions": int(len(pred_df)),
        "true_region_summary": summarize_regions(true_df, "match_type"),
        "pred_region_summary": summarize_regions(pred_df, "match_type"),
    }

    out_summary = Path(args.summary_json)
    out_summary.parent.mkdir(parents=True, exist_ok=True)
    with open(out_summary, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"[done] wrote true-region errors -> {out_true}")
    print(f"[done] wrote pred-region errors -> {out_pred}")
    print(f"[done] wrote region summary -> {out_summary}")


if __name__ == "__main__":
    main()
