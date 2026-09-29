import argparse
import json
from pathlib import Path
from typing import Dict, List

import pandas as pd


JOIN_KEYS = ["genome_id", "contig_id", "start", "end", "label", "neg_type"]


def derive_error_type(df: pd.DataFrame) -> pd.Series:
    pred = df["pred_label"].astype(int)
    label = df["label"].astype(int)
    out = pd.Series(index=df.index, dtype="object")
    out[(label == 1) & (pred == 1)] = "tp"
    out[(label == 1) & (pred == 0)] = "fn"
    out[(label == 0) & (pred == 1)] = "fp"
    out[(label == 0) & (pred == 0)] = "tn"
    return out


def maybe_float(value: object) -> float:
    return float(value) if pd.notna(value) else float("nan")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--predictions_csv", required=True)
    ap.add_argument("--protein_csv", required=True)
    ap.add_argument("--out_csv", required=True)
    ap.add_argument("--summary_json", required=True)
    args = ap.parse_args()

    preds = pd.read_csv(args.predictions_csv)
    prots = pd.read_csv(args.protein_csv)

    joined = preds.merge(
        prots,
        on=JOIN_KEYS,
        how="inner",
        validate="one_to_one",
    )
    joined["error_type"] = derive_error_type(joined)

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    joined.to_csv(out_csv, index=False)

    summary: Dict[str, object] = {
        "rows_predictions": int(len(preds)),
        "rows_protein": int(len(prots)),
        "rows_joined": int(len(joined)),
        "join_keys": JOIN_KEYS,
        "checkpoint_values": sorted(str(x) for x in joined["checkpoint"].dropna().unique().tolist()),
        "split_values": sorted(str(x) for x in joined["split"].dropna().unique().tolist()),
    }

    error_rows: List[Dict[str, object]] = []
    for error_type, grp in joined.groupby("error_type", sort=False):
        error_rows.append(
            {
                "error_type": error_type,
                "count": int(len(grp)),
                "prob_mean": maybe_float(grp["prob"].mean()),
                "any_virus_hallmark_rate": maybe_float(grp["protein_any_virus_hallmark_gene"].mean()),
                "virus_hallmark_gene_count_mean": maybe_float(
                    grp["protein_virus_hallmark_gene_count"].mean()
                ),
                "any_marker_rate": maybe_float(grp["protein_any_marker_gene"].mean()),
            }
        )
    summary["by_error_type"] = error_rows

    fp = joined[joined["error_type"] == "fp"]
    fn = joined[joined["error_type"] == "fn"]
    summary["headline_findings"] = {
        "fp_without_virus_hallmark_rate": maybe_float(
            1.0 - fp["protein_any_virus_hallmark_gene"].mean()
        )
        if len(fp)
        else float("nan"),
        "fn_with_virus_hallmark_rate": maybe_float(
            fn["protein_any_virus_hallmark_gene"].mean()
        )
        if len(fn)
        else float("nan"),
    }

    fp_by_neg_type: List[Dict[str, object]] = []
    for neg_type, grp in fp.groupby("neg_type", sort=False):
        fp_by_neg_type.append(
            {
                "neg_type": str(neg_type),
                "count": int(len(grp)),
                "no_virus_rate": maybe_float(1.0 - grp["protein_any_virus_hallmark_gene"].mean()),
                "prob_mean": maybe_float(grp["prob"].mean()),
            }
        )
    summary["fp_by_neg_type"] = fp_by_neg_type

    fn_by_neg_type: List[Dict[str, object]] = []
    for neg_type, grp in fn.groupby("neg_type", sort=False):
        fn_by_neg_type.append(
            {
                "neg_type": str(neg_type),
                "count": int(len(grp)),
                "virus_rate": maybe_float(grp["protein_any_virus_hallmark_gene"].mean()),
                "prob_mean": maybe_float(grp["prob"].mean()),
            }
        )
    summary["fn_by_neg_type"] = fn_by_neg_type

    summary_json = Path(args.summary_json)
    summary_json.parent.mkdir(parents=True, exist_ok=True)
    with open(summary_json, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"[done] wrote joined csv -> {out_csv}")
    print(f"[done] wrote summary json -> {summary_json}")


if __name__ == "__main__":
    main()
