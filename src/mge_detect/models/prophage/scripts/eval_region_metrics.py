"""Region-level recall/precision/boundary-error evaluation against geNomad.

Reusable, checked-in version of the ad hoc script used for the 18-genome
baseline evaluation in docs/model-quality-baseline.md (2026-09-08). Use this
to evaluate any new checkpoint (e.g. prophage_v1_positional_clean,
prophage_v1_positional_boundary) the same rigorous way, instead of trusting
a handful of hand-picked genomes.

Typical usage, once you have a checkpoint and its regions.tsv outputs for
each genome in eval_genome_list.txt (run `mge-detect predict` per genome,
one output dir per genome under some --predict_dir):

    python eval_region_metrics.py \
      --genomad_dir /path/to/genomad_out \
      --predict_dir /path/to/predict_outputs \
      --label "prophage_v1_positional_clean, single-tile"

--genomad_dir must contain, per genome accession `G`:
    G/G_find_proviruses/G_provirus.tsv   (geNomad's own output layout)
--predict_dir must contain, per genome accession `G`:
    G/regions.tsv                        (mge-detect predict's output)
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def load_genome_list(path: Path) -> list:
    with open(path) as f:
        return [line.strip() for line in f if line.strip()]


def load_true_regions(genomad_dir: Path, accession: str) -> pd.DataFrame:
    tsv = genomad_dir / accession / f"{accession}_find_proviruses" / f"{accession}_provirus.tsv"
    df = pd.read_csv(tsv, sep="\t")
    return df[["source_seq", "start", "end"]].rename(columns={"source_seq": "contig_id"})


def load_pred_regions(predict_dir: Path, accession: str) -> pd.DataFrame:
    tsv = predict_dir / accession / "regions.tsv"
    # predict_genome.py writes a bare newline (no header) when it finds zero
    # regions for a genome -- pd.read_csv chokes on that with EmptyDataError.
    if not tsv.exists() or tsv.stat().st_size <= 1:
        return pd.DataFrame(columns=["contig_id", "start", "end"])
    df = pd.read_csv(tsv, sep="\t")
    if df.empty:
        return pd.DataFrame(columns=["contig_id", "start", "end"])
    return df[["contig_id", "start", "end"]]


def overlap(a_start, a_end, b_start, b_end) -> float:
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def iou(a_start, a_end, b_start, b_end) -> float:
    inter = overlap(a_start, a_end, b_start, b_end)
    union = (a_end - a_start) + (b_end - b_start) - inter
    return inter / union if union > 0 else 0.0


def evaluate(genome_list: list, genomad_dir: Path, predict_dir: Path, label: str, strata: list | None = None) -> None:
    n_true_total = 0
    n_pred_total = 0
    n_true_matched = 0
    n_pred_matched = 0
    true_region_rows: list = []
    start_errors: list = []
    end_errors: list = []
    ious: list = []
    per_genome_rows = []

    for accession in genome_list:
        true_df = load_true_regions(genomad_dir, accession)
        pred_df = load_pred_regions(predict_dir, accession)
        n_true_total += len(true_df)
        n_pred_total += len(pred_df)

        true_matched_flags = [False] * len(true_df)
        pred_matched_flags = [False] * len(pred_df)

        pairs = []
        for ti, trow in true_df.reset_index(drop=True).iterrows():
            for pi, prow in pred_df.reset_index(drop=True).iterrows():
                if trow["contig_id"] != prow["contig_id"]:
                    continue
                ov = overlap(trow["start"], trow["end"], prow["start"], prow["end"])
                if ov > 0:
                    pairs.append((iou(trow["start"], trow["end"], prow["start"], prow["end"]), ti, pi))
        pairs.sort(reverse=True)
        used_true, used_pred = set(), set()
        for iou_val, ti, pi in pairs:
            if ti in used_true or pi in used_pred:
                continue
            used_true.add(ti)
            used_pred.add(pi)
            true_matched_flags[ti] = True
            pred_matched_flags[pi] = True
            trow = true_df.iloc[ti]
            prow = pred_df.iloc[pi]
            start_errors.append(prow["start"] - trow["start"])
            end_errors.append(prow["end"] - trow["end"])
            ious.append(iou_val)

        n_true_matched += sum(true_matched_flags)
        n_pred_matched += sum(pred_matched_flags)

        for ti, trow in true_df.reset_index(drop=True).iterrows():
            true_region_rows.append(
                {
                    "accession": accession,
                    "len_bp": trow["end"] - trow["start"],
                    "matched": bool(true_matched_flags[ti]),
                }
            )

        per_genome_rows.append(
            {
                "accession": accession,
                "n_true": len(true_df),
                "n_pred": len(pred_df),
                "n_matched": sum(true_matched_flags),
                "n_missed": len(true_df) - sum(true_matched_flags),
                "n_spurious": len(pred_df) - sum(pred_matched_flags),
            }
        )

    recall = n_true_matched / n_true_total if n_true_total else float("nan")
    precision = n_pred_matched / n_pred_total if n_pred_total else float("nan")
    start_errors_arr = np.array(start_errors)
    end_errors_arr = np.array(end_errors)
    abs_start = np.abs(start_errors_arr)
    abs_end = np.abs(end_errors_arr)

    print(f"\n===== {label} =====")
    print(f"total true regions: {n_true_total}, total predicted regions: {n_pred_total}")
    print(f"region-level recall: {recall:.3f}  ({n_true_matched}/{n_true_total})")
    print(f"region-level precision: {precision:.3f}  ({n_pred_matched}/{n_pred_total})")
    print(f"matched pairs: {len(ious)}")
    if len(ious):
        print(f"IoU: mean={np.mean(ious):.3f} median={np.median(ious):.3f} min={np.min(ious):.3f}")
        print(
            f"abs start error (bp): mean={abs_start.mean():.0f} median={np.median(abs_start):.0f} "
            f"p90={np.percentile(abs_start, 90):.0f} max={abs_start.max():.0f}"
        )
        print(
            f"abs end error (bp):   mean={abs_end.mean():.0f} median={np.median(abs_end):.0f} "
            f"p90={np.percentile(abs_end, 90):.0f} max={abs_end.max():.0f}"
        )
        print(
            f"start error sign: {(start_errors_arr < 0).sum()} undershoot(-) / "
            f"{(start_errors_arr > 0).sum()} overshoot(+) / {(start_errors_arr == 0).sum()} exact"
        )
        print(
            f"end error sign:   {(end_errors_arr < 0).sum()} undershoot(-) / "
            f"{(end_errors_arr > 0).sum()} overshoot(+) / {(end_errors_arr == 0).sum()} exact"
        )

    if strata and true_region_rows:
        tr = pd.DataFrame(true_region_rows)
        tr["kb"] = tr.len_bp / 1000.0
        print("\nrecall by element size class:")
        print(f"{'stratum':>10}{'n_true':>8}{'matched':>9}{'recall':>9}")
        for lo, hi in strata:
            sub = tr[(tr.kb >= lo) & (tr.kb < hi)]
            if sub.empty:
                continue
            lab = f"{lo:g}-{hi:g}kb" if hi != float("inf") else f"{lo:g}kb+"
            print(f"{lab:>10}{len(sub):>8}{int(sub.matched.sum()):>9}"
                  f"{sub.matched.sum() / len(sub) * 100:>8.1f}%")

    per_genome_df = pd.DataFrame(per_genome_rows)
    print("\nper-genome breakdown:")
    print(per_genome_df.to_string(index=False))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--genome_list",
        default=str(Path(__file__).parent / "eval_genome_list.txt"),
        help="One genome accession per line (default: the checked-in 18-genome baseline list).",
    )
    parser.add_argument("--genomad_dir", required=True)
    parser.add_argument("--predict_dir", required=True)
    parser.add_argument("--label", default="eval run")
    parser.add_argument(
        "--strata",
        default="",
        help=(
            "Comma-separated element-length breakpoints in kb (e.g. '5,8,10,15,30,60') to "
            "additionally report recall per size class. Use with the held-out size-stratified "
            "sample from build_eval_sample.py; the pooled recall above hides exactly the "
            "small-element behaviour that sample exists to measure."
        ),
    )
    args = parser.parse_args()

    strata = None
    if args.strata:
        edges = [float(x) for x in args.strata.split(",")]
        strata = list(zip(edges, edges[1:] + [float("inf")]))

    genome_list = load_genome_list(Path(args.genome_list))
    evaluate(genome_list, Path(args.genomad_dir), Path(args.predict_dir), args.label, strata)


if __name__ == "__main__":
    main()
