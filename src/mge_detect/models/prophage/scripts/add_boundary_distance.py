from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def load_positive_edges(positives_dir: str, genome_id: str) -> dict:
    """Return {contig_id: sorted array of true-region edge positions} for one genome.

    Reads a plain BED file (chrom, start, end, ...) of the genome's true
    positive intervals only -- e.g. positives/by_genome/<genome_id>.bed on the
    HPCC. Each interval contributes both its start and end as an "edge".

    The chrom column in these files is "<genome_id>|<contig_id>" (verified
    2026-09-08 against positives/by_genome/GCF_000005845.2_ASM584v2_genomic.bed),
    not the bare contig_id the manifest uses -- strip the genome_id prefix so
    lookups actually match. Getting this wrong silently produces the "no
    positive interval" sentinel for every window instead of an error: an
    earlier run of this script did exactly that for all 1,567 genomes that do
    have positives (boundary_dist_bp came out constant at the sentinel value
    for all ~10.4M rows), and it still exited 0.
    """
    bed_path = Path(positives_dir) / f"{genome_id}.bed"
    edges_by_contig: dict = {}
    if not bed_path.exists():
        return edges_by_contig
    with open(bed_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            contig = chrom.split("|", 1)[1] if "|" in chrom else chrom
            edges_by_contig.setdefault(contig, []).extend([start, end])
    return {contig: np.array(sorted(v), dtype=np.int64) for contig, v in edges_by_contig.items()}


def nearest_edge_distance(centers: np.ndarray, edges: np.ndarray) -> np.ndarray:
    """Unsigned distance from each center position to the nearest edge.

    Unsigned deliberately: the model this feeds (RCCNN with predict_boundary)
    combines forward and reverse-complement embeddings symmetrically, so a
    signed "which direction" target would flip meaning under that view and
    get scrambled. See model_rc_cnn.py's RCCNN docstring/comments.
    """
    if edges.size == 0:
        return np.full(len(centers), np.inf)
    idx = np.searchsorted(edges, centers)
    left_idx = np.clip(idx - 1, 0, len(edges) - 1)
    right_idx = np.clip(idx, 0, len(edges) - 1)
    dist_left = np.abs(centers - edges[left_idx])
    dist_right = np.abs(centers - edges[right_idx])
    return np.minimum(dist_left, dist_right)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Add a boundary_dist_bp column (unsigned distance in bp from each "
            "window's center to the nearest true region edge) to a prophage "
            "training manifest, for the optional RCCNN boundary-regression "
            "head (model.predict_boundary in train_prophage.py)."
        )
    )
    parser.add_argument("--manifest_csv", required=True)
    parser.add_argument(
        "--positives_dir",
        required=True,
        help=(
            "Directory of per-genome BED files with true positive-only "
            "intervals, one file per genome named <genome_id>.bed "
            "(e.g. positives/by_genome on the HPCC)."
        ),
    )
    parser.add_argument("--out_csv", required=True)
    args = parser.parse_args()

    print(f"[boundary] reading manifest: {args.manifest_csv}", flush=True)
    df = pd.read_csv(args.manifest_csv)
    boundary_dist_bp = np.full(len(df), np.inf, dtype=np.float64)

    genome_groups = df.groupby("genome_id", sort=False)
    n_genomes = df["genome_id"].nunique()
    genomes_with_no_bed = 0

    for genome_index, (genome_id, genome_df) in enumerate(genome_groups, start=1):
        if genome_index % 200 == 0 or genome_index == n_genomes:
            print(f"[boundary] {genome_index}/{n_genomes} genomes processed", flush=True)

        edges_by_contig = load_positive_edges(args.positives_dir, genome_id)
        if not edges_by_contig:
            genomes_with_no_bed += 1

        for contig_id, contig_df in genome_df.groupby("contig_id", sort=False):
            edges = edges_by_contig.get(contig_id, np.array([], dtype=np.int64))
            centers = ((contig_df["start"] + contig_df["end"]) / 2.0).to_numpy()
            boundary_dist_bp[contig_df.index.to_numpy()] = nearest_edge_distance(centers, edges)

    print(f"[boundary] genomes with no positives BED found: {genomes_with_no_bed}/{n_genomes}")

    # No positive interval anywhere on that contig/genome -> genuinely
    # infinite distance. Saturate to a large sentinel rather than leaving inf
    # in the CSV; downstream normalization (min(dist, boundary_max_dist_bp))
    # treats any large value identically anyway.
    boundary_dist_bp[~np.isfinite(boundary_dist_bp)] = 1_000_000
    df["boundary_dist_bp"] = boundary_dist_bp.astype(np.int64)

    print(f"[boundary] writing: {args.out_csv}", flush=True)
    df.to_csv(args.out_csv, index=False)
    print("[boundary] boundary_dist_bp summary:")
    print(df["boundary_dist_bp"].describe())
    print("[boundary] done")


if __name__ == "__main__":
    main()
