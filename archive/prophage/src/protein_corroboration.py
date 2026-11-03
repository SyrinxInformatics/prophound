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
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import pandas as pd


Interval = Tuple[int, int]


@dataclass
class GeneRecord:
    contig: str
    start: int
    end: int
    marker: int
    virus_hallmark: int
    plasmid_hallmark: int


def normalize_contig_id(raw: str) -> List[str]:
    out: List[str] = []
    token = raw.split()[0]
    if token:
        out.append(token)
    if "|" in token:
        last = token.split("|")[-1]
        if last and last not in out:
            out.append(last)
    if "." in token:
        base = token.rsplit(".", 1)[0]
        if base and base not in out:
            out.append(base)
    return out


def read_split_ids(splits_dir: str, split: str) -> List[str]:
    path = Path(splits_dir) / f"{split}.txt"
    with open(path, "r") as f:
        return [line.strip() for line in f if line.strip()]


def sniff_dialect(path: Path):
    with open(path, "r", newline="") as f:
        sample = f.read(4096)
    try:
        return csv.Sniffer().sniff(sample, delimiters="\t,;|")
    except Exception:

        class _D:
            delimiter = "\t"

        return _D()


def pick_column(fieldnames: Iterable[str], options: List[str]) -> Optional[str]:
    lowered = {name.lower().strip(): name for name in fieldnames}
    for opt in options:
        if opt in lowered:
            return lowered[opt]
    return None


def contig_from_gene_id(gene_id: str) -> str:
    gene_id = gene_id.strip()
    if not gene_id:
        return ""
    head, sep, tail = gene_id.rpartition("_")
    if sep and tail.isdigit() and head:
        return head
    return gene_id


def find_genes_tsv(genomad_root: Path, genome_id: str) -> Optional[Path]:
    exact = (
        genomad_root
        / genome_id
        / f"{genome_id}_annotate"
        / f"{genome_id}_genes.tsv"
    )
    if exact.exists():
        return exact

    candidates = list(genomad_root.glob("**/*_annotate/*_genes.tsv"))
    target = genome_id.lower()
    for cand in candidates:
        parent = cand.parent.name.replace("_annotate", "").lower()
        stem = cand.stem.replace("_genes", "").lower()
        above = cand.parent.parent.name.lower()
        if target in {parent, stem, above}:
            return cand
    return None


def parse_flag(raw: object) -> int:
    val = str(raw).strip().lower()
    if val in {"1", "true", "t", "yes", "y"}:
        return 1
    return 0


def load_gene_intervals(tsv_path: Path) -> Dict[str, List[GeneRecord]]:
    dialect = sniff_dialect(tsv_path)
    genes_by_contig: Dict[str, List[GeneRecord]] = {}

    with open(tsv_path, "r", newline="") as f:
        reader = csv.DictReader(f, dialect=dialect)
        if not reader.fieldnames:
            raise ValueError(f"{tsv_path} is empty or missing headers")

        k_contig = pick_column(
            reader.fieldnames,
            [
                "seq_name",
                "sequence",
                "contig",
                "source_seq",
                "chrom",
                "scaffold",
                "sequence_name",
            ],
        )
        k_start = pick_column(
            reader.fieldnames,
            ["start", "gene_start", "begin", "pos_start"],
        )
        k_end = pick_column(
            reader.fieldnames,
            ["end", "gene_end", "stop", "pos_end"],
        )
        k_gene = pick_column(reader.fieldnames, ["gene", "gene_id", "id"])
        k_marker = pick_column(reader.fieldnames, ["marker"])
        k_virus_hallmark = pick_column(reader.fieldnames, ["virus_hallmark"])
        k_plasmid_hallmark = pick_column(reader.fieldnames, ["plasmid_hallmark"])
        if not (k_start and k_end):
            raise ValueError(
                f"{tsv_path} missing start/end columns; headers={reader.fieldnames}"
            )
        if not (k_contig or k_gene):
            raise ValueError(
                f"{tsv_path} missing contig or gene identifier column; headers={reader.fieldnames}"
            )

        for row in reader:
            try:
                if k_contig:
                    contig = str(row[k_contig]).strip()
                else:
                    contig = contig_from_gene_id(str(row[k_gene]))
                start = int(float(str(row[k_start]).strip()))
                end = int(float(str(row[k_end]).strip()))
            except Exception:
                continue
            if not contig:
                continue
            # geNomad annotate tables are typically 1-based inclusive.
            start = start - 1
            if end <= start:
                continue
            marker = 0
            if k_marker and str(row[k_marker]).strip().upper() != "NA":
                marker = 1
            virus_hallmark = parse_flag(row[k_virus_hallmark]) if k_virus_hallmark else 0
            plasmid_hallmark = (
                parse_flag(row[k_plasmid_hallmark]) if k_plasmid_hallmark else 0
            )
            genes_by_contig.setdefault(contig, []).append(
                GeneRecord(
                    contig=contig,
                    start=start,
                    end=end,
                    marker=marker,
                    virus_hallmark=virus_hallmark,
                    plasmid_hallmark=plasmid_hallmark,
                )
            )

    for contig in list(genes_by_contig.keys()):
        genes_by_contig[contig].sort(key=lambda gene: (gene.start, gene.end))
    return genes_by_contig


def merge_intervals(intervals: List[Interval]) -> List[Interval]:
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [intervals[0]]
    for start, end in intervals[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))
    return merged


def interval_distance(window: Interval, gene: GeneRecord) -> int:
    w_start, w_end = window
    g_start, g_end = gene.start, gene.end
    if g_end <= w_start:
        return w_start - g_end
    if g_start >= w_end:
        return g_start - w_end
    return 0


def window_gene_features(
    window: Interval, genes: List[GeneRecord], start_idx: int
) -> Tuple[int, Dict[str, float]]:
    w_start, w_end = window
    idx = start_idx
    n_genes = len(genes)

    while idx < n_genes and genes[idx].end <= w_start:
        idx += 1

    overlaps: List[Interval] = []
    overlap_count = 0
    marker_count = 0
    virus_hallmark_count = 0
    plasmid_hallmark_count = 0
    scan_idx = idx
    while scan_idx < n_genes and genes[scan_idx].start < w_end:
        gene = genes[scan_idx]
        g_start, g_end = gene.start, gene.end
        left = max(w_start, g_start)
        right = min(w_end, g_end)
        if right > left:
            overlap_count += 1
            overlaps.append((left, right))
            marker_count += gene.marker
            virus_hallmark_count += gene.virus_hallmark
            plasmid_hallmark_count += gene.plasmid_hallmark
        scan_idx += 1

    merged = merge_intervals(overlaps)
    overlap_bp = sum(end - start for start, end in merged)
    window_len = max(w_end - w_start, 1)

    nearest = None
    if overlap_count > 0:
        nearest = 0
    else:
        candidates: List[int] = []
        if idx > 0:
            candidates.append(interval_distance(window, genes[idx - 1]))
        if idx < n_genes:
            candidates.append(interval_distance(window, genes[idx]))
        nearest = min(candidates) if candidates else -1

    feats = {
        "protein_any_gene_overlap": int(overlap_count > 0),
        "protein_overlap_gene_count": overlap_count,
        "protein_overlap_bp": overlap_bp,
        "protein_overlap_frac": overlap_bp / float(window_len),
        "protein_nearest_gene_bp": nearest,
        "protein_marker_gene_count": marker_count,
        "protein_virus_hallmark_gene_count": virus_hallmark_count,
        "protein_plasmid_hallmark_gene_count": plasmid_hallmark_count,
        "protein_any_marker_gene": int(marker_count > 0),
        "protein_any_virus_hallmark_gene": int(virus_hallmark_count > 0),
        "protein_any_plasmid_hallmark_gene": int(plasmid_hallmark_count > 0),
    }
    return idx, feats


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest_csv", required=True)
    ap.add_argument("--genomad_root", required=True)
    ap.add_argument("--out_csv", required=True)
    ap.add_argument("--splits_dir", default=None)
    ap.add_argument("--split", choices=["train", "val", "test"], default=None)
    ap.add_argument("--limit_rows", type=int, default=0)
    args = ap.parse_args()

    usecols = ["genome_id", "contig_id", "start", "end", "label", "neg_type"]
    manifest = pd.read_csv(args.manifest_csv, usecols=usecols)

    if args.split:
        if not args.splits_dir:
            raise ValueError("--splits_dir is required when --split is set")
        split_ids = set(read_split_ids(args.splits_dir, args.split))
        manifest = manifest[manifest["genome_id"].isin(split_ids)].copy()

    if args.limit_rows > 0:
        manifest = manifest.head(args.limit_rows).copy()

    manifest["start"] = manifest["start"].astype(int)
    manifest["end"] = manifest["end"].astype(int)
    manifest = manifest.sort_values(["genome_id", "contig_id", "start", "end"]).reset_index(drop=True)

    genomad_root = Path(args.genomad_root)
    out_path = Path(args.out_csv)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    rows_out: List[Dict[str, object]] = []

    for genome_id, genome_df in manifest.groupby("genome_id", sort=False):
        genes_tsv = find_genes_tsv(genomad_root, str(genome_id))
        genes_by_contig: Dict[str, List[GeneRecord]] = {}
        genes_found = int(genes_tsv is not None)
        if genes_tsv is not None:
            genes_by_contig = load_gene_intervals(genes_tsv)

        normalized_gene_map: Dict[str, List[GeneRecord]] = {}
        for contig, intervals in genes_by_contig.items():
            for cand in normalize_contig_id(contig):
                normalized_gene_map.setdefault(cand, intervals)

        for contig_id, contig_df in genome_df.groupby("contig_id", sort=False):
            genes = genes_by_contig.get(contig_id)
            if genes is None:
                for cand in normalize_contig_id(str(contig_id)):
                    if cand in normalized_gene_map:
                        genes = normalized_gene_map[cand]
                        break
            if genes is None:
                genes = []

            start_idx = 0
            for row in contig_df.itertuples(index=False):
                start_idx, feats = window_gene_features((row.start, row.end), genes, start_idx)
                rows_out.append(
                    {
                        "genome_id": row.genome_id,
                        "contig_id": row.contig_id,
                        "start": row.start,
                        "end": row.end,
                        "label": row.label,
                        "neg_type": row.neg_type,
                        "protein_genes_tsv_found": genes_found,
                        **feats,
                    }
                )

    pd.DataFrame(rows_out).to_csv(out_path, index=False)
    print(f"[done] wrote {len(rows_out):,} rows -> {out_path}")


if __name__ == "__main__":
    main()
