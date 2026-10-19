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
from __future__ import annotations

import argparse
import csv
import math
import os
from dataclasses import dataclass
from typing import Dict, Iterable, List, Tuple


@dataclass
class Segment:
    start: int
    end: int
    seg_type: str


def iter_fasta(path: str) -> Iterable[Tuple[str, int]]:
    """Yield `(contig_id, contig_length)` pairs from a FASTA file."""
    contig_id = None
    contig_length = 0
    with open(path, "r") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if contig_id is not None:
                    yield contig_id, contig_length
                contig_id = line[1:].split()[0]
                contig_length = 0
            else:
                contig_length += len(line.strip())
        if contig_id is not None:
            yield contig_id, contig_length


def read_bed(path: str) -> Dict[str, List[Segment]]:
    """Read a BED-like label map keyed by contig or genome ID."""
    segments_by_key: Dict[str, List[Segment]] = {}
    with open(path, "r") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 4:
                raise ValueError(f"BED line has <4 columns in {path}: {line}")
            key, start, end, segment_type = (
                parts[0],
                int(parts[1]),
                int(parts[2]),
                parts[3],
            )
            segments_by_key.setdefault(key, []).append(
                Segment(start, end, segment_type)
            )

    for key, segments in segments_by_key.items():
        segments.sort(key=lambda segment: segment.start)
    return segments_by_key


def seg_overlap_bp(
    segment: Segment,
    window_start: int,
    window_end: int,
    coord_mode: str,
) -> int:
    """Return the overlap in base pairs between one segment and one window."""
    if coord_mode == "inclusive":
        segment_start, segment_end = segment.start, segment.end + 1
    elif coord_mode == "halfopen":
        segment_start, segment_end = segment.start, segment.end
    else:
        raise ValueError("coord_mode must be 'inclusive' or 'halfopen'")

    left = max(window_start, segment_start)
    right = min(window_end, segment_end)
    return max(0, right - left)


def stem(path: str) -> str:
    base = os.path.basename(path)
    for ext in (".fasta", ".fa", ".fna", ".fas", ".bed", ".tsv", ".txt"):
        if base.endswith(ext):
            return base[: -len(ext)]
    return base


def list_files(directory: str, extensions: Tuple[str, ...]) -> List[str]:
    return sorted(
        os.path.join(directory, filename)
        for filename in os.listdir(directory)
        if filename.endswith(extensions)
    )


def normalize_contig_id(raw: str) -> List[str]:
    """Generate candidate contig IDs for BED and FASTA matching."""
    normalized_ids: List[str] = []
    token = raw.split()[0]
    if token:
        normalized_ids.append(token)

    if "|" in token:
        last_token = token.split("|")[-1]
        if last_token and last_token not in normalized_ids:
            normalized_ids.append(last_token)

    if "." in token:
        versionless_id = token.rsplit(".", 1)[0]
        if versionless_id and versionless_id not in normalized_ids:
            normalized_ids.append(versionless_id)

    return normalized_ids


def pick_segment_key(
    contig_id: str,
    genome_id: str,
    segments_by_key: Dict[str, List[Segment]],
) -> str | None:
    """Choose the best segment-map key for a FASTA contig."""
    if contig_id in segments_by_key:
        return contig_id

    for normalized_id in normalize_contig_id(contig_id):
        if normalized_id in segments_by_key:
            return normalized_id

    if genome_id in segments_by_key:
        return genome_id
    return None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--genomes_dir", required=True, help="Directory of FASTA files")
    parser.add_argument("--beds_dir", required=True, help="Directory of BED maps")
    parser.add_argument("--out_csv", required=True, help="Output manifest CSV path")
    parser.add_argument("--window_len", type=int, default=10_000)
    parser.add_argument("--stride", type=int, default=1_000)
    parser.add_argument(
        "--coord_mode",
        choices=["inclusive", "halfopen"],
        default="inclusive",
    )
    parser.add_argument(
        "--overlap_min_frac",
        type=float,
        default=0.50,
        help="Label by overlap: minimum fraction of window to assign a class.",
    )
    parser.add_argument(
        "--overlap_min_bp",
        type=int,
        default=0,
        help="Label by overlap: minimum base pairs to assign a class.",
    )
    parser.add_argument(
        "--include_near_neg",
        action="store_true",
        help="Include near-negative rows in the manifest.",
    )
    parser.add_argument(
        "--min_contig_len",
        type=int,
        default=10_000,
        help="Skip contigs shorter than this threshold.",
    )
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.out_csv) or ".", exist_ok=True)

    fasta_files = list_files(args.genomes_dir, (".fasta", ".fa", ".fna", ".fas"))
    if not fasta_files:
        raise SystemExit(f"No FASTA files found in {args.genomes_dir}")

    bed_files = list_files(args.beds_dir, (".bed", ".tsv", ".txt"))
    if not bed_files:
        raise SystemExit(f"No BED-like files found in {args.beds_dir}")

    bed_by_genome_id = {stem(path): path for path in bed_files}
    min_contig_length = max(args.min_contig_len, args.window_len)

    n_windows = 0
    n_genomes_no_bed = 0
    n_contigs_no_map = 0
    n_contigs_short = 0

    with open(args.out_csv, "w", newline="") as out_handle:
        writer = csv.DictWriter(
            out_handle,
            fieldnames=[
                "genome_id",
                "contig_id",
                "start",
                "end",
                "label",
                "neg_type",
                "fasta_path",
            ],
        )
        writer.writeheader()

        for fasta_path in fasta_files:
            genome_id = stem(fasta_path)
            bed_path = bed_by_genome_id.get(genome_id)
            if bed_path is None:
                n_genomes_no_bed += 1
                continue

            segments_by_key = read_bed(bed_path)

            for contig_id, contig_length in iter_fasta(fasta_path):
                if contig_length < min_contig_length:
                    n_contigs_short += 1
                    continue

                segment_key = pick_segment_key(contig_id, genome_id, segments_by_key)
                if segment_key is None:
                    n_contigs_no_map += 1
                    continue

                segments = segments_by_key[segment_key]
                min_label_bp = max(
                    args.overlap_min_bp,
                    int(math.ceil(args.window_len * args.overlap_min_frac)),
                )

                for start in range(
                    0,
                    contig_length - args.window_len + 1,
                    args.stride,
                ):
                    end = start + args.window_len
                    overlap_by_type = {"positive": 0, "near_neg": 0, "far_neg": 0}

                    for segment in segments:
                        if segment.seg_type not in overlap_by_type:
                            continue
                        overlap_by_type[segment.seg_type] += seg_overlap_bp(
                            segment,
                            start,
                            end,
                            args.coord_mode,
                        )

                    segment_type = None
                    if overlap_by_type["positive"] >= min_label_bp:
                        segment_type = "positive"
                    elif overlap_by_type["near_neg"] >= min_label_bp:
                        # Prefer boundary-adjacent negatives over background negatives
                        # when both overlap thresholds are satisfied.
                        segment_type = "near_neg"
                    elif overlap_by_type["far_neg"] >= min_label_bp:
                        segment_type = "far_neg"

                    if segment_type is None:
                        continue
                    if segment_type == "near_neg" and not args.include_near_neg:
                        continue

                    writer.writerow(
                        {
                            "genome_id": genome_id,
                            "contig_id": contig_id,
                            "start": start,
                            "end": end,
                            "label": 1 if segment_type == "positive" else 0,
                            "neg_type": segment_type,
                            "fasta_path": fasta_path,
                        }
                    )
                    n_windows += 1

    print(f"Wrote {n_windows:,} windows -> {args.out_csv}")
    if n_genomes_no_bed:
        print(
            f"WARNING: {n_genomes_no_bed} FASTA file(s) had no matching BED by filename stem."
        )
    if n_contigs_no_map:
        print(
            f"WARNING: {n_contigs_no_map} contig(s) had no matching BED key."
        )
    if n_contigs_short:
        print(
            f"INFO: Skipped {n_contigs_short} contig(s) shorter than {min_contig_length} bp."
        )


if __name__ == "__main__":
    main()
