#!/usr/bin/env python3
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
"""
Build per-accession BED files from FASTA inputs and geNomad provirus calls.

Each output BED labels intervals as `positive`, `near_neg`, or `far_neg` and
the run also writes a summary CSV plus a `beds_pos/` subset containing only
accessions with at least one positive region.
"""
import argparse
import csv
import io
import re
import shutil
import subprocess
import sys
from pathlib import Path

# Types
Interval = tuple[int, int]  # half-open [start, end)

# FASTA utilities


def _iter_fasta_lines_text(path: Path):
    """Yield FASTA lines from plain text or gzip-compressed files."""
    import gzip as _gzip
    import shutil as _shutil

    def _is_gzip(p: Path) -> bool:
        try:
            with open(p, "rb") as fh:
                sig = fh.read(2)
            return sig == b"\x1f\x8b" or str(p).endswith(".gz")
        except Exception:
            return str(p).endswith(".gz")

    if not _is_gzip(path):
        with open(path, encoding="utf-8", errors="ignore") as f:
            for line in f:
                yield line.rstrip("\n")
        return

    # Try native gzip first.
    try:
        with _gzip.open(str(path), "rt", encoding="utf-8", errors="ignore") as f:
            for line in f:
                yield line.rstrip("\n")
        return
    except Exception:
        # Fall back to external gzip if needed.
        if _shutil.which("gzip") is not None:
            proc = subprocess.Popen(
                ["gzip", "-dc", str(path)],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            buf = io.BufferedReader(proc.stdout)
            for raw in buf:
                try:
                    yield raw.decode("utf-8", "ignore").rstrip("\n")
                except Exception:
                    continue
            proc.stdout.close()
            proc.wait()
            return
        # Last resort: salvage a partial read from a truncated gzip stream.
        try:
            with open(path, "rb") as fh:
                gf = _gzip.GzipFile(fileobj=fh, mode="rb")
                chunk = gf.read(65536)
                carry = ""
                while chunk:
                    s = chunk.decode("utf-8", "ignore")
                    lines = (carry + s).split("\n")
                    carry = lines.pop()
                    for ln in lines:
                        yield ln
                    try:
                        chunk = gf.read(65536)
                    except EOFError:
                        break
                if carry:
                    yield carry
        except Exception:
            return


def read_fasta_lengths(fasta_path: Path) -> dict[str, int]:
    """Return {contig_name: length} for a (possibly multi-contig) FASTA.
    Supports plain text and .gz (even if truncated)."""
    lengths: dict[str, int] = {}
    name: str | None = None
    length = 0

    for line in _iter_fasta_lines_text(fasta_path):
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                lengths[name] = length
            name = line[1:].strip().split()[0]
            length = 0
        else:
            length += len(line.strip())
    if name is not None:
        lengths[name] = length
    return lengths


def infer_accession_from_filename(p: Path) -> str:
    """Define accession as the FASTA filename without common FASTA/gzip extensions."""
    name = p.name
    for suf in (".fa.gz", ".fna.gz", ".fasta.gz", ".fa", ".fna", ".fasta"):
        if name.endswith(suf):
            return name[: -len(suf)]
    # Fall back to Path.stem if no common FASTA suffix matches.
    return p.stem


# Interval helpers
def merge_intervals(intervals: list[Interval], min_gap: int = 0) -> list[Interval]:
    """Merge overlapping or nearly adjacent intervals (gap <= min_gap)."""
    if not intervals:
        return []
    ints = sorted(intervals)
    merged: list[Interval] = [ints[0]]
    for s, e in ints[1:]:
        last_s, last_e = merged[-1]
        if s <= last_e + min_gap:
            merged[-1] = (last_s, max(last_e, e))
        else:
            merged.append((s, e))
    return merged


def subtract_intervals(domain: Interval, blocks: list[Interval]) -> list[Interval]:
    """Return domain \\ union(blocks), all half-open."""
    (D0, D1) = domain
    result: list[Interval] = []
    cur = D0
    for s, e in merge_intervals([b for b in blocks if not (b[1] <= D0 or b[0] >= D1)]):
        if s > cur:
            result.append((cur, min(s, D1)))
        cur = max(cur, e)
        if cur >= D1:
            break
    if cur < D1:
        result.append((cur, D1))
    return result


def clip_interval(iv: Interval, bounds: Interval) -> Interval | None:
    """Clip an interval to bounds; return None if empty after clipping."""
    s = max(bounds[0], iv[0])
    e = min(bounds[1], iv[1])
    return (s, e) if e > s else None


# geNomad input parsing
def _normalize(s: str) -> str:
    """Normalize a string for fuzzy matching: lowercase alphanumeric only."""
    return re.sub(r"[^a-z0-9]", "", s.lower())


def find_genomad_provirus_tsv(
    genomad_root: Path, accession: str, debug: bool = False
) -> Path | None:
    """
    Try several strategies to find <ACC>_find_proviruses/<ACC>_provirus.tsv:
    1) Exact match: {root}/{ACC}_find_proviruses/{ACC}_provirus.tsv
    2) Fuzzy: search **/*_find_proviruses/*_provirus.tsv and match on normalized ACC.
    3) Fallback: first *_find_proviruses/*_provirus.tsv under a directory that startswith ACC.
    """
    acc_norm = _normalize(accession)

    # 1) exact
    exact = (
        genomad_root
        / f"{accession}"
        / f"{accession}_find_proviruses"
        / f"{accession}_provirus.tsv"
    )
    if exact.exists():
        if debug:
            print(f"[DEBUG] Using exact TSV: {exact}")
        return exact

    # 2) fuzzy scan
    candidates = list(genomad_root.glob("**/*_find_proviruses/*_provirus.tsv"))
    for cand in candidates:
        # Parent like ".../<ACC>_find_proviruses".
        parent = cand.parent.name
        # File like "<ACC>_provirus.tsv".
        file_stem = cand.stem.replace("_provirus", "")
        # Directory just above *_find_proviruses.
        above = cand.parent.parent.name
        # Try several candidate names for normalization-based matching.
        names = [parent.replace("_find_proviruses", ""), file_stem, above]
        if any(_normalize(n) == acc_norm for n in names):
            if debug:
                print(f"[DEBUG] Fuzzy matched TSV: {cand}")
            return cand

    # 3) startswith fallback
    for cand in candidates:
        above = cand.parent.parent.name
        if above.startswith(accession):
            if debug:
                print(f"[DEBUG] Startswith fallback TSV: {cand}")
            return cand

    if debug:
        print(f"[DEBUG] No TSV found for accession '{accession}'. Looked for:")
        print(f"        - {exact}")
        print(f"        - among {len(candidates)} candidates under {genomad_root}")
    return None


def load_provirus_tsv(
    tsv_path: Path, is_one_based: bool = True
) -> dict[str, list[Interval]]:
    """
    Load geNomad <ACC>_provirus.tsv into {contig: [(start,end), ...]} in 0-based, half-open.
    geNomad typically reports 1-based inclusive; set is_one_based=True by default.
    """
    calls: dict[str, list[Interval]] = {}
    with open(tsv_path, newline="") as f:
        sample = f.read(4096)
        f.seek(0)
        # Be tolerant to delimiters.
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters="\t,;|")
        except Exception:

            class _D:
                delimiter = "\t"

            dialect = _D()
        reader = csv.DictReader(f, dialect=dialect)
        if not reader.fieldnames:
            raise ValueError(
                "Empty or headerless TSV not supported; please ensure the file has headers."
            )
        hdr = {h.lower().strip(): h for h in reader.fieldnames}

        def pick(*opts: str) -> str | None:
            for o in opts:
                if o in hdr:
                    return hdr[o]
            return None

        k_contig = pick(
            "source_seq",
            "contig",
            "sequence",
            "seqname",
            "chrom",
            "scaffold",
            "seq_name",
            "sequence_name",
        )
        k_start = pick("start", "begin", "pos_start")
        k_end = pick("end", "stop", "pos_end")
        if not (k_contig and k_start and k_end):
            raise ValueError(
                f"[{tsv_path}] Missing contig/start/end columns (headers: {reader.fieldnames})"
            )

        for row in reader:
            ctg = str(row[k_contig]).strip()
            if ctg == "" or row[k_start] is None or row[k_end] is None:
                continue
            try:
                s = int(float(str(row[k_start]).strip()))
                e = int(float(str(row[k_end]).strip()))
            except Exception:
                continue
            # Convert 1-based inclusive to 0-based half-open.
            if is_one_based:
                s = s - 1
            if e <= s:
                continue
            calls.setdefault(ctg, []).append((s, e))
    return calls


# Core logic per accession and contig
def build_bed_for_accession(
    accession: str,
    contig_lengths: dict[str, int],
    calls_by_contig: dict[str, list[Interval]],
    flank_len: int = 5000,
) -> list[tuple[str, int, int, str]]:
    """
    Return BED rows [(contig, start, end, label)] for all contigs in this accession.
    Labels: 'positive', 'near_neg', 'far_neg'.
    """
    rows: list[tuple[str, int, int, str]] = []

    for contig, L in contig_lengths.items():
        domain = (0, L)
        # Merge overlapping positive calls.
        pos = merge_intervals(calls_by_contig.get(contig, []))
        # Clamp positives to the contig domain.
        pos = [clip_interval(iv, domain) for iv in pos]
        pos = [iv for iv in pos if iv is not None]

        # Build flanking near-negative intervals beside each positive block.
        near_blocks: list[Interval] = []
        for s, e in pos:
            left = clip_interval((s - flank_len, s), domain)
            right = clip_interval((e, e + flank_len), domain)
            if left:
                near_blocks.append(left)
            if right:
                near_blocks.append(right)

        # Remove any overlap with positives, then merge again.
        if near_blocks:
            cleaned: list[Interval] = []
            for nb in near_blocks:
                remaining = subtract_intervals(nb, pos)
                cleaned.extend(remaining)
            near = merge_intervals(cleaned)
        else:
            near = []

        # Far negatives are the rest of the contig.
        union_blocks = merge_intervals(pos + near)
        far = subtract_intervals(domain, union_blocks)

        # Stage BED rows for this contig.
        contig_rows = []
        for s, e in pos:
            contig_rows.append((contig, s, e, "positive"))
        for s, e in near:
            contig_rows.append((contig, s, e, "near_neg"))
        for s, e in far:
            contig_rows.append((contig, s, e, "far_neg"))

        # Enforce non-overlap and one-base separation with priority to positives.
        contig_rows = _separate_by_one(contig_rows, L)
        rows.extend(contig_rows)

    # Sort BED rows by contig and coordinate.
    rows.sort(key=lambda r: (r[0], r[1], r[2], r[3]))
    return rows


# File-level operations
def write_bed(path: Path, rows: list[tuple[str, int, int, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        for contig, s, e, label in rows:
            f.write(f"{contig}\t{s}\t{e}\t{label}\n")

def _separate_by_one(rows_for_contig, contig_len):
    """Keep neighboring label blocks separated by one base when possible."""
    priority = {"positive": 0, "near_neg": 1, "far_neg": 2}
    # Sort by start, then prefer positives on ties.
    rows_for_contig.sort(key=lambda r: (r[1], priority.get(r[3], 9)))

    adjusted = []
    prev = None
    for contig, s, e, label in rows_for_contig:
        # Clamp to the contig domain.
        s = max(0, s)
        e = min(contig_len, e)
        if e <= s:
            continue
        if prev is None:
            prev = [contig, s, e, label]
            continue
        # Fix touching or overlapping intervals according to label priority.
        if s <= prev[2]:
            # Preserve the positive block when it conflicts with a non-positive block.
            if label == "positive" and prev[3] != "positive":
                prev[2] = min(prev[2], s - 1)
                if prev[2] > prev[1]:
                    adjusted.append(tuple(prev))
                prev = [contig, s, e, label]
            else:
                # Otherwise shift the later block right by one base.
                s = prev[2] + 1
                if s < e:
                    prev_block = tuple(prev)
                    adjusted.append(prev_block)
                    prev = [contig, s, e, label]
        else:
            adjusted.append(tuple(prev))
            prev = [contig, s, e, label]
    if prev is not None and prev[2] > prev[1]:
        adjusted.append(tuple(prev))

    # Final clamp and sort.
    out = []
    for contig, s, e, label in adjusted:
        s = max(0, s)
        e = min(contig_len, e)
        if e > s:
            out.append((contig, s, e, label))
    out.sort(key=lambda r: (r[1], priority.get(r[3], 9)))
    return out


def write_summary_csv(path: Path, rows: list[dict[str, object]]) -> None:
    header = [
        "accession",
        "total_bp",
        "positive_bp",
        "near_neg_bp",
        "far_neg_bp",
        "positive_segments",
        "near_neg_segments",
        "far_neg_segments",
        "has_positive",
        "tsv_path",
    ]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=header)
        writer.writeheader()
        for r in rows:
            writer.writerow({k: r.get(k, "") for k in header})


# Main
def main():
    ap = argparse.ArgumentParser(
        description="Generate per-accession BEDs + summary CSV from FASTAs and geNomad provirus TSVs."
    )
    ap.add_argument(
        "--fasta_dir",
        required=True,
        help="Directory with source FASTA files (e.g., raw/genomes).",
    )
    ap.add_argument(
        "--genomad_root",
        required=True,
        help="Root dir containing per-accession outputs.",
    )
    ap.add_argument("--out_root", required=True, help="Output root directory.")
    ap.add_argument(
        "--flank",
        type=int,
        default=5000,
        help="Flank length for near_neg (bp). Default: 5000.",
    )
    ap.add_argument(
        "--calls_one_based",
        action="store_true",
        help="Force geNomad TSV coordinates as 1-based inclusive (default).",
    )
    ap.add_argument(
        "--calls_zero_based",
        action="store_true",
        help="Treat geNomad TSV coordinates as 0-based half-open.",
    )
    ap.add_argument(
        "--debug",
        action="store_true",
        help="Print extra info about TSV discovery per accession.",
    )
    args = ap.parse_args()

    # Determine coordinate mode.
    is_one_based = True
    if args.calls_zero_based:
        is_one_based = False
    elif args.calls_one_based:
        is_one_based = True

    fasta_dir = Path(args.fasta_dir)
    genomad_root = Path(args.genomad_root)
    out_root = Path(args.out_root)
    out_seq = out_root / "sequences"
    out_bed = out_root / "beds"
    out_bed_pos = out_root / "beds_pos"
    out_seq.mkdir(parents=True, exist_ok=True)
    out_bed.mkdir(parents=True, exist_ok=True)
    out_bed_pos.mkdir(parents=True, exist_ok=True)

    # Scan FASTAs, gather lengths, and copy them into sequences/.
    fasta_paths: list[Path] = sorted(
        [
            p
            for p in fasta_dir.iterdir()
            if p.is_file()
            and any(
                p.name.endswith(s)
                for s in (".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz", ".fasta.gz")
            )
        ]
    )
    if not fasta_paths:
        print(f"[ERROR] No FASTA files found in {fasta_dir}", file=sys.stderr)
        sys.exit(1)

    # Build outputs accession by accession.
    summary_rows: list[dict[str, object]] = []
    for fp in fasta_paths:
        contig_lengths = read_fasta_lengths(fp)
        # Copy FASTA to sequences/.
        dest = out_seq / fp.name
        shutil.copy2(fp, dest)

        accession = infer_accession_from_filename(fp)

        # Locate the matching provirus TSV.
        prov_tsv = find_genomad_provirus_tsv(genomad_root, accession, debug=args.debug)
        if prov_tsv is None:
            print(
                f"[WARN] Missing geNomad TSV for {accession}. Tried fuzzy search under: {genomad_root}",
                file=sys.stderr,
            )
            calls_for_acc: dict[str, list[Interval]] = {}
            tsv_path_str = ""
        else:
            try:
                calls_for_acc = load_provirus_tsv(prov_tsv, is_one_based=is_one_based)
                tsv_path_str = str(prov_tsv)
                if args.debug:
                    total_iv = sum(len(v) for v in calls_for_acc.values())
                    first_keys = list(calls_for_acc.keys())[:5]
                    print(
                        f"[DEBUG] Loaded {total_iv} positive intervals from TSV. Example contigs: {first_keys}"
                    )
            except Exception as e:
                print(f"[ERROR] Failed reading {prov_tsv}: {e}", file=sys.stderr)
                calls_for_acc = {}
                tsv_path_str = str(prov_tsv)

        # Keep only contigs present in the FASTA.
        filtered_calls = {
            ctg: calls_for_acc.get(ctg, []) for ctg in contig_lengths.keys()
        }
        if args.debug:
            kept = sum(len(v) for v in filtered_calls.values())
            print(
                f"[DEBUG] After contig-name filter (matching FASTA headers), kept {kept} intervals"
            )

        rows = build_bed_for_accession(
            accession, contig_lengths, filtered_calls, flank_len=args.flank
        )

        # Write BED output.
        out_bed_path = out_bed / f"{accession}.bed"
        write_bed(out_bed_path, rows)

        # Compute per-accession summary values.
        positive_bp = sum(e - s for _, s, e, lbl in rows if lbl == "positive")
        near_neg_bp = sum(e - s for _, s, e, lbl in rows if lbl == "near_neg")
        far_neg_bp = sum(e - s for _, s, e, lbl in rows if lbl == "far_neg")
        total_bp = positive_bp + near_neg_bp + far_neg_bp
        pos_segments = sum(1 for *_, lbl in rows if lbl == "positive")
        near_segments = sum(1 for *_, lbl in rows if lbl == "near_neg")
        far_segments = sum(1 for *_, lbl in rows if lbl == "far_neg")
        has_positive = pos_segments > 0

        summary_rows.append(
            {
                "accession": accession,
                "total_bp": total_bp,
                "positive_bp": positive_bp,
                "near_neg_bp": near_neg_bp,
                "far_neg_bp": far_neg_bp,
                "positive_segments": pos_segments,
                "near_neg_segments": near_segments,
                "far_neg_segments": far_segments,
                "has_positive": int(has_positive),
                "tsv_path": tsv_path_str,
            }
        )

        # Copy positive-bearing BEDs into beds_pos/.
        if has_positive:
            shutil.copy2(out_bed_path, out_bed_pos / out_bed_path.name)
            status = "[POS]"
        else:
            status = "[NO POS]"
        print(f"[OK] {accession}: wrote {out_bed_path} ({len(rows)} rows) {status}")

    # Write the global CSV summary.
    write_summary_csv(out_root / "summary.csv", summary_rows)

    print(f"[DONE] FASTAs copied to:   {out_seq}")
    print(f"[DONE] BEDs written to:    {out_bed}")
    print(f"[DONE] POS BEDs to:        {out_bed_pos}")
    print(f"[DONE] Summary CSV at:     {out_root / 'summary.csv'}")


if __name__ == "__main__":
    main()
