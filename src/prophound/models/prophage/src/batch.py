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
"""Several genomes in one `predict` run: which inputs, and the combined outputs.

A single `--input` file keeps the single-genome layout. A batch (several `--input`, a
directory, or `--input-list`) writes `outdir/<genome_id>/` per genome with the usual
files, and at the top combined calls, `batch_summary.tsv` and a batch
`run_metadata.json`.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from .inputs import COMPRESSION_SUFFIXES, genome_id_for
from .io_checks import InputError
from .options import PredictError

# genome files picked up from a directory (with or without a compression suffix)
GENOME_SUFFIXES = (".fa", ".fna", ".fasta", ".fas", ".fsa", ".gb", ".gbk", ".gbff", ".genbank", ".embl",
                   ".gff", ".gff3")
# files a batch writes at the top of the outdir (each genome's directory holds the single-genome set)
BATCH_FILES = ("regions.tsv", "regions.filtered.tsv", "regions.bed", "regions.gff3", "batch_summary.tsv",
               "run_metadata.json")
SUMMARY_COLUMNS = ["genome_id", "input", "status", "format", "contigs", "bp", "calls", "filtered", "elapsed_s",
                   "error"]


def _is_genome_file(p: Path) -> bool:
    name = p.name.lower()
    for suffix in COMPRESSION_SUFFIXES:
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    return p.is_file() and not p.name.startswith(".") and name.endswith(GENOME_SUFFIXES)


def collect_inputs(inputs: list[str] | None, input_list: str | None) -> tuple[list[str], bool]:
    """The genome files to scan, in order, and whether this is a batch.

    Every file must exist and genome IDs must be unique; both are checked before any work.
    """
    if not inputs and not input_list:
        raise PredictError("no input: give --input (repeatable, or a directory) or --input-list")
    paths: list[str] = []
    batch = bool(input_list) or len(inputs or []) > 1
    for item in inputs or []:
        p = Path(item)
        if p.is_dir():
            batch = True
            found = sorted(str(f) for f in p.iterdir() if _is_genome_file(f))
            if not found:
                raise InputError(f"no genome files in directory {item} (looked for {', '.join(GENOME_SUFFIXES)}, "
                                 "optionally compressed)")
            paths += found
        else:
            paths.append(item)
    if input_list:
        lst = Path(input_list)
        if not lst.is_file():
            raise InputError(f"--input-list not found: {input_list}")
        for line in lst.read_text().splitlines():
            line = line.split("#", 1)[0].strip()
            if line:
                p = Path(line).expanduser()
                paths.append(str(p if p.is_absolute() else lst.parent / p))
        if not paths:
            raise InputError(f"--input-list {input_list} names no files")
    missing = [p for p in paths if not Path(p).exists()]
    if missing:
        raise InputError(f"input not found: {', '.join(missing[:5])}" + (" ..." if len(missing) > 5 else ""))
    dup = sorted(g for g, n in Counter(genome_id_for(p) for p in paths).items() if n > 1)
    if dup:
        raise InputError(f"several inputs share a genome ID ({', '.join(dup[:5])}); each genome's outputs go to "
                         "outdir/<genome_id>, so rename one of each")
    return paths, batch


def batch_outputs_held(outdir: Path, genome_ids: list[str], held) -> list[str]:
    """Batch outputs already in `outdir`: top-level files, and genome directories that hold
    outputs (`held(dir)` lists those)."""
    if not outdir.is_dir():
        return []
    top = [name for name in BATCH_FILES if (outdir / name).exists()]
    return top + [f"{g}/" for g in genome_ids if held(outdir / g)]


@dataclass
class GenomeRow:
    """One line of batch_summary.tsv."""

    genome_id: str
    input: str
    status: str
    format: str = ""
    contigs: int | None = None
    bp: int | None = None
    calls: int | None = None
    filtered: int | None = None
    elapsed_s: float | None = None
    error: str = ""


def combined_calls(tables: list[pd.DataFrame]) -> pd.DataFrame:
    """Every genome's calls in input order; genome_id tells them apart, call_id is per genome."""
    tables = [t for t in tables if len(t)]
    return pd.concat(tables, ignore_index=True) if tables else pd.DataFrame(columns=["call_id", "genome_id",
                                                                                      "contig_id", "start", "end"])


def unique_names(calls: pd.DataFrame) -> pd.DataFrame:
    """For combined BED/GFF3, where names must be unique: call_id becomes `genome_id:call_id`."""
    out = calls.copy()
    if len(out):
        out["call_id"] = out["genome_id"].astype(str) + ":" + out["call_id"].astype(str)
    return out


def summary_table(rows: list[GenomeRow]) -> pd.DataFrame:
    table = pd.DataFrame([vars(r) for r in rows], columns=SUMMARY_COLUMNS)
    # counts stay integers (blank for a failed genome) instead of turning into floats
    return table.astype(dict.fromkeys(("contigs", "bp", "calls", "filtered"), "Int64"))
