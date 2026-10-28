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
"""Golden regression check: a fixed small genome set against stored regions.tsv files.

The genomes are in tests/golden/genomes.tsv (NCBI URL and SHA-256 of the uncompressed FASTA);
the stored outputs, made on CPU with the 0.1.0 defaults, are in tests/golden/expected/<genome>/.

  fetch  download the genomes into --fasta_dir and check their SHA-256
  make   run prophound (Default, CPU) and store regions.tsv and regions.detailed.tsv as golden
  check  run prophound and compare: byte for byte on CPU; the calls only on any other device,
         since devices agree in calls, not in every digit (docs/benchmarks.md)

  pixi run python scripts/golden_check.py fetch --fasta_dir <dir>
  pixi run python scripts/golden_check.py check --fasta_dir <dir> --out <dir> [--device cpu]
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

import pandas as pd

GOLDEN = Path(__file__).resolve().parents[1] / "tests" / "golden"
FILES = ["regions.tsv", "regions.detailed.tsv"]
CALL_COLUMNS = ["call_id", "contig_id", "start", "end", "boundary", "contig_edge", "seed_support"]


def genomes() -> pd.DataFrame:
    return pd.read_csv(GOLDEN / "genomes.tsv", sep="\t")


def fasta(row, fasta_dir: Path) -> Path:
    return fasta_dir / f"{row.genome_id}.fna"


def cmd_fetch(a) -> None:
    d = Path(a.fasta_dir); d.mkdir(parents=True, exist_ok=True)
    for row in genomes().itertuples():
        p = fasta(row, d)
        if not p.exists():
            with urllib.request.urlopen(row.fna_url) as r, open(p, "wb") as f:
                shutil.copyfileobj(gzip.GzipFile(fileobj=r), f)
        got = hashlib.sha256(p.read_bytes()).hexdigest()
        if got != row.sha256:
            sys.exit(f"[golden] FAIL: {p} sha256 {got} != {row.sha256}")
        print(f"[golden] {row.genome_id}: sha256 ok")


def predict(row, fasta_dir: Path, out: Path, device: str) -> Path:
    o = out / row.genome_id
    subprocess.run([sys.executable, "-c", "from prophound.cli import main; main()", "predict",
                    "--input", str(fasta(row, fasta_dir)), "--outdir", str(o), "--device", device, "-q", "--force"],
                   check=True, capture_output=True, text=True)
    return o


def cmd_make(a) -> None:
    for row in genomes().itertuples():
        o = predict(row, Path(a.fasta_dir), Path(a.out), "cpu")
        dst = GOLDEN / "expected" / row.genome_id; dst.mkdir(parents=True, exist_ok=True)
        for f in FILES:
            shutil.copyfile(o / f, dst / f)
        print(f"[golden] {row.genome_id}: stored")


def cmd_check(a) -> None:
    fail = False
    for row in genomes().itertuples():
        o = predict(row, Path(a.fasta_dir), Path(a.out), a.device)
        exp = GOLDEN / "expected" / row.genome_id
        if a.device == "cpu":
            bad = [f for f in FILES if (o / f).read_bytes() != (exp / f).read_bytes()]
            ok, what = not bad, "byte-identical" if not bad else f"differs in {', '.join(bad)}"
        else:
            g, e = (pd.read_csv(p / "regions.detailed.tsv", sep="\t")[CALL_COLUMNS] for p in (o, exp))
            ok = g.equals(e); what = "calls identical" if ok else "calls differ"
        fail |= not ok
        print(f"[golden] {row.genome_id}: {what}")
    print("[golden] FAIL" if fail else "[golden] PASS")
    sys.exit(1 if fail else 0)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    f = sub.add_parser("fetch"); f.add_argument("--fasta_dir", required=True)
    for n in ("make", "check"):
        p = sub.add_parser(n); p.add_argument("--fasta_dir", required=True); p.add_argument("--out", required=True)
        if n == "check":
            p.add_argument("--device", default="cpu")
    a = ap.parse_args()
    dict(fetch=cmd_fetch, make=cmd_make, check=cmd_check)[a.cmd](a)


if __name__ == "__main__":
    main()
