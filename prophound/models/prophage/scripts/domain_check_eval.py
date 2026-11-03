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
"""Check 8 (docs/robustness-checks.md): the non-prokaryote signals of src/domain.py on a labeled panel.

Run from the model directory:
  python scripts/domain_check_eval.py --out analysis/robustness/domain_check.tsv \\
      prokaryote=../../../arc_data/lockbox2_fasta prokaryote=../../../arc_data/raw/lambda \\
      groups=<dir with the panel FASTAs and groups.json>
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from domain import genome_signals  # noqa: E402


def read_fasta(path: Path) -> dict[str, str]:
    opener = gzip.open if path.suffix == ".gz" else open
    seqs, name, buf = {}, None, []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("sources", nargs="+", help="GROUP=DIR (every .fna in DIR), or groups=DIR (DIR/groups.json names each)")
    args = ap.parse_args()
    rows = []
    for spec in args.sources:
        group, _, d = spec.partition("=")
        d = Path(d)
        labels = json.loads((d / "groups.json").read_text()) if group == "groups" else {}
        for f in sorted(d.glob("*.fna")):
            g = labels.get(f.stem, "eukaryote" if f.stem.startswith("human") else group)
            if group == "groups" and g == "bacteria_hard":
                g = "prokaryote"
            sig = genome_signals(read_fasta(f))
            rows.append({"genome": f.stem, "group": g, "source": d.name, **sig})
            print(f"{g:11s} {f.stem[:40]:40s} orf_excess {sig['orf_excess']:.3f} telomere {sig['telomere_frac']:.2f}",
                  file=sys.stderr)
    pd.DataFrame(rows).to_csv(args.out, sep="\t", index=False, float_format="%.4f")


if __name__ == "__main__":
    main()
