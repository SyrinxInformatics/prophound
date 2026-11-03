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
"""Inputs for the independent confirmation run (analysis/confirm18/).

1. Element table for the original 18-genome sample (eval_genome_list.txt), in
   the same format as eval_sizestrat_elements.csv, from eval_genomad_truth/.
2. A seeded sample of test-split genomes where geNomad called no provirus, for
   a call count on (geNomad-)phage-free sequence. Calls there are either false
   positives or candidates geNomad missed; this sample cannot tell which.
   Stratified by genome-size tertile so small and large genomes are both present.

Example:
  python scripts/build_confirm_inputs.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_eval_sample import label_stratum  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def elements_18(genome_list: Path, truth_dir: Path) -> pd.DataFrame:
    rows = []
    for g in genome_list.read_text().split():
        tsv = truth_dir / g / f"{g}_find_proviruses" / f"{g}_provirus.tsv"
        df = pd.read_csv(tsv, sep="\t")
        for r in df.itertuples():
            rows.append((g, r.source_seq, int(r.start), int(r.end), int(r.end) - int(r.start), int(r.n_genes)))
    e = pd.DataFrame(rows, columns=["accession", "contig", "start", "end", "len_bp", "n_genes"])
    e["kb"] = e.len_bp / 1000.0
    e["stratum"] = [label_stratum(k) or "<5kb" for k in e.kb]
    return e


def phagefree_sample(summary: Path, test_split: Path, exclude: set, n: int, seed: int) -> list:
    s = pd.read_csv(summary)
    test = set(test_split.read_text().split())
    pool = s[(s.has_positive == 0) & s.accession.isin(test) & ~s.accession.isin(exclude)].copy()
    pool["tertile"] = pd.qcut(pool.total_bp, 3, labels=False)
    per = [n // 3 + (1 if i < n % 3 else 0) for i in range(3)]
    picks = [pool[pool.tertile == t].sample(k, random_state=seed) for t, k in enumerate(per)]
    return sorted(pd.concat(picks).accession)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n_phagefree", type=int, default=20)
    ap.add_argument("--seed", type=int, default=20260928)
    args = ap.parse_args()

    e = elements_18(HERE / "eval_genome_list.txt", HERE / "eval_genomad_truth")
    e.to_csv(HERE / "eval_confirm18_elements.csv", index=False)
    print(f"[done] {len(e)} elements in {e.accession.nunique()} genomes; strata {e.stratum.value_counts().to_dict()}")

    used = set((HERE / "eval_genome_list.txt").read_text().split()) | \
        set((HERE / "eval_genome_list_sizestrat.txt").read_text().split())
    pf = phagefree_sample(ROOT / "summary.csv", ROOT / "splits/test.txt", used, args.n_phagefree, args.seed)
    (HERE / "eval_phagefree_list.txt").write_text("\n".join(pf) + "\n")
    print(f"[done] {len(pf)} phage-free genomes -> eval_phagefree_list.txt")


if __name__ == "__main__":
    main()
