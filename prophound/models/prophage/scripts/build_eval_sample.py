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




"""Build a held-out, element-size-stratified evaluation sample from geNomad calls.

Written 2026-09-20. The original 18-genome sample (`eval_genome_list.txt`) was
stratified by *genome* size, which left element sizes to chance: 1 of its 52
regions is under 10 kb and 5 are under 15 kb. Every coverage-builder setting
measured so far (`--region-builder coverage`, `--coverage-frac`, operating mode,
replicate ensembling) was selected on that same 52-region set, so none of it is
held out, and the earlier claim that elements below
~10 kb are structurally uncallable has never been testable on it.

This script selects a second sample that fixes both problems:

  * **Held out** -- drawn only from test-split genomes (`splits/test.txt`, never
    trained on) that are NOT in `eval_genome_list.txt`, so builder settings tuned
    on the original 52 can be confirmed on regions that had no part in choosing
    them.
  * **Stratified by element length**, taking every available element in the
    scarce small strata rather than sampling them.

Note on what geNomad can and cannot supply (measured over all 3,970 elements in
`genomad_out/`): the smallest call in the whole dataset is 3.50 kb, and only
2.9% of elements are under 10 kb. A <5 kb stratum is empty in the held-out pool,
so the smallest testable class here is 5-8 kb. Evaluating anything below that
needs a different truth source (PHASTEST -- analysis/phastest_cmp/, or the
literature-curated coordinates in analysis/litpilot/), not a larger sample.

Selection is deterministic: scarce strata are taken whole, then genomes are
added by a greedy rule that maximises elements contributed to still-unfilled
strata, ties broken by accession so reruns are identical.

Usage:
    python build_eval_sample.py --genomad_dir /path/to/genomad_out \
        --out_list eval_genome_list_sizestrat.txt \
        --out_manifest eval_sizestrat_elements.csv
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
PROPHAGE = HERE.parent

# (label, low_kb, high_kb, target_elements). `None` target = take everything
# available: the small strata are too scarce to sample from.
STRATA = [
    ("5-8kb", 5.0, 8.0, None),
    ("8-10kb", 8.0, 10.0, None),
    ("10-15kb", 10.0, 15.0, 30),
    ("15-30kb", 15.0, 30.0, 30),
    ("30-60kb", 30.0, 60.0, 35),
    ("60kb+", 60.0, float("inf"), 20),
]


def load_all_elements(genomad_dir: Path) -> pd.DataFrame:
    rows = []
    for d in sorted(genomad_dir.iterdir()):
        if not d.is_dir():
            continue
        tsv = d / f"{d.name}_find_proviruses" / f"{d.name}_provirus.tsv"
        if not tsv.exists() or tsv.stat().st_size < 10:
            continue
        try:
            df = pd.read_csv(tsv, sep="\t")
        except Exception:
            continue
        for _, r in df.iterrows():
            rows.append(
                (d.name, r["source_seq"], int(r["start"]), int(r["end"]),
                 int(r["end"]) - int(r["start"]), int(r.get("n_genes", 0)))
            )
    e = pd.DataFrame(rows, columns=["accession", "contig", "start", "end", "len_bp", "n_genes"])
    e["kb"] = e.len_bp / 1000.0
    return e


def label_stratum(kb: float) -> str | None:
    for label, lo, hi, _ in STRATA:
        if lo <= kb < hi:
            return label
    return None  # below the smallest stratum (nothing in the held-out pool)


def select(pool: pd.DataFrame, max_genomes: int) -> list:
    """Greedy: scarce strata whole, then fill the rest by marginal contribution."""
    targets = {lab: t for lab, _, _, t in STRATA}
    avail = {lab: int((pool.stratum == lab).sum()) for lab, _, _, _ in STRATA}
    # `None` target means take every element in that stratum.
    need = {lab: (avail[lab] if targets[lab] is None else min(targets[lab], avail[lab]))
            for lab in targets}

    chosen: list = []
    got = dict.fromkeys(targets, 0)

    # Phase 1 -- every genome holding an element in a take-everything stratum.
    must = sorted(pool[pool.stratum.isin([l for l, _, _, t in STRATA if t is None])].accession.unique())
    for acc in must:
        chosen.append(acc)
    for lab in targets:
        got[lab] = int(pool[pool.accession.isin(chosen) & (pool.stratum == lab)].shape[0])

    # Phase 2 -- greedy fill of the remaining strata.
    remaining = sorted(set(pool.accession.unique()) - set(chosen))
    while len(chosen) < max_genomes:
        deficits = {lab: max(0, need[lab] - got[lab]) for lab in targets}
        if not any(deficits.values()):
            break
        best, best_gain = None, 0
        for acc in remaining:
            sub = pool[pool.accession == acc]
            gain = sum(min(int((sub.stratum == lab).sum()), deficits[lab]) for lab in targets)
            if gain > best_gain:
                best, best_gain = acc, gain
        if best is None:
            break
        chosen.append(best)
        remaining.remove(best)
        for lab in targets:
            got[lab] += int(pool[(pool.accession == best) & (pool.stratum == lab)].shape[0])
    return sorted(chosen)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genomad_dir", required=True, help="Directory of per-genome geNomad output dirs.")
    ap.add_argument("--splits_test", default=str(PROPHAGE / "splits/test.txt"))
    ap.add_argument("--exclude_list", default=str(HERE / "eval_genome_list.txt"),
                    help="Genomes to hold out of the pool (the original eval sample).")
    ap.add_argument("--out_list", default=str(HERE / "eval_genome_list_sizestrat.txt"))
    ap.add_argument("--out_manifest", default=str(HERE / "eval_sizestrat_elements.csv"))
    ap.add_argument("--max_genomes", type=int, default=45)
    args = ap.parse_args()

    e = load_all_elements(Path(args.genomad_dir))
    test = {l.strip() for l in open(args.splits_test) if l.strip()}
    used = {l.strip() for l in open(args.exclude_list) if l.strip()}
    e["stratum"] = e.kb.map(label_stratum)

    pool = e[e.accession.isin(test) & ~e.accession.isin(used) & e.stratum.notna()].copy()
    print(f"pool: {pool.accession.nunique()} held-out genomes, {len(pool)} elements")

    chosen = select(pool, args.max_genomes)
    sel = pool[pool.accession.isin(chosen)].copy().sort_values(["accession", "contig", "start"])

    print(f"\nselected {len(chosen)} genomes, {len(sel)} elements")
    print(f"{'stratum':>10}{'selected':>10}{'available':>11}{'current-18':>12}")
    cur = e[e.accession.isin(used)]
    for lab, lo, hi, t in STRATA:
        print(f"{lab:>10}{int((sel.stratum == lab).sum()):>10}"
              f"{int((pool.stratum == lab).sum()):>11}{int((cur.stratum == lab).sum()):>12}")

    Path(args.out_list).write_text("\n".join(chosen) + "\n")
    sel.to_csv(args.out_manifest, index=False)
    print(f"\nwrote {args.out_list}\nwrote {args.out_manifest}")


if __name__ == "__main__":
    main()
