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
"""Scores the Sakai phage-homology tests exactly as preregistered (docs/sakai-candidates.md).

Protein test (DIAMOND blastx against NCBI RefSeq viral proteins): per region, the fraction
of its bases covered by hits, and the number of distinct loci (merged overlapping hit
segments) whose best hit's product is phage-specific (head and packaging, tail, lysis). A
region passes with at least 2 such loci.

Nucleotide test (megablast against core_nt viruses): per region, the fraction of its bases
covered by virus hits. A region passes at 10% or more.

Query names are `<id>|<contig>:<start>-<end>` (1-based, as written by the query builder).

  python scripts/sakai_phage_homology.py --dir analysis/sakai_novel
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

SETS = {"default_novel": "default", "q_exploratory_only": "exploratory_only",
        "q_positive_control": "positive_control", "q_negative_control": "negative_control"}
PHAGE_SPECIFIC = re.compile(
    r"capsid|portal|terminase|prohead|scaffold|major head|head-tail|connector|neck|"
    r"tail|baseplate|tape measure|sheath|tube|holin|endolysin|lysozyme|spanin|\bRz\b", re.I)
SHARED = re.compile(r"integrase|recombinase|transposase|resolvase|invertase", re.I)
PROTEIN_MIN_LOCI = 2
NUCLEOTIDE_MIN_FRAC = 0.10


def region_lengths(fna: Path) -> dict[str, int]:
    out = {}
    for line in open(fna):
        if line.startswith(">"):
            q = line[1:].split()[0]
            a, b = q.rsplit(":", 1)[1].split("-")
            out[q] = int(b) - int(a) + 1
    return out


def intervals(df: pd.DataFrame) -> pd.DataFrame:
    lo = df[["qstart", "qend"]].min(axis=1)
    hi = df[["qstart", "qend"]].max(axis=1)
    return df.assign(lo=lo, hi=hi)


def merged(df: pd.DataFrame) -> list[tuple[int, int, pd.DataFrame]]:
    loci, cur, rows = [], None, []
    for r in df.sort_values("lo").itertuples():
        if cur and r.lo <= cur[1]:
            cur = (cur[0], max(cur[1], r.hi)); rows.append(r.Index)
        else:
            if cur:
                loci.append((cur[0], cur[1], df.loc[rows]))
            cur, rows = (r.lo, r.hi), [r.Index]
    if cur:
        loci.append((cur[0], cur[1], df.loc[rows]))
    return loci


def product(stitle: str) -> str:
    t = re.sub(r"^\S+\s+", "", str(stitle))
    return re.sub(r"\s*\[[^\]]*\]\s*$", "", t)


def organism(stitle: str) -> str:
    m = re.search(r"\[([^\]]*)\]\s*$", str(stitle))
    return m.group(1) if m else ""


def score_protein(d: Path) -> pd.DataFrame:
    cols = ["qseqid", "sseqid", "pident", "length", "qstart", "qend", "evalue", "bitscore", "stitle"]
    rows = []
    for stem, setname in SETS.items():
        lengths = region_lengths(d / f"{stem}.fna")
        f = d / f"diamond_{stem}.tsv"
        hits = pd.read_csv(f, sep="\t", names=cols) if f.stat().st_size else pd.DataFrame(columns=cols)
        hits = intervals(hits)
        for q, L in lengths.items():
            h = hits[hits.qseqid == q]
            loci = merged(h) if len(h) else []
            cov = sum(b - a + 1 for a, b, _ in loci) / L
            spec, shared, names = 0, 0, []
            for a, b, g in loci:
                best = g.loc[g.bitscore.idxmax()]
                p = product(best.stitle)
                if PHAGE_SPECIFIC.search(p):
                    spec += 1; names.append(f"{p} ({organism(best.stitle)}, {best.pident:.0f}%)")
                elif SHARED.search(p):
                    shared += 1
            rows.append(dict(set=setname, region=q.split("|")[0], where=q.split("|")[1], kb=round(L / 1000, 1),
                             protein_cov=round(cov, 3), loci=len(loci), phage_specific_loci=spec,
                             integrase_transposase_loci=shared, protein_pass=spec >= PROTEIN_MIN_LOCI,
                             phage_specific_hits="; ".join(names)))
    return pd.DataFrame(rows)


def score_nucleotide(d: Path) -> pd.DataFrame | None:
    files = {"default_novel": d / "blast_viruses.tsv", "q_exploratory_only": d / "blast_viruses_exploratory_only.tsv",
             "q_positive_control": d / "blast_viruses_positive_control.tsv",
             "q_negative_control": d / "blast_viruses_negative_control.tsv"}
    if not all(f.exists() for f in files.values()):
        return None
    cols = ["qseqid", "sacc", "pident", "length", "qstart", "qend", "evalue", "bitscore", "staxids", "sscinames", "stitle"]
    rows = []
    for stem, f in files.items():
        lengths = region_lengths(d / f"{stem}.fna")
        hits = pd.read_csv(f, sep="\t", names=cols) if f.stat().st_size else pd.DataFrame(columns=cols)
        hits = intervals(hits)
        for q, L in lengths.items():
            h = hits[hits.qseqid == q]
            loci = merged(h) if len(h) else []
            cov = sum(b - a + 1 for a, b, _ in loci) / L
            top = h.loc[h.length.idxmax()] if len(h) else None
            rows.append(dict(set=SETS[stem], region=q.split("|")[0], nt_cov=round(cov, 3),
                             longest_aln=int(top.length) if top is not None else 0,
                             longest_pident=round(float(top.pident), 1) if top is not None else None,
                             longest_hit=str(top.stitle)[:90] if top is not None else "",
                             nt_pass=cov >= NUCLEOTIDE_MIN_FRAC))
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", required=True, type=Path)
    a = ap.parse_args()
    if not (a.dir / "default_novel.fna").exists():
        raise SystemExit(f"no default_novel.fna in {a.dir}")
    # the per-set query files share the region names of controls_and_exploratory.fna
    p = score_protein(a.dir)
    n = score_nucleotide(a.dir)
    t = p.merge(n, on=["set", "region"], how="left") if n is not None else p
    t.to_csv(a.dir / "phage_homology.regions.tsv", sep="\t", index=False)
    agg = {"n": ("region", "size"), "protein_pass": ("protein_pass", "sum"), "median_protein_cov": ("protein_cov", "median")}
    if n is not None:
        agg |= {"nt_pass": ("nt_pass", "sum"), "median_nt_cov": ("nt_cov", "median")}
    s = t.groupby("set", sort=False).agg(**agg)
    s.to_csv(a.dir / "phage_homology.summary.tsv", sep="\t")
    with pd.option_context("display.width", 250, "display.max_colwidth", 160):
        print(s.to_string()); print()
        print(t.drop(columns=[c for c in ("longest_hit",) if c in t.columns]).to_string(index=False))


if __name__ == "__main__":
    main()
