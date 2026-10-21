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


"""Clade-holdout affordability and clade-novel recall, from the clustering output.

Consumes `analysis/clades/elements_clades.csv` (written by
`cluster_prophage_elements.py` on the HPCC) and answers two questions locally,
without retraining (docs/scientific-validity.md, item 2a):

1. `split` -- how leaky is the existing genome-level split at the clustering
   threshold (fraction of val/test elements sharing a clade with a train
   genome), and is a clade-held-out split affordable? Genomes and clades form
   a bipartite graph; holding out whole connected components strands nothing
   by construction, so the component-size distribution (not the summary's
   `genomes_with_multiple_clades`, which counts singleton clades too) is what
   decides affordability.

2. `recall` -- covered recall on eval-sample elements that share a clade with
   a train genome ("leaky") versus those that do not ("clade-novel"), for one
   or more directories of per-genome `regions.tsv` (as written by `predict`).
   If the >=95% contract holds on the clade-novel subset, the headline recall
   is not carried by near-duplicates of training elements.

"Clade-novel" is relative to the clustering threshold (`min_jaccard`=0.2 is
roughly <=5% divergent); it does not exclude more distant homologs, so this is
a within-corpus generalization check, not the novelty claim proper.

Examples (from the repo root):
  python scripts/clade_holdout_checks.py split \
      --clades analysis/clades/elements_clades.csv --splits_dir splits
  python scripts/clade_holdout_checks.py recall \
      --clades analysis/clades/elements_clades.csv --splits_dir splits \
      --elements scripts/eval_sizestrat_elements.csv \
      --predict_dir label1=/path/to/regions_dir1 label2=/path/to/regions_dir2
"""

from __future__ import annotations

import argparse
import glob
import os
from pathlib import Path

import numpy as np
import pandas as pd
from evalkit import calls_file


def load_split_ids(splits_dir: Path, name: str) -> set:
    return {l.strip() for l in open(splits_dir / f"{name}.txt") if l.strip()}


def train_clade_set(cl: pd.DataFrame, splits_dir: Path) -> set:
    return set(cl[cl.genome_id.isin(load_split_ids(splits_dir, "train"))].clade)


def _read_clades(args: argparse.Namespace) -> pd.DataFrame:
    cl = pd.read_csv(args.clades)
    if args.clade_col != "clade":
        if args.clade_col not in cl:
            raise SystemExit(f"{args.clades} has no column {args.clade_col!r}; available: "
                             f"{[c for c in cl.columns if c.startswith('clade')]}")
        cl = cl.drop(columns=[c for c in cl.columns if c.startswith("clade") and c != args.clade_col])
        cl = cl.rename(columns={args.clade_col: "clade"})
    return cl

# ---------------------------------------------------------------- split ----
def components(cl: pd.DataFrame) -> pd.DataFrame:
    """Connected components of the genome-clade bipartite graph (union-find)."""
    parent: dict = {}

    def find(x):
        while parent.setdefault(x, x) != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for g, c in zip(cl.genome_id, cl.clade):
        parent[find(("g", g))] = find(("c", c))
    comp = pd.Series([find(("g", g)) for g in cl.genome_id], index=cl.index)
    return (cl.assign(comp=comp).groupby("comp")
            .agg(n_el=("clade", "size"), n_gen=("genome_id", "nunique"), n_cl=("clade", "nunique"))
            .sort_values("n_el", ascending=False))


def cmd_split(args: argparse.Namespace) -> None:
    cl = _read_clades(args)
    splits = Path(args.splits_dir)
    sp = {n: load_split_ids(splits, n) for n in ("train", "val", "test")}
    cl["split"] = cl.genome_id.map(lambda g: next((n for n in sp if g in sp[n]), "none"))
    print("elements per split:", cl.split.value_counts().to_dict())
    train_clades = set(cl[cl.split == "train"].clade)
    for n in ("val", "test"):
        sub = cl[cl.split == n]
        leak = sub.clade.isin(train_clades)
        print(f"[{n}] {len(sub)} elements, {leak.sum()} ({leak.mean() * 100:.1f}%) share a clade with a train genome")

    comp = components(cl)
    n_gen = cl.genome_id.nunique()
    print(f"\n[components] {len(comp)}; largest {int(comp.iloc[0].n_el)} elements / {int(comp.iloc[0].n_gen)} genomes "
          f"({comp.iloc[0].n_el / len(cl) * 100:.1f}% of elements, {comp.iloc[0].n_gen / n_gen * 100:.1f}% of genomes)")
    print("top components:\n" + comp.head(8).to_string())
    order = comp.sample(frac=1, random_state=args.seed)
    cum = order.n_el.cumsum()
    for frac in (0.10, 0.15, 0.20):
        pick = order[cum <= frac * len(cl)]
        print(f"[holdout] {frac:.0%} of elements by whole components: {int(pick.n_el.sum())} elements, "
              f"{int(pick.n_gen.sum())} genomes, {len(pick)} components, 0 stranded")
    rng = np.random.default_rng(args.seed)
    clades = cl.clade.unique()
    for frac in (0.10, 0.20):
        held = set(rng.choice(clades, int(frac * len(clades)), replace=False))
        per_g = cl.groupby("genome_id").clade.apply(lambda s: (s.isin(held).all(), s.isin(held).any()))
        test = sum(a for a, _ in per_g)
        strand = sum((not a) and b for a, b in per_g)
        print(f"[random-clade holdout] {frac:.0%} of clades: {test} genomes fully test, {strand} stranded, "
              f"{len(per_g) - test - strand} fully train")


# --------------------------------------------------------------- recall ----
def tag_elements(ev: pd.DataFrame, cl: pd.DataFrame, train_clades: set) -> pd.DataFrame:
    """Attach the best-overlapping clade row to each eval element (genome+contig+overlap)."""
    out = []
    for r in ev.itertuples():
        c = cl[(cl.genome_id == r.accession) & (cl.contig_id == r.contig)]
        ov = c[(c.end > r.start) & (c.start < r.end)]
        if len(ov) == 0:
            out.append((np.nan, False, np.nan))
            continue
        best = ov.iloc[(np.minimum(ov.end, r.end) - np.maximum(ov.start, r.start)).argmax()]
        out.append((best.clade, best.clade in train_clades, int((cl.clade == best.clade).sum())))
    return ev.assign(**dict(zip(("clade", "leaky", "clade_size"), zip(*out))))


def load_regions(d: str) -> pd.DataFrame:
    dfs = []
    for f in sorted({os.path.dirname(x) for x in glob.glob(os.path.join(d, "*", "regions.tsv"))}):
        f = str(calls_file(f))   # regions.detailed.tsv since 0.1.0
        try:
            dfs.append(pd.read_csv(f, sep="\t"))
        except pd.errors.EmptyDataError:  # header-only file: a genome with no calls
            pass
    return pd.concat(dfs, ignore_index=True)


def covered(ev: pd.DataFrame, reg: pd.DataFrame):
    cov, iou = [], []
    for r in ev.itertuples():
        g = reg[(reg.genome_id == r.accession) & (reg.contig_id == r.contig)]
        o = g[(g.end > r.start) & (g.start < r.end)]
        cov.append(len(o) > 0)
        if len(o):
            inter = np.minimum(o.end, r.end) - np.maximum(o.start, r.start)
            union = np.maximum(o.end, r.end) - np.minimum(o.start, r.start)
            iou.append(float((inter / union).max()))
        else:
            iou.append(0.0)
    return np.asarray(cov), np.asarray(iou)


def cmd_recall(args: argparse.Namespace) -> None:
    cl = _read_clades(args)
    ev = tag_elements(pd.read_csv(args.elements), cl, train_clade_set(cl, Path(args.splits_dir)))
    unmatched = ev.clade.isna().sum()
    leaky = ev.leaky.astype(bool).values
    print(f"eval elements: {len(ev)} ({unmatched} not in the clade table); leaky: {leaky.sum()}, clade-novel: {(~leaky).sum()}")
    if "stratum" in ev:
        print("leaky by stratum:", ev.groupby("stratum").leaky.sum().astype(int).to_dict())
    rows = []
    for spec in args.predict_dir:
        label, _, d = spec.rpartition("=")
        label = label or os.path.basename(d.rstrip("/"))
        cov, iou = covered(ev, load_regions(d))
        rows.append(dict(run=label,
                         recall_leaky=cov[leaky].mean(), n_leaky=int(leaky.sum()),
                         recall_novel=cov[~leaky].mean(), n_novel=int((~leaky).sum()),
                         iou_leaky=float(np.median(iou[leaky])), iou_novel=float(np.median(iou[~leaky])),
                         missed_leaky=int((~cov[leaky]).sum()), missed_novel=int((~cov[~leaky]).sum())))
        if args.list_missed:
            cols = [c for c in ("accession", "contig", "start", "end", "kb", "stratum", "leaky", "clade_size") if c in ev]
            print(f"\n[{label}] missed:\n" + ev[~cov][cols].to_string(index=False))
    pd.set_option("display.width", 200)
    print("\n" + pd.DataFrame(rows).round(3).to_string(index=False))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name, fn in (("split", cmd_split), ("recall", cmd_recall)):
        p = sub.add_parser(name)
        p.add_argument("--clades", required=True, help="elements_clades.csv from cluster_prophage_elements.py")
        p.add_argument("--splits_dir", required=True)
        p.add_argument("--seed", type=int, default=0)
        p.add_argument("--clade_col", default="clade",
                       help="which clade column to use, e.g. clade_0.05 from a --thresholds run")
        p.set_defaults(fn=fn)
    r = sub.choices["recall"]
    r.add_argument("--elements", required=True, help="eval element manifest (accession,contig,start,end[,stratum])")
    r.add_argument("--predict_dir", nargs="+", required=True, help="[label=]dir holding <genome>/regions.tsv, repeatable")
    r.add_argument("--list_missed", action="store_true")
    args = ap.parse_args()
    args.fn(args)



if __name__ == "__main__":
    main()
