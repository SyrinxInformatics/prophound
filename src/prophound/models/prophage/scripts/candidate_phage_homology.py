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
"""Phage homology of off-tool candidates on held-out genomes (docs/candidate-phage-homology.md).

  build  candidate, positive-control and negative-control regions from a run's calls and truth
  score  apply the preregistered protein call (scripts/sakai_phage_homology.py's rule) to DIAMOND
         blastx output, with a genome-bootstrap interval for candidates minus negative controls

  python scripts/candidate_phage_homology.py build --runs analysis/lockbox2/prophound \\
      --elements analysis/lockbox2/lockbox_elements.csv --genome_list analysis/lockbox2/lockbox_genome_list.txt \\
      --fasta_dir ../../../../arc_data/lockbox2_fasta --out analysis/candidate_homology
  diamond blastx ... -q <out>/regions.fna -o <out>/diamond.tsv   (settings in the doc)
  python scripts/candidate_phage_homology.py score --out analysis/candidate_homology
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sakai_phage_homology import PHAGE_SPECIFIC, PROTEIN_MIN_LOCI, SHARED, intervals, merged, organism, product

SEED = 20261027
PAD = 5000


def read_fasta(path: Path) -> dict[str, str]:
    seqs, name, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def frac_inside(s: int, e: int, ivs: list[tuple[int, int]]) -> float:
    pts = sorted((max(s, a), min(e, b)) for a, b in ivs if b > s and a < e)
    tot, cur = 0, None
    for a, b in pts:
        if cur and a <= cur[1]:
            cur = (cur[0], max(cur[1], b))
        else:
            if cur:
                tot += cur[1] - cur[0]
            cur = (a, b)
    if cur:
        tot += cur[1] - cur[0]
    return tot / max(1, e - s)


def cmd_build(a) -> None:
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    truth = pd.read_csv(a.elements)
    genomes = [l.strip() for l in open(a.genome_list) if l.strip()]
    rng = np.random.default_rng(SEED)
    rows, missing_neg = [], 0
    with open(out / "regions.fna", "w") as fa:
        for g in genomes:
            seqs = read_fasta(Path(a.fasta_dir) / f"{g}.fna")
            calls = {m: pd.read_csv(Path(a.runs) / m / g / "regions.detailed.tsv", sep="\t") for m in ("default", "exploratory", "conservative")}
            tr = truth[truth.accession == g]
            busy = {c: [(int(r.start), int(r.end)) for r in tr[tr.contig == c].itertuples()] for c in seqs}
            for m, d in calls.items():
                for r in d.itertuples():
                    busy.setdefault(r.contig_id, []).append((int(r.start), int(r.end)))
            def emit(rid, kind, c, s, e, extra):
                fa.write(f">{rid}|{c}:{s+1}-{e}\n{seqs[c][s:e]}\n")
                rows.append(dict(region=rid, kind=kind, genome=g, contig=c, start=s, end=e, kb=round((e - s) / 1000, 1)) | extra)
            dflt = calls["default"]
            for r in dflt.itertuples():
                ivs = [(int(t.start), int(t.end)) for t in tr[tr.contig == r.contig_id].itertuples()]
                inside_truth = frac_inside(int(r.start), int(r.end), ivs)
                cover_el = max((frac_inside(int(t.start), int(t.end), [(int(r.start), int(r.end))]) for t in tr[tr.contig == r.contig_id].itertuples()), default=0)
                extra = dict(support=int(r.seed_support), boundary=r.boundary, max_prob=float(r.max_prob))
                if inside_truth == 0:
                    rid = f"cand_d_{len(rows)}"
                    emit(rid, "candidate_default", r.contig_id, int(r.start), int(r.end), extra)
                    n = int(r.end - r.start)
                    order = [r.contig_id] + [c for c in seqs if c != r.contig_id]
                    placed = False
                    for c in order:
                        L = len(seqs[c])
                        if L < n + 2 * PAD:
                            continue
                        for _ in range(2000):
                            s = int(rng.integers(0, L - n)); e = s + n
                            if all(e + PAD <= x or s - PAD >= y for x, y in busy.get(c, [])):
                                busy.setdefault(c, []).append((s, e))
                                emit(f"neg_{rid}", "negative_control", c, s, e, dict(paired=rid)); placed = True
                                break
                        if placed:
                            break
                    missing_neg += not placed
                elif cover_el >= 0.9 and inside_truth >= 0.9:
                    emit(f"pos_{len(rows)}", "positive_control", r.contig_id, int(r.start), int(r.end), extra)
            dl = [(r.contig_id, int(r.start), int(r.end)) for r in dflt.itertuples()]
            for r in calls["exploratory"].itertuples():
                ivs = [(int(t.start), int(t.end)) for t in tr[tr.contig == r.contig_id].itertuples()]
                if frac_inside(int(r.start), int(r.end), ivs) > 0:
                    continue
                if any(c == r.contig_id and e > r.start and s < r.end for c, s, e in dl):
                    continue
                emit(f"cand_e_{len(rows)}", "candidate_exploratory_only", r.contig_id, int(r.start), int(r.end),
                     dict(support=int(r.seed_support), boundary=r.boundary, max_prob=float(r.max_prob)))
    t = pd.DataFrame(rows); t.to_csv(out / "regions.tsv", sep="\t", index=False)
    s = t.groupby("kind").agg(n=("region", "size"), Mb=("kb", lambda x: round(x.sum() / 1000, 2)), genomes=("genome", "nunique"))
    print(s.to_string()); print(f"negative controls that could not be placed: {missing_neg}")


def cmd_score(a) -> None:
    out = Path(a.out)
    reg = pd.read_csv(out / "regions.tsv", sep="\t")
    cols = ["qseqid", "sseqid", "pident", "length", "qstart", "qend", "evalue", "bitscore", "stitle"]
    hits = intervals(pd.read_csv(out / "diamond.tsv", sep="\t", names=cols))
    hits["region"] = hits.qseqid.str.split("|").str[0]
    res = []
    for r in reg.itertuples():
        h = hits[hits.region == r.region]
        loci = merged(h) if len(h) else []
        spec, shared, names = 0, 0, []
        for lo, hi, g in loci:
            best = g.loc[g.bitscore.idxmax()]
            p = product(best.stitle)
            if PHAGE_SPECIFIC.search(p):
                spec += 1; names.append(f"{p} ({organism(best.stitle)}, {best.pident:.0f}%)")
            elif SHARED.search(p):
                shared += 1
        res.append(dict(region=r.region, protein_cov=round(sum(hi - lo + 1 for lo, hi, _ in loci) / ((r.end - r.start) or 1), 3),
                        loci=len(loci), phage_specific_loci=spec, integrase_transposase_loci=shared,
                        protein_pass=spec >= PROTEIN_MIN_LOCI, phage_specific_hits="; ".join(names)))
    t = reg.merge(pd.DataFrame(res), on="region")
    t.to_csv(out / "scored.tsv", sep="\t", index=False)
    summ = t.groupby("kind").agg(n=("region", "size"), passed=("protein_pass", "sum"), rate=("protein_pass", "mean"),
                                 median_cov=("protein_cov", "median"))
    # primary: candidates minus their negative controls, genome bootstrap
    c = t[t.kind == "candidate_default"]; n = t[t.kind == "negative_control"]
    by_g = pd.DataFrame({"c_pass": c.groupby("genome").protein_pass.sum(), "c_n": c.groupby("genome").size(),
                         "n_pass": n.groupby("genome").protein_pass.sum(), "n_n": n.groupby("genome").size()}).fillna(0)
    diff = c.protein_pass.mean() - n.protein_pass.mean()
    rng = np.random.default_rng(20261027); B = []
    gs = by_g.index.to_numpy()
    for _ in range(10000):
        s = by_g.loc[rng.choice(gs, len(gs))]
        if s.c_n.sum() and s.n_n.sum():
            B.append(s.c_pass.sum() / s.c_n.sum() - s.n_pass.sum() / s.n_n.sum())
    lo, hi = np.percentile(B, [2.5, 97.5])
    primary = dict(candidate_rate=float(c.protein_pass.mean()), negative_rate=float(n.protein_pass.mean()),
                   difference=float(diff), genome_boot95=[float(lo), float(hi)], n_candidates=len(c), n_negatives=len(n),
                   n_genomes_with_candidates=int(c.genome.nunique()))
    secondary = {}
    for col, bins in [("support", None), ("boundary", None)]:
        secondary[col] = c.groupby(col).protein_pass.agg(["size", "sum"]).to_dict(orient="index")
    c = c.assign(size=pd.cut(c.kb, [0, 5, 10, 15, 30, 60, 1e9], labels=["<5", "5-10", "10-15", "15-30", "30-60", "60+"]))
    secondary["size_kb"] = {str(k): v for k, v in c.groupby("size", observed=True).protein_pass.agg(["size", "sum"]).to_dict(orient="index").items()}
    json.dump(dict(primary=primary, by_set=summ.reset_index().to_dict(orient="records"), secondary=secondary),
              open(out / "results.json", "w"), indent=2, default=float)
    with pd.option_context("display.width", 220, "display.max_colwidth", 150):
        print(summ.to_string()); print(); print(json.dumps(primary, indent=1)); print(json.dumps(secondary, indent=1, default=float))
        p = t[(t.kind.str.startswith("candidate")) & t.protein_pass]
        print(f"\npassing candidates: {len(p)}")
        print(p[["region", "kind", "genome", "contig", "start", "end", "kb", "support", "boundary", "phage_specific_loci", "phage_specific_hits"]].to_string(index=False))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build")
    b.add_argument("--runs", required=True); b.add_argument("--elements", required=True)
    b.add_argument("--genome_list", required=True); b.add_argument("--fasta_dir", required=True); b.add_argument("--out", required=True)
    s = sub.add_parser("score"); s.add_argument("--out", required=True)
    a = ap.parse_args()
    {"build": cmd_build, "score": cmd_score}[a.cmd](a)


if __name__ == "__main__":
    main()
