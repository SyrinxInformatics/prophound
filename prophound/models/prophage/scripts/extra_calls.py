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
"""What prophound's extra calls are made of (docs/extra-calls.md), on the second lockbox.

Every region (off-geNomad Default calls, Exploratory-only calls, length- and composition-matched
host stretches) gets the multi-label signal table and one exclusive bin, in the preregistered
order: 0 known to geNomad (raw provirus call), 1 phage-like (homology v2.1), 2 other mobile
element (CONJscan, plasmid or virus hallmark, integrase-family description, att boundary),
3 island of unknown function (tetranucleotide atypicality >= the genome's 95th percentile for
windows of that length), 4 no detectable signal. Shares carry 95% genome-bootstrap intervals.
The open look ranks gene descriptions and themes by enrichment over the length-matched stretches.

  python scripts/extra_calls.py --v21 analysis/candidate_homology/v2_1/scored.tsv \\
      --genomad ../../../arc_data/genomad_out_lockbox2 --fasta_dir ../../../arc_data/lockbox2_fasta \\
      --out analysis/extra_calls
"""

from __future__ import annotations

import argparse
import glob
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from candidate_phage_homology_v2 import jsd, kmer_index, profile, read_fasta

SEED = 20261028
N_WIN = 50
MOBILE = re.compile(r"\b(integrase|recombinase|transposase|resolvase|invertase)\b", re.I)
THEMES = {
    "biosynthetic cluster": r"polyketide|nonribosomal|\bNRPS\b|\bPKS\b|lantibiotic|bacteriocin",
    "defense": r"restriction|\bCRISPR\b|\bCas\d|abortive|antitoxin|\btoxin\b|retron|\bBREX\b|DISARM",
    "secretion or pilus": r"secretion|\bpil[iu]s?\b|\bpil[A-Z]\b|type IV|type VI|\bT[46]SS\b|fimbri",
    "surface polysaccharide": r"glycosyl ?transferase|polysaccharide|capsul|O-antigen|lipopolysaccharide|\bwz[a-z]\b",
    "metabolism": r"dehydrogenase|reductase|oxidase|synthase|kinase|transporter|permease",
    "phage-like": r"phage|capsid|portal|terminase|tail|holin|lysozyme",
}
KINDS = ["candidate_default", "candidate_exploratory_only", "negative_control", "composition_control"]


def boot_share(df: pd.DataFrame, col: str, B: int = 2000) -> list:
    g = df.groupby("genome")[col].agg(["sum", "size"])
    rng = np.random.default_rng(20261027); out = []
    for _ in range(B):
        s = g.iloc[rng.integers(0, len(g), len(g))]
        out.append(s["sum"].sum() / s["size"].sum())
    return [round(float(x), 3) for x in np.percentile(out, [2.5, 97.5])]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--v21", required=True); ap.add_argument("--genomad", required=True)
    ap.add_argument("--fasta_dir", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    reg = pd.read_csv(a.v21, sep="\t", low_memory=False)
    reg = reg[reg.kind.isin(KINDS) & (reg.start >= 0)].copy()
    reg["start"] = reg.start.astype(int); reg["end"] = reg.end.astype(int)
    rng = np.random.default_rng(SEED)
    rows, desc_rows = [], []
    for g, rg in reg.groupby("genome"):
        gd = Path(a.genomad) / g
        genes = pd.read_csv(next(gd.glob("*_annotate/*_genes.tsv")), sep="\t")
        genes["contig"] = genes.gene.str.rsplit("_", n=1).str[0]
        genes["mid"] = (genes.start + genes.end) // 2
        prov = pd.read_csv(next(gd.glob("*_find_proviruses/*_provirus.tsv")), sep="\t")
        seqs = read_fasta(Path(a.fasta_dir) / f"{g}.fna")
        idx = {c: kmer_index(s) for c, s in seqs.items()}
        gprof = profile(np.concatenate(list(idx.values())))
        for r in rg.itertuples():
            x = genes[(genes.contig == r.contig) & (genes.mid > r.start) & (genes.mid <= r.end)]
            desc = x.annotation_description.fillna("").astype(str)
            n = r.end - r.start
            at = jsd(profile(idx[r.contig][r.start:max(r.start, r.end - 3)]), gprof)
            pool = [(c, len(s)) for c, s in seqs.items() if len(s) > n]
            null = []
            for _ in range(N_WIN):
                c, L = pool[int(rng.integers(len(pool)))]
                s0 = int(rng.integers(0, L - n))
                null.append(jsd(profile(idx[c][s0:s0 + n - 3]), gprof))
            p = prov[(prov.source_seq == r.contig) & (prov.end > r.start) & (prov.start <= r.end)]
            sig = dict(genomad_provirus=len(p) > 0, phage_v21=bool(r.v2_pass), phage_v1=bool(r.v1_pass) if pd.notna(r.v1_pass) else False,
                       conjscan=int(x.annotation_conjscan.notna().sum()), plasmid_hallmark=int(x.plasmid_hallmark.sum()),
                       virus_hallmark=int(x.virus_hallmark.sum()), integrase_family=int(desc.str.contains(MOBILE).sum()),
                       att=str(getattr(r, "boundary", "")) == "att", amr=int(x.annotation_amr.notna().sum()),
                       atypicality_pct=float((np.array(null) < at).mean()), genes=len(x),
                       genes_per_kb=round(len(x) / max(n / 1000, 1e-9), 2))
            sig["mobile_any"] = bool(sig["conjscan"] or sig["plasmid_hallmark"] or sig["virus_hallmark"] or sig["integrase_family"] or sig["att"])
            sig["atypical"] = sig["atypicality_pct"] >= 0.95
            b = (0 if sig["genomad_provirus"] else 1 if sig["phage_v21"] else 2 if sig["mobile_any"] else 3 if sig["atypical"] else 4)
            for k, pat in THEMES.items():
                sig[f"theme_{k}"] = int(desc.str.contains(pat, flags=re.I, regex=True).sum())
            rows.append(dict(region=r.region, kind=r.kind, genome=g, kb=r.kb, bin=b) | sig)
            for d in set(desc[desc != ""]):
                desc_rows.append((r.kind, b, d))
        print(f"[{g}] {len(rg)} regions", flush=True)
    t = pd.DataFrame(rows); t.to_csv(out / "regions.tsv", sep="\t", index=False)
    names = {0: "known to geNomad", 1: "phage-like", 2: "other mobile element", 3: "island, function unknown", 4: "no detectable signal"}
    res = {}
    for k in KINDS:
        x = t[t.kind == k]
        res[k] = {"n": len(x), "bins": {}}
        for b in range(5):
            col = f"is_{b}"; x = x.assign(**{col: x.bin == b})
            res[k]["bins"][names[b]] = dict(share=round(float(x[col].mean()), 3), ci=boot_share(x, col),
                                            bp_share=round(float(x[x[col]].kb.sum() / x.kb.sum()), 3))
    # multi-label shares and themes, candidates against length-matched host stretches
    lab = ["phage_v21", "phage_v1", "mobile_any", "atypical", "conjscan", "plasmid_hallmark", "virus_hallmark", "integrase_family", "att", "amr"]
    multi = {k: {c: round(float((t[t.kind == k][c] > 0).mean()), 3) for c in lab + [f"theme_{h}" for h in THEMES]} for k in KINDS}
    dd = pd.DataFrame(desc_rows, columns=["kind", "bin", "desc"])
    cand4 = dd[(dd.kind == "candidate_default") & (dd.bin == 4)].desc.value_counts()
    neg = dd[dd.kind == "negative_control"].desc.value_counts()
    nc4 = int(((t.kind == "candidate_default") & (t.bin == 4)).sum()); nn = int((t.kind == "negative_control").sum())
    enr = pd.DataFrame({"bin4_candidates": cand4, "random": neg}).fillna(0)
    enr["share_bin4"] = enr.bin4_candidates / max(nc4, 1); enr["share_random"] = enr.random / max(nn, 1)
    enr["ratio"] = (enr.share_bin4 + 1e-3) / (enr.share_random + 1e-3)
    enr = enr[enr.bin4_candidates >= 10].sort_values("ratio", ascending=False)
    enr.to_csv(out / "bin4_descriptions.tsv", sep="\t")
    json.dump(dict(bins=res, multi_label=multi), open(out / "results.json", "w"), indent=1)
    for k in KINDS:
        print(f"\n== {k} (n={res[k]['n']})")
        for b, v in res[k]["bins"].items():
            print(f"   {b:28s} {v['share']:.3f} {v['ci']}  (bp share {v['bp_share']})")
    print("\n== multi-label shares (candidates | length-matched | composition-matched)")
    for c in lab + [f"theme_{h}" for h in THEMES]:
        print(f"   {c:28s} {multi['candidate_default'][c]:.3f} | {multi['negative_control'][c]:.3f} | {multi['composition_control'][c]:.3f}")
    print("\n== bin-4 candidate descriptions most enriched over random stretches (>=10 candidates)")
    print(enr.head(20).round(3).to_string())


if __name__ == "__main__":
    main()
