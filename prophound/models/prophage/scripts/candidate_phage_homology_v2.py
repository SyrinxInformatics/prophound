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
"""Phage-homology test v2 (docs/candidate-phage-homology.md "Refined test, v2"). Offline.

Stages, each cached under --out so a re-run resumes:
  regions   v1's regions plus composition-matched negative controls (tetranucleotide
            atypicality, seed 20261028) and geNomad-missed curated positives
  genes     pyrodigal (metagenomic mode, table 11) on every region; proteins >= 60 aa
  search    pyhmmer hmmsearch: stage 1 against the phage-specific structural families,
            stage 2 re-searching stage-1 hits against every VOGDB family; a gene counts
            only if its best family overall (full-sequence bitscore, E <= 1e-5) is one of
            the phage-specific structural families
  score     pass = at least 2 such genes; the preregistered endpoints with genome bootstraps

Run with pyhmmer and pyrodigal available (not in prophound's own environment), from
prophound/models/prophage:
  pixi exec -c conda-forge -c bioconda -s python=3.12 -s pyhmmer -s pyrodigal -s pandas -s numpy -- \\
      python scripts/candidate_phage_homology_v2.py --out analysis/candidate_homology/v2 ...
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

E_MAX = 1e-5
MIN_AA = 60
MIN_GENES = 2
PAD = 5000
SEED = 20261028
N_DRAWS = 50
TOL = 0.25
_LUT = np.full(256, 255, dtype=np.uint8)
for _i, _b in enumerate(b"ACGT"):
    _LUT[_b] = _i; _LUT[ord(chr(_b).lower())] = _i


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


def canon_map() -> np.ndarray:
    comp = {0: 3, 1: 2, 2: 1, 3: 0}
    m = np.zeros(256, dtype=np.int64); seen = {}
    for k in range(256):
        d = [(k >> (2 * (3 - i))) & 3 for i in range(4)]
        rc = sum(comp[x] << (2 * i) for i, x in enumerate(d))
        c = min(k, rc)
        m[k] = seen.setdefault(c, len(seen))
    return m


CANON = canon_map()


def kmer_index(seq: str) -> np.ndarray:
    a = _LUT[np.frombuffer(seq.encode(), dtype=np.uint8)].astype(np.int64)
    if len(a) < 4:
        return np.full(0, -1)
    idx = (a[:-3] << 6) | (a[1:-2] << 4) | (a[2:-1] << 2) | a[3:]
    bad = (a[:-3] == 255) | (a[1:-2] == 255) | (a[2:-1] == 255) | (a[3:] == 255)
    idx = CANON[np.where(bad, 0, idx)]
    idx[bad] = -1
    return idx


def profile(idx: np.ndarray) -> np.ndarray:
    v = np.bincount(idx[idx >= 0], minlength=136).astype(float) + 0.5
    return v / v.sum()


def jsd(p: np.ndarray, q: np.ndarray) -> float:
    m = (p + q) / 2
    return float(0.5 * (p * np.log2(p / m)).sum() + 0.5 * (q * np.log2(q / m)).sum())


def stage_regions(a, out: Path) -> pd.DataFrame:
    f = out / "regions.tsv"
    if f.exists():
        return pd.read_csv(f, sep="\t")
    v1 = pd.read_csv(a.v1_regions, sep="\t")
    v1["source"] = "lockbox2"
    truth = pd.read_csv(a.elements)
    rng = np.random.default_rng(SEED)
    rows = []
    for g, cand in v1[v1.kind == "candidate_default"].groupby("genome"):
        seqs = read_fasta(Path(a.fasta_dir) / f"{g}.fna")
        idx = {c: kmer_index(s) for c, s in seqs.items()}
        gprof = profile(np.concatenate(list(idx.values())))
        busy = {c: [(int(r.start), int(r.end)) for r in truth[(truth.accession == g) & (truth.contig == c)].itertuples()] for c in seqs}
        for m in ("default", "exploratory", "conservative"):
            d = pd.read_csv(Path(a.runs) / m / g / "regions.detailed.tsv", sep="\t")
            for r in d.itertuples():
                busy.setdefault(r.contig_id, []).append((int(r.start), int(r.end)))
        for r in cand.itertuples():
            n = int(r.end - r.start)
            target = jsd(profile(idx[r.contig][r.start:max(r.start, r.end - 3)]), gprof)
            best = None
            contigs = [c for c in seqs if len(seqs[c]) >= n + 2 * PAD]
            for _ in range(N_DRAWS):
                if not contigs:
                    break
                c = contigs[int(rng.integers(len(contigs)))]
                s = int(rng.integers(0, len(seqs[c]) - n)); e = s + n
                if not all(e + PAD <= x or s - PAD >= y for x, y in busy.get(c, [])):
                    continue
                at = jsd(profile(idx[c][s:max(s, e - 3)]), gprof)
                if best is None or abs(at - target) < abs(best[3] - target):
                    best = (c, s, e, at)
            row = dict(region=f"comp_{r.region}", kind="composition_control", genome=g, paired=r.region,
                       candidate_atypicality=round(target, 5), source="lockbox2")
            if best and abs(best[3] - target) <= TOL * target:
                c, s, e, at = best
                busy.setdefault(c, []).append((s, e))
                rows.append(row | dict(contig=c, start=s, end=e, kb=round(n / 1000, 1), atypicality=round(at, 5), matched=True))
            else:
                rows.append(row | dict(contig="", start=-1, end=-1, kb=round(n / 1000, 1),
                                       atypicality=round(best[3], 5) if best else np.nan, matched=False))
    comp = pd.DataFrame(rows)
    cur = pd.read_csv(a.curated, sep="\t")
    cur = cur.assign(region=[f"gmiss_{i}" for i in range(len(cur))], kind="positive_genomad_missed",
                     kb=(cur.end - cur.start) / 1000)[["region", "kind", "genome", "contig", "start", "end", "kb", "source", "name", "cls", "genomad_cov"]]
    t = pd.concat([v1, comp, cur], ignore_index=True)
    t.to_csv(f, sep="\t", index=False)
    return t


def fasta_for(row, a) -> Path:
    if row.source == "LAMBDA":
        return Path(a.lambda_fasta) / f"{row.genome}.fna"
    if row.source == "literature":
        return Path(a.lit_fasta) / f"{row.genome}.fna"
    return Path(a.fasta_dir) / f"{row.genome}.fna"


def stage_genes(a, out: Path, reg: pd.DataFrame) -> Path:
    f = out / "proteins.faa"
    if f.exists():
        return f
    import pyrodigal
    gf = pyrodigal.GeneFinder(meta=True)
    cache, n = {}, 0
    with open(f, "w") as fh, open(out / "genes_per_region.tsv", "w") as gt:
        gt.write("region\tgenes\n")
        for r in reg[reg.start >= 0].itertuples():
            p = fasta_for(r, a)
            if p not in cache:
                cache.clear(); cache[p] = read_fasta(p)
            seq = cache[p][r.contig][int(r.start):int(r.end)]
            k = 0
            for i, gene in enumerate(gf.find_genes(seq.encode())):
                prot = gene.translate(translation_table=11).rstrip("*")
                if len(prot) >= MIN_AA:
                    fh.write(f">{r.region}|{i}\n{prot}\n"); k += 1
            gt.write(f"{r.region}\t{k}\n"); n += k
    print(f"[genes] {n} proteins", flush=True)
    return f


def load_hmms(paths: list[Path]):
    import pyhmmer
    hmms = []
    for p in paths:
        with pyhmmer.plan7.HMMFile(p) as h:
            hmms.extend(h)
    return hmms


def stage_search(a, out: Path, faa: Path) -> pd.DataFrame:
    f = out / f"best_family.{a.family_column}.tsv"
    if f.exists():
        return pd.read_csv(f, sep="\t")
    import pyhmmer
    cls = pd.read_csv(a.vog_classes, sep="\t")
    spec = set(cls[cls[a.family_column]].vog)
    hmm_dir = Path(a.hmm_dir)
    alphabet = pyhmmer.easel.Alphabet.amino()
    with pyhmmer.easel.SequenceFile(faa, digital=True, alphabet=alphabet) as sf:
        seqs = sf.read_block()
    s1 = load_hmms([hmm_dir / f"{v}.hmm" for v in sorted(spec)])
    hit = set()
    for top in pyhmmer.hmmsearch(s1, seqs, cpus=a.cpus, E=E_MAX):
        for h in top:
            if h.evalue <= E_MAX:
                hit.add(h.name if isinstance(h.name, str) else h.name.decode())
    print(f"[search] stage 1: {len(hit)} proteins hit a phage-specific structural family", flush=True)
    sub = pyhmmer.easel.DigitalSequenceBlock(alphabet, [s for s in seqs if (s.name if isinstance(s.name, str) else s.name.decode()) in hit])
    best = {}
    allp = sorted(hmm_dir.glob("VOG*.hmm"))
    for i in range(0, len(allp), 2000):
        block = load_hmms(allp[i:i + 2000])
        for top in pyhmmer.hmmsearch(block, sub, cpus=a.cpus, E=E_MAX):
            q = top.query.name if isinstance(top.query.name, str) else top.query.name.decode()
            for h in top:
                if h.evalue > E_MAX:
                    continue
                n = h.name if isinstance(h.name, str) else h.name.decode()
                if n not in best or h.score > best[n][1]:
                    best[n] = (q, h.score, h.evalue)
        print(f"[search] stage 2: {min(i + 2000, len(allp))}/{len(allp)} families", flush=True)
    rows = [dict(protein=p, region=p.rsplit("|", 1)[0], best_vog=v, bitscore=s, evalue=e, phage_specific=v in spec)
            for p, (v, s, e) in best.items()]
    t = pd.DataFrame(rows, columns=["protein", "region", "best_vog", "bitscore", "evalue", "phage_specific"])
    t.to_csv(f, sep="\t", index=False)
    return t


def boot_diff(a_df: pd.DataFrame, b_df: pd.DataFrame, col: str, B: int = 10000) -> tuple[float, list]:
    g = pd.DataFrame({"ap": a_df.groupby("genome")[col].sum(), "an": a_df.groupby("genome").size(),
                      "bp": b_df.groupby("genome")[col].sum(), "bn": b_df.groupby("genome").size()}).fillna(0)
    rng = np.random.default_rng(20261027); out = []
    for _ in range(B):
        s = g.loc[rng.choice(g.index.to_numpy(), len(g))]
        if s.an.sum() and s.bn.sum():
            out.append(s.ap.sum() / s.an.sum() - s.bp.sum() / s.bn.sum())
    return float(a_df[col].mean() - b_df[col].mean()), [float(x) for x in np.percentile(out, [2.5, 97.5])]


def stage_score(a, out: Path, reg: pd.DataFrame, best: pd.DataFrame) -> None:
    n_spec = best[best.phage_specific].groupby("region").size()
    reg = reg.copy()
    reg["v2_genes"] = reg.region.map(n_spec).fillna(0).astype(int)
    reg["v2_pass"] = reg.v2_genes >= MIN_GENES
    v1 = pd.read_csv(a.v1_scored, sep="\t")[["region", "protein_pass"]].rename(columns={"protein_pass": "v1_pass"})
    reg = reg.merge(v1, on="region", how="left")
    reg.to_csv(out / "scored.tsv", sep="\t", index=False)
    ok = reg[(reg.kind != "composition_control") | (reg.matched == True)]
    summ = ok.groupby("kind").agg(n=("region", "size"), v2_pass=("v2_pass", "sum"), v2_rate=("v2_pass", "mean"),
                                  v1_rate=("v1_pass", "mean"))
    cand = ok[ok.kind == "candidate_default"]; comp = ok[ok.kind == "composition_control"]; neg = ok[ok.kind == "negative_control"]
    d_comp, ci_comp = boot_diff(cand, comp, "v2_pass")
    d_len, ci_len = boot_diff(cand, neg, "v2_pass")
    unmatched = int((reg.kind == "composition_control").sum() - len(comp))
    named = ["cand_d_1787", "cand_d_1911", "cand_d_6761", "cand_d_6921", "cand_d_723", "cand_d_709", "cand_d_716", "cand_d_5281",
             "cand_d_4466", "cand_d_1830", "cand_d_1250", "cand_d_1557", "cand_d_5874", "cand_d_2922"]
    res = dict(primary=dict(difference=d_comp, genome_boot95=ci_comp, candidate_rate=float(cand.v2_pass.mean()),
                            composition_control_rate=float(comp.v2_pass.mean()), n_candidates=len(cand),
                            n_composition_controls=len(comp), unmatched_composition_controls=unmatched),
               secondary=dict(vs_length_matched=dict(difference=d_len, genome_boot95=ci_len),
                              by_set=summ.reset_index().to_dict(orient="records"),
                              named=reg[reg.region.isin(named)][["region", "v2_genes", "v2_pass", "v1_pass"]].to_dict(orient="records")))
    json.dump(res, open(out / "results.json", "w"), indent=1, default=float)
    print(summ.to_string()); print(json.dumps(res["primary"], indent=1)); print(json.dumps(res["secondary"]["vs_length_matched"]))
    print(pd.DataFrame(res["secondary"]["named"]).to_string(index=False))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k in ("out", "v1_regions", "v1_scored", "elements", "runs", "fasta_dir", "curated", "lambda_fasta", "lit_fasta", "vog_classes", "hmm_dir"):
        ap.add_argument(f"--{k}", required=True)
    ap.add_argument("--cpus", type=int, default=6)
    ap.add_argument("--family_column", default="phage_specific_structural",
                    help="vog_classes.tsv column naming the counted families (v2: phage_specific_structural; v2.1: v21)")
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    reg = stage_regions(a, out)
    print(reg.groupby("kind").size().to_string(), flush=True)
    faa = stage_genes(a, out, reg)
    best = stage_search(a, out, faa)
    stage_score(a, out, reg, best)


if __name__ == "__main__":
    main()
