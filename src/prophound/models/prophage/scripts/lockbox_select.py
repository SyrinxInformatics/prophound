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
"""Select the locked test set (docs/lockbox.md): genus-disjoint new genomes.

  select    From NCBI's assembly_summary_refseq.txt, take bacterial RefSeq
            `Complete Genome` assemblies (latest version) whose genus is absent
            from the training corpus -- which is one genome per genus, so genus
            disjointness is the strongest cheap hold-out available. One assembly
            per genus (reference > representative > other, then earliest
            accession), then a seeded random sample of genera. Writes the list,
            a selection table with provenance, and the list's SHA-256.
  manifest  After geNomad has run, write the element manifest (accession,
            contig, start, end, len_bp, stratum) in eval_sizestrat_elements.csv's
            format, strata as in build_eval_sample.py.

Nothing here reads a model output. Example:
  python scripts/lockbox_select.py select --summary assembly_summary_refseq.txt \\
      --corpus_lists splits/train.txt,splits/val.txt,splits/test.txt \\
      --n 100 --seed 20260929 --out_dir analysis/lockbox
  python scripts/lockbox_select.py manifest --genomad_dir ../genomad_out_lockbox \\
      --genome_list analysis/lockbox/lockbox_genome_list.txt --out analysis/lockbox/lockbox_elements.csv
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

STRATA = [("5-8kb", 5.0, 8.0), ("8-10kb", 8.0, 10.0), ("10-15kb", 10.0, 15.0),
          ("15-30kb", 15.0, 30.0), ("30-60kb", 30.0, 60.0), ("60kb+", 60.0, float("inf"))]
CATEGORY_RANK = {"reference genome": 0, "representative genome": 1}
BAD_GENUS = {"uncultured", "bacterium", "unclassified", "candidate", "unidentified"}


def read_summary(path: str) -> pd.DataFrame:
    with open(path) as fh:
        header = None
        for line in fh:
            if line.startswith("#") and "assembly_accession" in line:
                header = line.lstrip("#").strip().split("\t")
                break
    if header is None:
        raise SystemExit(f"{path}: no '#assembly_accession' header line")
    return pd.read_csv(path, sep="\t", comment="#", header=None, names=header, dtype=str, quoting=3)


def genus_of(organism: str) -> str:
    toks = str(organism).replace("[", "").replace("]", "").split()
    if toks and toks[0] == "Candidatus" and len(toks) > 1:
        return f"Candidatus {toks[1]}"
    return toks[0] if toks else ""


def cmd_select(a) -> None:
    s = read_summary(a.summary)
    need = {"assembly_accession", "organism_name", "assembly_level", "version_status", "ftp_path",
            "refseq_category", "seq_rel_date", "asm_name"}
    missing = need - set(s.columns)
    if missing:
        raise SystemExit(f"summary lacks columns {sorted(missing)}")
    group = "group" if "group" in s.columns else None
    s["genus"] = s.organism_name.map(genus_of)

    corpus_ids = set()
    for p in a.corpus_lists.split(","):
        corpus_ids |= {x.strip() for x in Path(p).read_text().split() if x.strip()}
    corpus_acc = {"_".join(g.split("_")[:2]) for g in corpus_ids}
    hit = s[s.assembly_accession.isin(corpus_acc)]
    corpus_genera = set(hit.genus)
    print(f"[corpus] {len(corpus_acc)} accessions, {len(hit)} found in the summary, {len(corpus_genera)} genera",
          flush=True)
    if len(hit) < 0.9 * len(corpus_acc):
        # a superseded accession drops out of the current summary; too many means a wrong file
        raise SystemExit("fewer than 90% of corpus accessions are in the summary -- wrong or truncated file?")

    pool = s[(s.assembly_level == "Complete Genome") & (s.version_status == "latest")
             & (s.ftp_path.str.startswith("https://") | s.ftp_path.str.startswith("ftp://"))]
    if group:
        pool = pool[pool[group] == "bacteria"]
    pool = pool[~pool.genus.isin(corpus_genera) & (pool.genus != "")
                & ~pool.genus.str.lower().isin(BAD_GENUS) & ~pool.assembly_accession.isin(corpus_acc)]
    pool = pool.assign(rank=pool.refseq_category.map(CATEGORY_RANK).fillna(2))
    per_genus = pool.sort_values(["rank", "assembly_accession"]).groupby("genus", sort=True).head(1)
    print(f"[pool] {len(pool)} complete assemblies in {len(per_genus)} genera absent from the corpus", flush=True)
    n = min(a.n, len(per_genus))   # every eligible genus if fewer than requested; recorded below
    if n < a.n:
        print(f"[warn] only {len(per_genus)} eligible genera for n={a.n}: taking all of them", flush=True)
    rng = np.random.default_rng(a.seed)
    pick = per_genus.iloc[np.sort(rng.choice(len(per_genus), size=n, replace=False))].copy()
    pick["ftp_https"] = pick.ftp_path.str.replace("^ftp://", "https://", regex=True)
    base = pick.ftp_https.str.rstrip("/").str.split("/").str[-1]
    pick["genome_id"] = base + "_genomic"
    pick["fna_url"] = pick.ftp_https.str.rstrip("/") + "/" + base + "_genomic.fna.gz"

    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    lst = out / "lockbox_genome_list.txt"
    lst.write_text("".join(f"{g}\n" for g in pick.genome_id))
    cols = ["genome_id", "assembly_accession", "organism_name", "genus", "refseq_category", "seq_rel_date", "fna_url"]
    pick[cols].to_csv(out / "lockbox_selection.tsv", sep="\t", index=False)
    meta = {"summary_file": a.summary, "summary_sha256": sha256(a.summary), "seed": a.seed,
            "n_requested": a.n, "n": int(n),
            "corpus_genera": len(corpus_genera), "eligible_genera": int(len(per_genus)),
            "list_sha256": sha256(str(lst)),
            "rule": "bacteria, RefSeq Complete Genome, latest, genus absent from corpus; one per genus "
                    "(reference > representative > other, then accession); seeded sample of genera"}
    (out / "lockbox_selection.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(f"[done] {n} genomes -> {lst} (sha256 {meta['list_sha256'][:16]}...)", flush=True)


def cmd_manifest(a) -> None:
    rows, missing = [], []
    for g in [x.strip() for x in Path(a.genome_list).read_text().split() if x.strip()]:
        tsv = Path(a.genomad_dir) / g / f"{g}_find_proviruses" / f"{g}_provirus.tsv"
        if not tsv.exists():
            missing.append(g)
            continue
        if tsv.stat().st_size < 10:
            continue
        for r in pd.read_csv(tsv, sep="\t").itertuples():
            # geNomad coordinates are 1-based inclusive; the eval manifests store them as given
            rows.append({"accession": g, "contig": r.source_seq, "start": int(r.start), "end": int(r.end),
                         "len_bp": int(r.end) - int(r.start), "n_genes": int(getattr(r, "n_genes", 0))})
    if missing:
        raise SystemExit(f"{len(missing)} genomes have no geNomad output, e.g. {missing[:3]}")
    e = pd.DataFrame(rows, columns=["accession", "contig", "start", "end", "len_bp", "n_genes"])
    e["kb"] = kb = e.len_bp / 1000.0
    e["stratum"] = "<5kb"
    for lab, lo, hi in STRATA:
        e.loc[(kb >= lo) & (kb < hi), "stratum"] = lab
    e.to_csv(a.out, index=False)
    print(f"[manifest] {len(e)} elements in {e.accession.nunique()} genomes; by stratum "
          f"{e.stratum.value_counts().to_dict()}", flush=True)


def sha256(path: str) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("select")
    s.add_argument("--summary", required=True)
    s.add_argument("--corpus_lists", required=True, help="comma list of genome-id lists (all splits)")
    s.add_argument("--n", type=int, default=100)
    s.add_argument("--seed", type=int, default=20260929)
    s.add_argument("--out_dir", required=True)
    m = sub.add_parser("manifest")
    m.add_argument("--genomad_dir", required=True)
    m.add_argument("--genome_list", required=True)
    m.add_argument("--out", required=True)
    a = ap.parse_args()
    cmd_select(a) if a.cmd == "select" else cmd_manifest(a)


if __name__ == "__main__":
    main()
