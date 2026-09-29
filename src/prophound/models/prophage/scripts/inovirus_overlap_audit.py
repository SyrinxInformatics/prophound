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
"""Where do inoviruses already sit in our corpus, and how are they labelled?

docs/v1-scope.md A.1: v1 retraining adds real small elements (integrated
inoviruses) as positives. Before that, check the corpus for label conflicts:
an inovirus that geNomad did not call is currently a `near_neg`/`far_neg`
interval, and adding the same kind of sequence as a positive would put it in
both classes. The same pass finds inoviruses in test/eval genomes, which must
stay out of training.

Inovirus loci are found two independent ways, from geNomad `annotate` output
that already exists for every genome (no gene calling here):
  - `tax`: genes whose geNomad taxonomy falls in Faserviricetes/Tubulavirales
    (Inoviridae and relatives). Same tool as the labels, so it asks where
    geNomad sees inovirus genes vs where its provirus caller drew regions.
  - `pI`:  hmmsearch of the Roux et al. 2019 morphogenesis-protein (pI) HMMs
    (Inovirus_detector's Final_marker_morph.hmm) against geNomad's proteins.
    Independent of geNomad, so it can see inoviruses geNomad missed. Optional:
    skipped if --pi_hmm is not given.
Tagged genes on one contig within --merge_gap bp are merged into a locus.
A pI hit alone is Roux's step-1 seed, not a confirmed inovirus; loci are
reported by evidence level so the count can be read as a range.
"""

import argparse
import json
import shutil
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

INO_TAX = ("Faserviricetes", "Tubulavirales", "Inoviridae")
LABELS = ("positive", "near_neg", "far_neg")


def contig_of(gene_id: str) -> str:
    return gene_id.rsplit("_", 1)[0]


def read_bed(path: Path) -> dict:
    out = defaultdict(list)
    for line in open(path):
        parts = line.split()
        if len(parts) < 4 or line.startswith("#"):
            continue
        out[parts[0]].append((int(parts[1]), int(parts[2]), parts[3]))
    return out


def bed_key(contig: str, bed: dict):
    if contig in bed:
        return contig
    base = contig.split("|")[-1]
    for k in (base, base.rsplit(".", 1)[0]):
        if k in bed:
            return k
    for k in bed:
        if k.rsplit(".", 1)[0] == base.rsplit(".", 1)[0]:
            return k
    return None


def run_hmmsearch(hmm: Path, faa: Path, tmp: Path, evalue: float, cpus: int) -> set:
    tbl = tmp / (faa.stem + ".tbl")
    subprocess.run(["hmmsearch", "--cpu", str(cpus), "-E", str(evalue), "--noali",
                    "--tblout", str(tbl), "-o", "/dev/null", str(hmm), str(faa)], check=True)
    hits = set()
    for line in open(tbl):
        if not line.startswith("#"):
            hits.add(line.split()[0])
    tbl.unlink()
    return hits


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genomad_dir", required=True, help="genomad_out/ (sibling of data_refactored/)")
    ap.add_argument("--beds_dir", required=True, help="beds_all/ (positive + near_neg + far_neg)")
    ap.add_argument("--splits_dir", required=True)
    ap.add_argument("--eval_list", action="append", default=[], help="genome lists to flag (repeatable)")
    ap.add_argument("--pi_hmm", default=None, help="Final_marker_morph.hmm; omit to run tax-only")
    ap.add_argument("--evalue", type=float, default=1e-5)
    ap.add_argument("--merge_gap", type=int, default=10000)
    ap.add_argument("--cpus", type=int, default=4)
    ap.add_argument("--limit", type=int, default=0, help="first N genomes only (smoke test)")
    ap.add_argument("--out_dir", required=True)
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    tmp = out / "_tmp"
    tmp.mkdir(exist_ok=True)
    if args.pi_hmm and not shutil.which("hmmsearch"):
        raise SystemExit("--pi_hmm given but hmmsearch is not on PATH")

    split_of = {}
    for s in ("train", "val", "test"):
        p = Path(args.splits_dir) / f"{s}.txt"
        if p.exists():
            for l in open(p):
                if l.strip():
                    split_of[l.strip()] = s
    evalset = {l.strip() for f in args.eval_list for l in open(f) if l.strip()}

    gdirs = sorted(d for d in Path(args.genomad_dir).iterdir() if d.is_dir())
    if args.limit:
        gdirs = gdirs[:args.limit]
    print(f"[scan] {len(gdirs)} genomad_out genomes, pI HMM: {bool(args.pi_hmm)}", flush=True)

    loci = []
    n_missing_annot = n_no_bed = n_contig_unmatched = 0
    pi_skipped = []  # genomes whose pI search could not run: counted and named, never silently dropped
    for i, gd in enumerate(gdirs, 1):
        g = gd.name
        genes_tsv = gd / f"{g}_annotate" / f"{g}_genes.tsv"
        faa = gd / f"{g}_annotate" / f"{g}_proteins.faa"
        if not genes_tsv.exists():
            n_missing_annot += 1
            continue
        genes = pd.read_csv(genes_tsv, sep="\t", usecols=["gene", "start", "end", "taxname"])
        tax_hit = genes.taxname.astype(str).str.contains("|".join(INO_TAX))
        pi_ids = set()
        if args.pi_hmm:
            # 817355 died at genome 1,600+ on one empty geNomad proteins.faa
            if not faa.exists() or faa.stat().st_size == 0:
                pi_skipped.append((g, "empty_or_missing_faa"))
            else:
                try:
                    pi_ids = run_hmmsearch(Path(args.pi_hmm), faa, tmp, args.evalue, args.cpus)
                except subprocess.CalledProcessError as e:
                    pi_skipped.append((g, f"hmmsearch_exit_{e.returncode}"))
        genes["tax"] = tax_hit
        genes["pI"] = genes.gene.isin(pi_ids)
        tagged = genes[genes.tax | genes.pI].copy()
        if tagged.empty:
            continue
        tagged["contig"] = tagged.gene.map(contig_of)

        bed_path = Path(args.beds_dir) / f"{g}.bed"
        bed = read_bed(bed_path) if bed_path.exists() else None
        if bed is None:
            n_no_bed += 1

        first_new = len(loci)
        for contig, grp in tagged.sort_values("start").groupby("contig"):
            cur = None
            for r in grp.itertuples():
                if cur and r.start - cur["end"] <= args.merge_gap:
                    cur["end"] = max(cur["end"], r.end)
                    cur["n_tax"] += int(r.tax)
                    cur["n_pI"] += int(r.pI)
                else:
                    if cur:
                        loci.append(cur)
                    cur = {"genome_id": g, "contig": contig, "start": int(r.start), "end": int(r.end),
                           "n_tax": int(r.tax), "n_pI": int(r.pI)}
            loci.append(cur)

        # annotate this genome's new loci with BED overlap
        for L in loci[first_new:]:
            ov = {k: 0 for k in LABELS}
            pos_len = 0
            key = bed_key(L["contig"], bed) if bed else None
            if bed and key is None:
                n_contig_unmatched += 1
            for s, e, lab in (bed.get(key, []) if key else []):
                o = min(e, L["end"]) - max(s, L["start"])
                if o > 0 and lab in ov:
                    ov[lab] += o
                    if lab == "positive":
                        pos_len = max(pos_len, e - s)
            L.update({f"ov_{k}": v for k, v in ov.items()})
            L["label"] = ("positive" if ov["positive"] else "near_neg" if ov["near_neg"]
                          else "far_neg" if ov["far_neg"] else "unlabelled" if key else "no_bed")
            L["positive_element_len"] = pos_len
            L["split"] = split_of.get(g, "none")
            L["in_eval"] = g in evalset
            L["evidence"] = "both" if L["n_tax"] and L["n_pI"] else "tax_only" if L["n_tax"] else "pI_only"
        if i % 200 == 0:
            print(f"[scan] {i}/{len(gdirs)}  loci so far {len(loci)}", flush=True)
    tmp.rmdir()

    df = pd.DataFrame(loci)
    df.to_csv(out / "loci.tsv", sep="\t", index=False)

    def tab(sub):
        return {f"{s}/{l}": int(n) for (s, l), n in sorted(Counter(zip(sub.split, sub.label)).items())}

    res = {"n_genomes_scanned": len(gdirs), "missing_annotate": n_missing_annot, "no_bed": n_no_bed,
           "contigs_unmatched_to_bed": n_contig_unmatched, "pI_hmm_used": bool(args.pi_hmm),
           "pI_skipped_genomes": {"n": len(pi_skipped), "list": [f"{g}:{why}" for g, why in pi_skipped]},
           "merge_gap": args.merge_gap, "n_loci": len(df)}
    if len(df):
        res["genomes_with_loci"] = int(df.genome_id.nunique())
        res["by_evidence"] = {k: int(v) for k, v in df.evidence.value_counts().items()}
        for ev in ["all", "both", "tax_only", "pI_only"]:
            sub = df if ev == "all" else df[df.evidence == ev]
            res[f"split_x_label[{ev}]"] = tab(sub)
        conflict = df[(df.split == "train") & df.label.isin(["near_neg", "far_neg"])]
        res["train_conflicts"] = {"n": len(conflict),
                                  "by_evidence": {k: int(v) for k, v in conflict.evidence.value_counts().items()}}
        res["loci_in_test_or_eval"] = int(((df.split == "test") | df.in_eval).sum())
        pos = df[df.label == "positive"]
        res["positive_element_len_quantiles"] = (
            {str(q): float(pos.positive_element_len.quantile(q)) for q in (0.1, 0.25, 0.5, 0.75, 0.9)}
            if len(pos) else {})
        res["inovirus_positive_elements_under_10kb"] = int((pos.positive_element_len.between(1, 9999)).sum())
    res["interpretation"] = (
        "train_conflicts are inovirus loci inside training genomes that the current labels call negative: "
        "adding inoviruses as positives without relabelling these puts the same kind of sequence in both "
        "classes. Read by evidence level -- 'both' is the confident count, 'pI_only' is Roux's step-1 seed and "
        "includes non-inovirus Zot-like ATPases. loci_in_test_or_eval must be excluded from any training "
        "source and kept homology-disjoint from it.")
    (out / "results.json").write_text(json.dumps(res, indent=2) + "\n")
    print(json.dumps(res, indent=2), flush=True)
    print(f"\n[done] wrote {out / 'results.json'}", flush=True)


if __name__ == "__main__":
    main()
