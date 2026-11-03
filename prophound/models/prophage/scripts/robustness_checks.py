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
"""Robustness checks 2-4 (analysis/robustness/), on the second lockbox and out-of-domain input.

  gc_phylum  covered50, calls per genome and no-signal share by genome GC and by phylum
  plasmids   extrachromosomal against chromosomal contigs: calls, bp called, whole-replicon calls
  ood        calls per Mb on out-of-domain inputs against the bacterial floor

Run from prophound/models/prophage:
  python scripts/robustness_checks.py gc_phylum --fasta_dir ../../../arc_data/lockbox2_fasta --out analysis/robustness
  python scripts/robustness_checks.py plasmids --fasta_dir ../../../arc_data/lockbox2_fasta --out analysis/robustness
  python scripts/robustness_checks.py ood --ood_runs analysis/robustness/ood --fasta_dir ../../../arc_data/lockbox2_fasta --out analysis/robustness
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from evalkit import wilson

L2 = Path("analysis/lockbox2")


def headers_and_lengths(fna: Path) -> list[tuple[str, str, int, int]]:
    out, name, desc, n, gc = [], None, "", 0, 0
    for line in open(fna):
        if line.startswith(">"):
            if name:
                out.append((name, desc, n, gc))
            parts = line[1:].strip().split(None, 1)
            name, desc, n, gc = parts[0], parts[1] if len(parts) > 1 else "", 0, 0
        else:
            s = line.strip().upper(); n += len(s); gc += s.count("G") + s.count("C")
    if name:
        out.append((name, desc, n, gc))
    return out


def default_calls(g: str) -> pd.DataFrame:
    return pd.read_csv(L2 / "prophound" / "default" / g / "regions.detailed.tsv", sep="\t")


def cmd_gc_phylum(a) -> None:
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    tax = pd.read_csv(L2 / "lockbox_taxonomy.tsv", sep="\t").set_index("genome_id")
    el = pd.read_csv(L2 / "results.elements.csv")
    el = el[el["mode"] == "default"]
    xc = pd.read_csv("analysis/extra_calls/regions.tsv", sep="\t")
    xc = xc[xc.kind == "candidate_default"]
    genomes = [l.strip() for l in open(L2 / "lockbox_genome_list.txt") if l.strip()]
    truth_g = set(pd.read_csv(L2 / "lockbox_elements.csv").accession)
    rows = []
    for g in genomes:
        hl = headers_and_lengths(Path(a.fasta_dir) / f"{g}.fna")
        bp = sum(x[2] for x in hl); gc = sum(x[3] for x in hl) / bp
        rows.append(dict(genome=g, gc=gc, phylum=tax.loc[g, "phylum"] if isinstance(tax.loc[g, "phylum"], str) and tax.loc[g, "phylum"] else "unassigned",
                         calls=len(default_calls(g)), element_free=g not in truth_g, mb=bp / 1e6))
    gt = pd.DataFrame(rows)
    gt["gc_bin"] = pd.cut(gt.gc, [0, 0.40, 0.55, 0.65, 1], labels=["<40%", "40-55%", "55-65%", ">65%"])
    big = gt.phylum.value_counts(); keep = set(big[big >= 5].index)
    gt["phylum_bin"] = np.where(gt.phylum.isin(keep), gt.phylum, "other")
    el = el.merge(gt[["genome", "gc_bin", "phylum_bin"]], left_on="accession", right_on="genome")
    xc = xc.merge(gt[["genome", "gc_bin", "phylum_bin"]], on="genome")
    pooled_free = gt[gt.element_free].calls.mean()
    res = {}
    for col in ("gc_bin", "phylum_bin"):
        tab = []
        for k, sub in gt.groupby(col, observed=True):
            e = el[el[col] == k]; x = xc[xc[col] == k]; free = sub[sub.element_free]
            lo, hi = wilson(int(e.covered50.sum()), len(e)) if len(e) else (np.nan, np.nan)
            warn = (len(e) >= 10 and hi < 0.90) or (len(free) >= 5 and free.calls.mean() >= 2 * pooled_free)
            tab.append(dict(stratum=str(k), genomes=len(sub), elements=len(e), covered50=round(e.covered50.mean(), 3) if len(e) else None,
                            wilson=[round(lo, 3), round(hi, 3)] if len(e) else None, calls_per_genome=round(sub.calls.mean(), 1),
                            element_free_genomes=len(free), calls_per_free_genome=round(free.calls.mean(), 1) if len(free) else None,
                            extra_calls=len(x), no_signal_share=round(float((x.bin == 4).mean()), 3) if len(x) else None, warning=bool(warn)))
        res[col] = tab
        print(f"\n== by {col} (pooled element-free calls per genome {pooled_free:.1f})")
        print(pd.DataFrame(tab).to_string(index=False))
    json.dump(dict(pooled_calls_per_free_genome=pooled_free, strata=res), open(out / "gc_phylum.json", "w"), indent=1, default=str)


def cmd_plasmids(a) -> None:
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    genomes = [l.strip() for l in open(L2 / "lockbox_genome_list.txt") if l.strip()]
    truth = pd.read_csv(L2 / "lockbox_elements.csv")
    el = pd.read_csv(L2 / "results.elements.csv"); el = el[el["mode"] == "default"]
    rows = []
    for g in genomes:
        d = default_calls(g)
        for name, desc, n, _ in headers_and_lengths(Path(a.fasta_dir) / f"{g}.fna"):
            low = desc.lower()
            kind = "plasmid" if "plasmid" in low else "chromosome" if ("chromosome" in low or "complete genome" in low) else "unclassified"
            c = d[d.contig_id == name]
            iv = sorted(zip(c.start, c.end)); tot, cur = 0, None
            for s, e in iv:
                if cur and s <= cur[1]:
                    cur = (cur[0], max(cur[1], e))
                else:
                    if cur:
                        tot += cur[1] - cur[0]
                    cur = (s, e)
            tot += (cur[1] - cur[0]) if cur else 0
            te = truth[(truth.accession == g) & (truth.contig == name)]
            rows.append(dict(genome=g, contig=name, kind=kind, kb=n / 1000, calls=len(c), bp_called=tot,
                             frac_called=tot / n if n else 0, max_call_frac=float(((c.end - c.start) / n).max()) if len(c) else 0.0,
                             truth_elements=len(te)))
    t = pd.DataFrame(rows); t.to_csv(out / "replicons.tsv", sep="\t", index=False)
    summ = t.groupby("kind").agg(contigs=("contig", "size"), with_call=("calls", lambda x: round((x > 0).mean(), 3)),
                                 bp_called_share=("bp_called", "sum"), kb=("kb", "sum"),
                                 mostly_called=("frac_called", lambda x: round((x >= 0.5).mean(), 3)),
                                 one_call_covers_most=("max_call_frac", lambda x: round((x >= 0.5).mean(), 3)),
                                 truth_elements=("truth_elements", "sum"), median_kb=("kb", "median"))
    summ["bp_called_share"] = (summ.bp_called_share / (summ.kb * 1000)).round(3)
    pl = t[t.kind == "plasmid"]
    pe = el.merge(t[t.kind == "plasmid"][["genome", "contig"]], left_on=["accession", "contig"], right_on=["genome", "contig"])
    chrom = summ.loc["chromosome", "bp_called_share"] if "chromosome" in summ.index else np.nan
    plas = summ.loc["plasmid", "bp_called_share"] if "plasmid" in summ.index else np.nan
    warn = bool((plas >= 2 * chrom) or ((pl.frac_called >= 0.5).mean() > 0.10 if len(pl) else False))
    print(summ.drop(columns=["kb"]).to_string()); print(f"\nplasmid truth elements: {len(pe)}, covered50 {pe.covered50.mean() if len(pe) else float('nan'):.3f}")
    print(f"plasmid bp called {plas} vs chromosome {chrom} (ratio {plas / chrom:.2f}); warning: {warn}")
    json.dump(dict(summary=summ.reset_index().to_dict(orient="records"), plasmid_truth_elements=len(pe),
                   plasmid_truth_covered50=float(pe.covered50.mean()) if len(pe) else None, warning=warn), open(out / "plasmids.json", "w"), indent=1, default=float)


def cmd_ood(a) -> None:
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    genomes = [l.strip() for l in open(L2 / "lockbox_genome_list.txt") if l.strip()]
    truth_g = set(pd.read_csv(L2 / "lockbox_elements.csv").accession)
    free = [g for g in genomes if g not in truth_g]
    mb = sum(sum(x[2] for x in headers_and_lengths(Path(a.fasta_dir) / f"{g}.fna")) for g in free) / 1e6
    floor = sum(len(default_calls(g)) for g in free) / mb
    labels = {"GCF_000091665.1": "archaea: Methanocaldococcus jannaschii", "GCF_000012285.1": "archaea: Sulfolobus acidocaldarius",
              "GCF_000025685.1": "archaea: Haloferax volcanii", "GCF_000007065.1": "archaea: Methanosarcina mazei",
              "GCF_000009965.1": "archaea: Thermococcus kodakarensis", "GCF_000146045.2": "eukaryote: Saccharomyces cerevisiae",
              "human_chr21_slice": "eukaryote: human chr21 10 Mb", "NC_006998.1": "eukaryotic virus: vaccinia",
              "NC_001806.2": "eukaryotic virus: herpes simplex 1", "NC_014649.1": "eukaryotic virus: mimivirus"}
    rows = []
    for d in sorted(Path(a.ood_runs).iterdir()):
        if not (d / "run_metadata.json").exists():
            continue
        m = json.load(open(d / "run_metadata.json")); bp = m["input"]["total_bp"]
        c = pd.read_csv(d / "regions.detailed.tsv", sep="\t")
        rows.append(dict(input=labels.get(d.name, d.name), mb=round(bp / 1e6, 2), calls=len(c), calls_per_mb=round(len(c) / (bp / 1e6), 2),
                         frac_called=round(((c.end - c.start).sum()) / bp, 3), median_call_kb=round(((c.end - c.start) / 1000).median(), 1) if len(c) else None))
    t = pd.DataFrame(rows)
    t["vs_floor"] = (t.calls_per_mb / floor).round(2)
    t.to_csv(out / "ood.tsv", sep="\t", index=False)
    print(f"bacterial floor (element-free lockbox-2 genomes): {floor:.2f} Default calls per Mb\n"); print(t.to_string(index=False))
    euk = t[t.input.str.startswith("eukaryote")]
    json.dump(dict(floor_calls_per_mb=floor, rows=rows, eukaryote_warning=bool((euk.calls_per_mb >= floor).any())), open(out / "ood.json", "w"), indent=1)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for n in ("gc_phylum", "plasmids", "ood"):
        p = sub.add_parser(n); p.add_argument("--fasta_dir", required=True); p.add_argument("--out", required=True)
        if n == "ood":
            p.add_argument("--ood_runs", required=True)
    a = ap.parse_args()
    {"gc_phylum": cmd_gc_phylum, "plasmids": cmd_plasmids, "ood": cmd_ood}[a.cmd](a)


if __name__ == "__main__":
    main()
