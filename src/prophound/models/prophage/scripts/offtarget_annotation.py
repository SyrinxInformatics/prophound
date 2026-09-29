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
"""What are the off-target calls? (docs/error-analysis.md)

For every call in --calls, and for length-matched random control intervals on
the same genome (clear of every call and every geNomad element), summarise the
genes geNomad's `annotate` step found inside it:

  - replicon: plasmid or chromosome, from the FASTA header
  - gene density and GC shift (mean gene GC minus the genome's)
  - functional terms in annotation_description: phage structural, lysis,
    integrase, transposase, conjugation, plasmid replication, generic phage
  - geNomad's own flags where the table has them (virus_hallmark,
    plasmid_hallmark, annotation_conjscan, annotation_amr, marker)

Functional terms come first because geNomad's markers partly restate the labels
the model was trained on (docs/scientific-validity.md, item 1b caveat). Each
interval gets one first-match category. Groups compared: on-target calls,
off-target calls by seed support, and controls. A category is evidence about
what a call resembles, not proof of what it is.

Example:
  python scripts/offtarget_annotation.py --calls analysis/offtarget/calls42.csv \\
      --elements scripts/eval_sizestrat_elements.csv --genomad_dir ../genomad_out \\
      --genomes_dir ../raw/genomes --out analysis/offtarget/annotation.json
"""

from __future__ import annotations

import argparse
import json
import warnings
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

TERMS = {
    "structural": r"terminase|portal|capsid|major head|head protein|prohead|tail|baseplate|head-tail|neck",
    "lysis": r"holin|endolysin|lysozyme|spanin|\blysis\b",
    "integrase": r"integrase|site-specific recombinase|tyrosine recombinase|serine recombinase|excisionase",
    "transposase": r"transposase|insertion element|insertion sequence",
    "conjugation": r"relaxase|conjugal|conjugative|type iv secretion|\bt4ss\b|mobilization|\btra[a-z]\b|\btrb[a-z]\b|\bvirb|\bvird4",
    "plasmid_rep": r"plasmid replication|replication initiat|\brepa\b|partition protein|\bpar[ab]\b",
    "phage_generic": r"phage|prophage|viral|virion",
}
FLAG_COLS = ["virus_hallmark", "plasmid_hallmark", "annotation_conjscan", "annotation_amr", "marker"]
CATEGORIES = ["phage_structural", "conjugative", "plasmid_replicon", "integrase_or_transposase",
              "phage_term_only", "no_mge_hallmark", "no_genes"]


def fasta_contigs(path: Path) -> dict:
    """contig -> (length, is_plasmid) by streaming the FASTA once."""
    out, name, desc, n = {}, None, "", 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    out[name] = (n, "plasmid" in desc.lower())
                parts = line[1:].strip().split(None, 1)
                name, desc, n = parts[0], (parts[1] if len(parts) > 1 else ""), 0
            else:
                n += len(line.strip())
    if name is not None:
        out[name] = (n, "plasmid" in desc.lower())
    return out


def flag_on(v) -> bool:
    if isinstance(v, str):
        return v.strip() not in ("", "NA", "nan", "0", "False")
    return bool(v) and not (isinstance(v, float) and np.isnan(v))


def features(genes: pd.DataFrame, contigs: dict, genome_gc: float, contig: str, s: int, e: int) -> dict:
    g = genes[(genes.contig == contig) & (genes.mid >= s) & (genes.mid < e)]
    desc = g.annotation_description.fillna("").str.lower() if "annotation_description" in g else pd.Series([], dtype=str)
    f = {"replicon_plasmid": bool(contigs.get(contig, (0, False))[1]), "n_genes": int(len(g)),
         "genes_per_10kb": len(g) / max(1.0, (e - s) / 1e4),
         "gc_shift": float(g.gc_content.mean() - genome_gc) if len(g) and "gc_content" in g else np.nan}
    for k, rx in TERMS.items():
        f[k] = int(desc.str.contains(rx, regex=True).sum()) if len(desc) else 0
    for c in FLAG_COLS:
        if c in g:
            f[c] = int(g[c].map(flag_on).sum())
    if "marker" in g:
        f["marker_suffix"] = dict(Counter(str(m).rsplit(".", 1)[-1] for m in g.marker if flag_on(m)))
    if f["n_genes"] == 0:
        cat = "no_genes"
    elif f["structural"] >= 2:
        cat = "phage_structural"
    elif f["conjugation"] >= 2 or f.get("annotation_conjscan", 0) >= 1:
        cat = "conjugative"
    elif f["replicon_plasmid"]:
        cat = "plasmid_replicon"
    elif f["integrase"] + f["transposase"] >= 1:
        cat = "integrase_or_transposase"
    elif f["structural"] + f["phage_generic"] >= 1:
        cat = "phage_term_only"
    else:
        cat = "no_mge_hallmark"
    f["category"] = cat
    return f


def controls(rng, contigs: dict, busy: dict, length: int, n: int, tries: int = 200) -> list:
    """n random [s, e) of `length` on contigs long enough, clear of `busy` intervals."""
    pool = [(c, L) for c, (L, _) in contigs.items() if L > length]
    out = []
    if not pool:
        return out
    w = np.array([L - length for _, L in pool], dtype=float)
    for _ in range(n):
        for _ in range(tries):
            c, L = pool[rng.choice(len(pool), p=w / w.sum())]
            s = int(rng.integers(0, L - length))
            if not any(bs < s + length and be > s for bs, be in busy.get(c, [])):
                out.append((c, s, s + length))
                break
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--calls", required=True, help="sample,mode,genome_id,contig_id,start,end,seed_support,on_target,...")
    ap.add_argument("--elements", required=True, help="geNomad elements (controls avoid them)")
    ap.add_argument("--genomad_dir", required=True)
    ap.add_argument("--genomes_dir", required=True)
    ap.add_argument("--n_controls", type=int, default=3, help="controls per call")
    ap.add_argument("--control_busy", choices=("calls_and_elements", "elements_only", "none"),
                    default="calls_and_elements",
                    help="what controls must avoid. The default avoids every call, which conditions the "
                         "background on the model not calling there; 'elements_only' avoids geNomad "
                         "elements only (model-independent); 'none' is fully random "
                         "(docs/research-integrity.md F8)")
    ap.add_argument("--seed", type=int, default=20260928)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    warnings.filterwarnings("ignore", message="Mean of empty slice")   # control groups with no genes

    calls = pd.read_csv(args.calls)
    ev = pd.read_csv(args.elements)
    rng = np.random.default_rng(args.seed)
    rows, missing, columns_seen = [], [], None
    for g, cg in calls.groupby("genome_id"):
        gt = Path(args.genomad_dir) / g / f"{g}_annotate" / f"{g}_genes.tsv"
        fa = Path(args.genomes_dir) / f"{g}.fna"
        if not gt.exists() or not fa.exists():
            missing.append(g)
            continue
        genes = pd.read_csv(gt, sep="\t")
        columns_seen = columns_seen or list(genes.columns)
        genes["contig"] = genes.gene.str.rsplit("_", n=1).str[0]
        genes["mid"] = (genes.start + genes.end) // 2
        contigs = fasta_contigs(fa)
        genome_gc = float(genes.gc_content.mean()) if "gc_content" in genes else np.nan
        busy = {}
        avoid = {"calls_and_elements": [cg[["contig_id", "start", "end"]].values,
                                        ev[ev.accession == g][["contig", "start", "end"]].values],
                 "elements_only": [ev[ev.accession == g][["contig", "start", "end"]].values],
                 "none": []}[args.control_busy]
        for x in avoid:
            for c, s, e in x:
                busy.setdefault(c, []).append((int(s), int(e)))
        for r in cg.itertuples():
            grp = "on_target" if r.on_target else ("off_3seed" if r.seed_support >= 3 else "off_1_2seed")
            base = {"sample": r.sample, "mode": r.mode, "genome_id": g, "group": grp}
            rows.append({**base, "contig_id": r.contig_id, "start": int(r.start), "end": int(r.end),
                         "seed_support": int(r.seed_support),
                         **features(genes, contigs, genome_gc, r.contig_id, int(r.start), int(r.end))})
            for c, s, e in controls(rng, contigs, busy, int(r.end - r.start), args.n_controls):
                rows.append({**base, "group": f"control_{grp}", "contig_id": c, "start": s, "end": e,
                             "seed_support": -1, **features(genes, contigs, genome_gc, c, s, e)})
        print(f"[genome] {g}: {len(cg)} calls", flush=True)

    df = pd.DataFrame(rows)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.drop(columns=["marker_suffix"], errors="ignore").to_csv(out.with_suffix(".csv"), index=False)
    summary = {"genes_columns": columns_seen, "missing_genomes": missing, "control_busy": args.control_busy,
               "groups": {}}
    feat = [k for k in ["replicon_plasmid", *TERMS, *FLAG_COLS] if k in df]
    for (mode, grp), sub in df.groupby(["mode", "group"]):
        s = {"n": int(len(sub)),
             "category": {c: round(float((sub.category == c).mean()), 4) for c in CATEGORIES},
             "frac_with": {k: round(float((sub[k] > 0).mean()), 4) for k in feat},
             "genes_per_10kb_median": float(sub.genes_per_10kb.median()),
             "gc_shift_median": float(sub.gc_shift.median())}
        if "marker_suffix" in sub:
            tot = Counter()
            for d in sub.marker_suffix.dropna():
                tot.update(d)
            s["marker_suffix_total"] = dict(tot.most_common(12))
        summary["groups"][f"{mode}/{grp}"] = s
        print(f"[{mode}/{grp}] n={s['n']} " + " ".join(f"{c}={v:.2f}" for c, v in s["category"].items() if v)
              + f" | gc_shift {s['gc_shift_median']:.4f} genes/10kb {s['genes_per_10kb_median']:.1f}", flush=True)
    out.write_text(json.dumps(summary, indent=2) + "\n")
    print(f"[done] wrote {out} and {out.with_suffix('.csv')}; missing genomes: {len(missing)}", flush=True)


if __name__ == "__main__":
    main()
