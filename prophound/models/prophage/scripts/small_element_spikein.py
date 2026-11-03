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
"""Sub-10 kb spike-ins with real small prophage-class elements (release checklist §3).

geNomad truth cannot certify recall below 10 kb (docs/release-scope-0.1.0.md), so real
small elements are implanted at known coordinates into genomes of the
42-genome sample and scored through the full pipeline. docs/release-scope-0.1.0.md A.1
constraints: prophage class (complete inovirus genomes, not IS elements), real
sequence, homology-disjoint from training.

  select    complete NCBI Inoviridae genomes within a size range, near-duplicates
            collapsed by canonical 21-mer containment (greedy, longest first);
            flags which encode an integrase/recombinase.         [runs anywhere]
  disjoint  containment of each candidate's 21-mers in every training-split
            positive element and every inovirus locus found in a training
            genome; candidates above --max_containment are dropped.  [ARC]
  build     inserts each kept element (random orientation) into --replicates
            different host genomes at quiet sites (every window within
            --quiet_flank of the site below --quiet_prob in the host's clean
            ensemble track, and --element_gap from any geNomad element), sites
            --spacing apart; writes spiked FASTAs and truth. Also writes the
            in-situ truth: confident inovirus loci in val/test genomes, which
            are scored where they sit (detection only: their extent is a
            hallmark-gene span).                                   [ARC]
  eval      rebuilds calls from the score tracks per variant (mode x seed vote
            x refinement) and reports detection (>=1 bp), covered50 and best
            IoU per size bin; in-situ by geNomad label.

Insertion context is artificial (no att site, host composition around it), so
this measures what the pipeline can resolve at a size, not field recall.

Example:
  python scripts/small_element_spikein.py select --genomes_fna ino/genomic.fna \\
      --proteins_faa ino/protein.faa --out_dir scripts/spikein_small
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from evalkit import canon_kmers

SIZE_BINS = [(0, 3000, "<3kb"), (3000, 5000, "3-5kb"), (5000, 8000, "5-8kb"), (8000, 10001, "8-10kb")]


def read_fasta(path) -> list:
    out, name, desc, buf = [], None, "", []
    for line in open(path):
        if line.startswith(">"):
            if name is not None:
                out.append((name, desc, "".join(buf).upper()))
            name, desc, buf = line[1:].split()[0], line[1:].strip(), []
        else:
            buf.append(line.strip())
    if name is not None:
        out.append((name, desc, "".join(buf).upper()))
    return out


def write_fasta(path, recs) -> None:
    with open(path, "w") as fh:
        for desc, seq in recs:
            fh.write(f">{desc}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")


def revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def size_bin(n: int) -> str:
    return next(b for lo, hi, b in SIZE_BINS if lo <= n < hi)


# ---------------------------------------------------------------- select
def source_of(desc: str) -> str:
    """Isolate phage vs metagenome-assembled, from the record's own description."""
    return "metagenome" if re.search(r"\b(MAG|isolate ct|Inoviridae sp\.|Inovirus sp\.|uncultured"
                                     r"|virus \d+ isolate)", desc, re.I) else "isolate"


def cmd_select(a) -> None:
    recs, dropped = [], {"unverified": 0, "construct": 0, "short_metagenome": 0}
    for n, d, s in read_fasta(a.genomes_fna):
        if not a.min_len <= len(s) <= a.max_len:
            continue
        if "UNVERIFIED" in d:            # GenBank's own flag for questionable sequence
            dropped["unverified"] += 1
        elif re.search(r"plasmid|transgenic|vector|construct|synthetic", d, re.I):
            dropped["construct"] += 1    # cloning artefacts, not integrated elements
        elif len(s) < a.min_metagenome_len and source_of(d) == "metagenome":
            dropped["short_metagenome"] += 1   # below any known complete inovirus; likely partial
        else:
            recs.append((n, d, s))
    print(f"[select] dropped {dropped}")
    integ = set()
    for line in open(a.proteins_faa):
        if line.startswith(">") and re.search(r"integrase|recombinase|resolvase", line, re.I):
            m = re.search(r"\[organism=([^\]]+)\]", line)
            if m:
                integ.add(m.group(1).strip())
    recs.sort(key=lambda r: -len(r[2]))
    kept, kmers, members = [], [], []
    for n, d, s in recs:
        km = canon_kmers(s)
        if km.size == 0:
            continue
        best, bi = 0.0, -1
        for i, kk in enumerate(kmers):
            c = np.intersect1d(km, kk, assume_unique=True).size / km.size
            if c > best:
                best, bi = c, i
        if best >= a.dedupe:
            members[bi].append(n)
            continue
        kept.append((n, d, s))
        kmers.append(km)
        members.append([n])
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for (n, d, s), mem in zip(kept, members):
        org = re.sub(r",? (complete|partial) (genome|sequence).*$", "", d.split(" ", 1)[1] if " " in d else d)
        rows.append({"element": n, "name": org, "len_bp": len(s), "gc": round((s.count("G") + s.count("C")) / len(s), 4),
                     "has_integrase": org in integ, "source": source_of(d), "cluster_size": len(mem), "cluster_members": ";".join(mem),
                     "size_bin": size_bin(len(s))})
    df = pd.DataFrame(rows)
    df.to_csv(out / "candidates.tsv", sep="\t", index=False)
    write_fasta(out / "candidates.fna", [(n, s) for n, _, s in kept])
    print(f"[select] {len(recs)} genomes in {a.min_len}-{a.max_len} bp -> {len(df)} representatives "
          f"(containment < {a.dedupe}); by size {df.size_bin.value_counts().to_dict()}; "
          f"with integrase/recombinase {int(df.has_integrase.sum())}")


# ---------------------------------------------------------------- disjoint
def cmd_disjoint(a) -> None:
    from pyfaidx import Fasta
    out = Path(a.candidates_dir)
    df = pd.read_csv(out / "candidates.tsv", sep="\t")
    seqs = {n: s for n, _, s in read_fasta(out / "candidates.fna")}
    codes = [canon_kmers(seqs[e]) for e in df.element]
    owner = np.concatenate([np.full(c.size, i) for i, c in enumerate(codes)])
    allc = np.concatenate(codes)
    order = np.argsort(allc)
    allc, owner = allc[order], owner[order]
    train = set(open(a.train_split).read().split())
    loci = pd.read_csv(a.inovirus_loci, sep="\t")
    loci = loci[(loci.split == "train") & (loci.evidence == "both")]
    seen = {"pos": np.zeros(allc.size, bool), "ino": np.zeros(allc.size, bool)}
    n_el = 0
    for bed in sorted(Path(a.beds_dir).glob("*")):
        g = re.sub(r"\.(bed|tsv|txt)$", "", bed.name)
        if g not in train:
            continue
        fa_path = Path(a.sequences_dir) / f"{g}.fna"
        if not fa_path.exists():
            continue
        ivs = [("pos", p[0], int(p[1]), int(p[2])) for p in
               (l.split() for l in open(bed) if l.strip() and not l.startswith("#")) if len(p) >= 4 and p[3] == "positive"]
        ivs += [("ino", r.contig, int(r.start), int(r.end)) for r in loci[loci.genome_id == g].itertuples()]
        if not ivs:
            continue
        fa = Fasta(str(fa_path), as_raw=True, sequence_always_upper=True)
        try:
            for kind, c, s, e in ivs:
                try:
                    km = canon_kmers(fa[c][max(0, s):e])
                except KeyError:
                    continue
                seen[kind] |= np.isin(allc, km, assume_unique=False)
                n_el += 1
        finally:
            fa.close()
    for kind, col in (("pos", "cont_train_positive"), ("ino", "cont_train_inovirus")):
        hit = np.bincount(owner[seen[kind]], minlength=len(df))
        df[col] = np.round(hit / np.array([c.size for c in codes]), 4)
    df["kept"] = (df.cont_train_positive <= a.max_containment) & (df.cont_train_inovirus <= a.max_containment)
    df.to_csv(out / "candidates_disjoint.tsv", sep="\t", index=False)
    print(f"[disjoint] {n_el} training intervals scanned; kept {int(df.kept.sum())} of {len(df)} "
          f"(max containment {a.max_containment}); by size {df[df.kept].size_bin.value_counts().to_dict()}", flush=True)


# ---------------------------------------------------------------- build
def quiet_sites(w: pd.DataFrame, contig_len: dict, elements: pd.DataFrame, a, rng) -> list:
    """Implant sites. site_rule 'model' (821410) also excludes bins near any window the
    host's clean track scores >= quiet_prob, which conditions the experiment on the
    model; 'genomad' uses geNomad elements and spacing only (docs/research-integrity.md F7)."""
    sites = []
    model_rule = getattr(a, "site_rule", "model") == "model"
    for c, L in contig_len.items():
        if 4 * a.quiet_flank > L:
            continue
        x = w[w.contig_id == c]
        hot = np.zeros(L // 1000 + 1, bool)       # 1 kb bins touched by a window above quiet_prob
        if model_rule:
            for s, e in zip(x.start[x.prob >= a.quiet_prob], x.end[x.prob >= a.quiet_prob]):
                hot[s // 1000:e // 1000 + 1] = True
        for r in elements[elements.contig == c].itertuples():
            hot[max(0, (r.start - a.element_gap) // 1000):(r.end + a.element_gap) // 1000 + 1] = True
        scanned = np.zeros_like(hot)               # bins with no window at all are not quiet either
        for s, e in zip(x.start, x.end):
            scanned[s // 1000:e // 1000 + 1] = True
        f = a.quiet_flank // 1000
        cand = [b * 1000 for b in range(f, len(hot) - f - 1)
                if scanned[b - f:b + f + 1].all() and not hot[b - f:b + f + 1].any()]
        rng.shuffle(cand)
        for p in cand:
            if all(c != c2 or abs(p - p2) >= a.spacing for c2, p2 in sites):
                sites.append((c, p))
    return sites


def cmd_build(a) -> None:
    rng = np.random.default_rng(a.seed)
    cd = Path(a.candidates_dir)
    cand = pd.read_csv(cd / "candidates_disjoint.tsv", sep="\t")
    cand = cand[cand.kept].reset_index(drop=True)
    seqs = {n: s for n, _, s in read_fasta(cd / "candidates.fna")}
    hosts = [h for h in open(a.hosts).read().split() if h]
    ev = pd.read_csv(a.host_elements)
    host_sites, host_seqs = {}, {}
    for h in hosts:
        tw = Path(a.host_tracks) / h / "windows.tsv"
        fa = Path(a.genomes_dir) / f"{h}.fna"
        if not tw.exists() or not fa.exists():
            continue
        recs = read_fasta(fa)
        host_seqs[h] = recs
        w = pd.read_csv(tw, sep="\t", usecols=["contig_id", "start", "end", "prob"])
        host_sites[h] = quiet_sites(w, {n: len(s) for n, _, s in recs}, ev[ev.accession == h], a, rng)
    cap = sum(len(v) for v in host_sites.values())
    print(f"[build] site rule {a.site_rule}; {len(host_sites)} hosts, {cap} quiet sites for up to {len(cand) * a.replicates} implants "
          f"({len(cand)} elements x {a.replicates})", flush=True)
    if cap < len(cand):
        raise SystemExit(f"only {cap} quiet sites for {len(cand)} elements: add hosts or lower --spacing")
    free = {h: list(v) for h, v in host_sites.items()}
    plan = []
    # one pass per replicate, so every element gets a first context before any gets a second
    for rep in range(a.replicates):
        for r in cand.itertuples():
            used = {p["host"] for p in plan if p["element"] == r.element}
            pool = [h for h in free if free[h] and h not in used]
            if not pool:
                continue
            h = max(pool, key=lambda x: len(free[x]))            # drain the roomiest host first
            c, p = free[h].pop()
            plan.append({"host": h, "contig": c, "pos": p, "element": r.element, "replicate": rep,
                         "orient": "+" if rng.random() < 0.5 else "-"})
    plan = pd.DataFrame(plan)
    got = plan.groupby("element").size()
    print(f"[build] contexts per element: {got.value_counts().sort_index().to_dict()} "
          f"(elements with none: {len(cand) - got.size})", flush=True)
    out = Path(a.out_dir)
    (out / "fasta").mkdir(parents=True, exist_ok=True)
    truth = []
    for h, ph in plan.groupby("host"):
        gid = f"spk_{h}"
        new = []
        for n, d, s in host_seqs[h]:
            pc = ph[ph.contig == n].sort_values("pos")
            parts, prev, shift = [], 0, 0
            for r in pc.itertuples():
                el = seqs[r.element] if r.orient == "+" else revcomp(seqs[r.element])
                parts += [s[prev:r.pos], el]
                start = r.pos + shift + 1                       # 1-based inclusive in the spiked contig
                truth.append({"genome_id": gid, "contig": n, "start": start, "end": start + len(el) - 1,
                              "element": r.element, "replicate": r.replicate, "orient": r.orient, "host": h,
                              "len_bp": len(el)})
                shift += len(el)
                prev = r.pos
            parts.append(s[prev:])
            new.append((d, "".join(parts)))
        write_fasta(out / "fasta" / f"{gid}.fna", new)
    t = pd.DataFrame(truth).merge(cand[["element", "name", "gc", "has_integrase", "size_bin"]], on="element")
    t.to_csv(out / "spike_truth.tsv", sep="\t", index=False)
    # in-situ arm: confident inovirus loci outside training, scored where they sit
    loci = pd.read_csv(a.inovirus_loci, sep="\t")
    ins = loci[(loci.evidence == "both") & (loci.n_tax >= 3) & (loci.split != "train")].copy()
    ins["len_bp"] = ins.end - ins.start
    ins = ins.rename(columns={"label": "genomad_label"})
    ins["start"] = ins.start + 1                                # loci are 0-based half-open; truth is 1-based
    ins[["genome_id", "contig", "start", "end", "len_bp", "genomad_label", "split", "in_eval"]].to_csv(
        out / "insitu_truth.tsv", sep="\t", index=False)
    (out / "insitu_genomes.txt").write_text("\n".join(sorted(ins.genome_id.unique())) + "\n")
    (out / "spike_genomes.txt").write_text("\n".join(sorted(t.genome_id.unique())) + "\n")
    print(f"[build] {len(t)} implants in {t.genome_id.nunique()} spiked genomes; by size "
          f"{t.size_bin.value_counts().to_dict()}; in-situ {len(ins)} loci in {ins.genome_id.nunique()} genomes", flush=True)


# ---------------------------------------------------------------- eval
def cmd_eval(a) -> None:
    from evalkit import MODE_THRESHOLDS, calls, fine_regions_by_genome, load_tracks, refine_all

    runs = dict(s.split("=", 1) for s in a.run)
    tracks = {k: load_tracks(v) for k, v in runs.items()}
    refines = [x for x in a.refine.split(",") if x]
    # fine tracks exist only around the calls of the mode they were scored in
    # (821410: Default at 0.995), so refinement variants are only meaningful then
    fine = pd.concat([pd.read_csv(p, sep="\t", usecols=["genome_id", "contig_id", "start", "end", "prob"])
                      for p in Path(runs["ens"]).glob("*/fine_windows.tsv")], ignore_index=True) \
        if any(x != "off" for x in refines) else None

    def overlaps_any(reg, g, c, s, e):
        return ((reg.genome_id == g) & (reg.contig_id == c) & (reg.end > s) & (reg.start < e)).any()

    def score(truth: pd.DataFrame, reg: pd.DataFrame) -> pd.DataFrame:
        out = []
        for r in truth.itertuples():
            s0, e0 = r.start - 1, r.end
            x = reg[(reg.genome_id == r.genome_id) & (reg.contig_id == r.contig) & (reg.end > s0) & (reg.start < e0)]
            iv = sorted((max(int(p), s0), min(int(q), e0)) for p, q in zip(x.start, x.end))
            bp, cs, ce = 0, None, None
            for p, q in iv:
                if ce is None or p > ce:
                    bp += (ce - cs) if ce is not None else 0
                    cs, ce = p, q
                else:
                    ce = max(ce, q)
            bp += (ce - cs) if ce is not None else 0
            L = e0 - s0
            iou = max([(min(int(q), e0) - max(int(p), s0)) / (max(int(q), e0) - min(int(p), s0))
                       for p, q in zip(x.start, x.end)], default=0.0)
            out.append({"cov": bp / L, "iou": iou})
        return pd.DataFrame(out, index=truth.index)

    spike = pd.read_csv(a.spike_truth, sep="\t")
    insitu = pd.read_csv(a.insitu_truth, sep="\t")
    res, per = {}, []
    fine_cache = {float(t): fine_regions_by_genome(fine, float(t), 0.1) for t in refines if t != "off"}
    for mode in ("exploratory", "default", "conservative"):
        base = calls(tracks["ens"], mode)
        seeds = {s: calls(tracks[s], mode) for s in ("r1", "r2", "r3") if s in tracks}
        for refine in refines:
            if mode != "default" and refine != "off":
                continue
            reg = base if refine == "off" else refine_all(base.assign(fasta_path="x"), fine_cache[float(refine)],
                                                          0.6, 4000, 0.05, (15000, 60000))
            sup = np.array([sum(overlaps_any(sr, r.genome_id, r.contig_id, r.start, r.end) for sr in seeds.values())
                            for r in reg.itertuples()]) if seeds else np.full(len(reg), 3)
            for k in (1, 2, 3):
                if k > 1 and not seeds:
                    continue
                v = reg[sup >= k]
                # the threshold is in the name: modes were re-derived 2026-09-29, and
                # 821410's names (0.95 / 0.995 / 0.999) must not be read as today's
                name = f"{mode}@{MODE_THRESHOLDS[mode]:g}/refine={refine}/vote>={k}"
                sc, si = score(spike, v), score(insitu, v)
                d = {"spike": {}, "insitu": {}}
                for b, idx in spike.groupby("size_bin").groups.items():
                    s = sc.loc[idx]
                    d["spike"][b] = {"n": int(len(s)), "detected": float((s["cov"] > 0).mean()),
                                     "covered50": float((s["cov"] >= 0.5).mean()), "median_iou": float(s["iou"].median())}
                d["spike"]["all"] = {"n": int(len(sc)), "detected": float((sc["cov"] > 0).mean()),
                                     "covered50": float((sc["cov"] >= 0.5).mean()), "median_iou": float(sc["iou"].median())}
                for lab, idx in insitu.groupby("genomad_label").groups.items():
                    s = si.loc[idx]
                    d["insitu"][lab] = {"n": int(len(s)), "detected": float((s["cov"] > 0).mean()),
                                        "covered50": float((s["cov"] >= 0.5).mean())}
                d["calls_per_spiked_genome"] = float(v[v.genome_id.isin(spike.genome_id)].groupby("genome_id").size().mean())
                res[name] = d
                per.append(spike.assign(variant=name, cov=sc["cov"].values, iou=sc["iou"].values))
                al = d["spike"]["all"]
                print(f"[{name}] spike detected {al['detected']:.3f} covered50 {al['covered50']:.3f} "
                      f"IoU {al['median_iou']:.3f} | " + " ".join(f"{b} {d['spike'][b]['detected']:.2f}/{d['spike'][b]['covered50']:.2f}"
                                                         for b in sorted(d['spike']) if b != 'all')
                      + " | insitu " + " ".join(f"{l} {x['detected']:.2f}" for l, x in d["insitu"].items()), flush=True)
    out = Path(a.out)
    out.write_text(json.dumps(res, indent=2) + "\n")
    pd.concat(per).to_csv(out.with_suffix(".elements.csv"), index=False)
    print(f"[done] wrote {out}", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("select")
    s.add_argument("--genomes_fna", required=True)
    s.add_argument("--proteins_faa", required=True)
    s.add_argument("--out_dir", required=True)
    s.add_argument("--min_len", type=int, default=1500)
    s.add_argument("--max_len", type=int, default=10000)
    s.add_argument("--dedupe", type=float, default=0.5, help="containment at or above which a genome is redundant")
    s.add_argument("--min_metagenome_len", type=int, default=4000,
                   help="metagenome-derived 'complete' genomes below this are implausible and dropped")
    d = sub.add_parser("disjoint")
    d.add_argument("--candidates_dir", required=True)
    d.add_argument("--beds_dir", required=True)
    d.add_argument("--sequences_dir", required=True)
    d.add_argument("--train_split", required=True)
    d.add_argument("--inovirus_loci", required=True)
    d.add_argument("--max_containment", type=float, default=0.1)
    b = sub.add_parser("build")
    b.add_argument("--candidates_dir", required=True)
    b.add_argument("--hosts", required=True, help="genome list; each needs <host_tracks>/<g>/windows.tsv")
    b.add_argument("--host_tracks", required=True, help="clean ensemble score tracks (10 kb)")
    b.add_argument("--host_elements", required=True, help="geNomad elements of the hosts (accession,contig,start,end)")
    b.add_argument("--genomes_dir", required=True)
    b.add_argument("--inovirus_loci", required=True)
    b.add_argument("--out_dir", required=True)
    b.add_argument("--replicates", type=int, default=2, help="host contexts per element (a maximum; "
                   "capacity permitting, every element gets one before any gets two)")
    # 80 kb apart keeps the gap between two implants' calls above DEFAULT_MAX_GAP_BP (50 kb),
    # so the fragment-merge stage cannot fuse two implants into one call; --element_gap does
    # the same for a real geNomad element. Measured 2026-09-28 on the 42-genome tracks:
    # 20 kb/80 kb gives 410 quiet sites, against 148 at 30 kb/150 kb.
    b.add_argument("--spacing", type=int, default=80000)
    b.add_argument("--quiet_flank", type=int, default=20000)
    b.add_argument("--quiet_prob", type=float, default=0.1)
    b.add_argument("--site_rule", choices=("model", "genomad"), default="model",
                   help="'model' (821410): sites must also be quiet on the host's score track; "
                        "'genomad': geNomad elements and spacing only, no model scores")
    b.add_argument("--element_gap", type=int, default=60000)
    b.add_argument("--seed", type=int, default=20260928)
    e = sub.add_parser("eval")
    e.add_argument("--spike_truth", required=True)
    e.add_argument("--insitu_truth", required=True)
    e.add_argument("--run", action="append", required=True, metavar="LABEL=DIR", help="r1..r3 and ens (ens needs fine tracks)")
    e.add_argument("--refine", default="off,0.8,0.95",
                   help="comma list of fine thresholds to score ('off' = none); only valid for the "
                        "mode the fine tracks were scored in")
    e.add_argument("--out", required=True)
    a = ap.parse_args()
    {"select": cmd_select, "disjoint": cmd_disjoint, "build": cmd_build, "eval": cmd_eval}[a.cmd](a)


if __name__ == "__main__":
    main()
