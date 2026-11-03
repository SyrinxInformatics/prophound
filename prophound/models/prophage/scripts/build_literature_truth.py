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
"""Literature-curated prophage truth for the 0.1.0 pilot (docs/final-benchmark.md).

Every element comes from the primary literature or from the original authors'
own curation in the genome record, never from a detection tool. The basis is
recorded per element:

  author_feature  a misc_feature the genome paper's authors placed in the
                  original GenBank record (EDL933 CP-933*/BP-933W, Sakai
                  Sp/SpLE, SF370 attP-L..attP-R)
  author_genes    the span of the densest cluster of genes the authors
                  annotated as belonging to the element (EDL933 notes, B.
                  subtilis 168 reannotation notes)
  locus_tags      a locus-tag range published in the literature (PAO1, N16961)
  phage_genome    the published phage genome anchored on the chromosome by
                  exact 31-mers (S. aureus 8325 phi11/12/13, Iandolo et al. 2002)

class: `prophage` (literature calls it a prophage), `cryptic_prophage` (Perna
et al.'s term for EDL933's CP-933 set), `remnant` (described as a remnant,
defective or phage-like element), `phage_derived` (phage-derived but not a
prophage: the TLC satellite and the R2/F2 pyocin, a phage-tail bacteriocin).
Classes were split this way 2026-09-29 (user; docs/research-integrity.md C4):
the release headline uses `prophage` only and reports the others apart.

Source coordinates are lifted onto the corpus RefSeq assemblies by sequence:
the build checks that each source record is identical to the RefSeq contig
(or anchors the source boundaries uniquely when it is not), and stores 30 bp
flanks at each boundary so the evaluation can verify its own FASTAs.

Inputs are downloaded by the caller into --src (NCBI GenBank records, RefSeq
FASTA/GFF from NCBI Datasets, phage genomes); see docs/final-benchmark.md.

Example:
  python scripts/build_literature_truth.py --src /tmp/lit --out scripts/literature_truth.tsv
"""

from __future__ import annotations

import argparse
import glob
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd

GENOMES = {  # accession: (corpus genome_id, split, organism, source GenBank record or None)
    "GCF_000006665.1": ("GCF_000006665.1_ASM666v1_genomic", "test", "E. coli O157:H7 EDL933", "AE005174"),
    "GCF_000006785.2": ("GCF_000006785.2_ASM678v2_genomic", "test", "S. pyogenes SF370", "AE004092"),
    "GCF_000013425.1": ("GCF_000013425.1_ASM1342v1_genomic", "test", "S. aureus NCTC 8325", None),
    "GCF_000006745.1": ("GCF_000006745.1_ASM674v1_genomic", "test", "V. cholerae N16961", None),
    "GCF_000008865.2": ("GCF_000008865.2_ASM886v2_genomic", "val", "E. coli O157:H7 Sakai", None),
    "GCF_000009045.1": ("GCF_000009045.1_ASM904v1_genomic", "val", "B. subtilis 168", "AL009126"),
    "GCF_000006765.1": ("GCF_000006765.1_ASM676v1_genomic", "val", "P. aeruginosa PAO1", None),
}
REF = {
    "EDL933": "Perna et al. 2001 Nature 409:529 (GenBank AE005174 author annotation)",
    "SF370": "Ferretti et al. 2001 PNAS 98:4658 (AE004092 attP features); Canchaya et al. 2002 Virology 302:245",
    "8325": "Iandolo et al. 2002 Gene 289:109 (phage genomes AF424781-3)",
    "N16961": "Waldor & Mekalanos 1996 Science 272:1910; Rubin et al. 1998 Mol Microbiol 28:1247 (TLC = VC1466-VC1470)",
    "Sakai": "Hayashi et al. 2001 DNA Res 8:11 (author misc_features in BA000007/NC_002695)",
    "Bsub": "Kunst et al. 1997 Nature 390:249; Barbe et al. 2009 Microbiology 155:1758 (AL009126 notes)",
    "PAO1_Pf4": "Webb et al. 2004 J Bacteriol 186:8066 (Pf4 = PA0715-PA0729)",
    "PAO1_pyo": "Nakayama et al. 2000 Mol Microbiol 38:213 (R2/F2 pyocin cluster PA0610-PA0648)",
}


# ---------------------------------------------------------------- readers
def fasta(path) -> dict:
    seqs, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        else:
            seqs[name].append(line.strip().upper())
    return {k: "".join(v) for k, v in seqs.items()}


def gb_features(path):
    """(type, start, end, qualifiers) for every feature, plus the sequence."""
    feats, seq, ftype, loc, quals, key, mode = [], [], None, "", {}, None, None

    def flush():
        nums = [int(x) for x in re.findall(r"\d+", loc)]
        if ftype and nums:
            feats.append((ftype, min(nums), max(nums), dict(quals)))

    for line in open(path):
        if line.startswith("FEATURES"):
            mode = "feat"
            continue
        if line.startswith("ORIGIN"):
            flush()
            ftype, mode = None, "seq"
            continue
        if mode == "seq":
            seq.append(re.sub(r"[^a-zA-Z]", "", line).upper())
        elif mode == "feat":
            if len(line) > 5 and line[5] != " ":
                flush()
                ftype, loc, quals, key = line[5:21].strip(), line[21:].strip(), {}, None
            elif line[21:22] == "/" and re.match(r"/(\w+)(=|$)", line[21:].strip()):
                m = re.match(r"/(\w+)=?(.*)", line[21:].strip())
                key = m.group(1)
                quals[key] = m.group(2).strip('"')
            elif key is None:
                loc += line.strip()
            else:
                quals[key] = (quals[key] + " " + line.strip().strip('"')).strip()
    return feats, "".join(seq)


def gff_genes(path) -> pd.DataFrame:
    rows = []
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] not in ("gene", "pseudogene"):
            continue
        a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        tags = {a.get("locus_tag", "")} | set(a.get("old_locus_tag", "").replace("%2C", ",").split(","))
        rows.append({"contig": f[0], "start": int(f[3]), "end": int(f[4]), "tags": tags})
    return pd.DataFrame(rows)


def tag_span(genes: pd.DataFrame, first: str, last: str):
    a = genes[genes.tags.map(lambda t: first in t)]
    b = genes[genes.tags.map(lambda t: last in t)]
    if len(a) != 1 or len(b) != 1 or a.contig.iloc[0] != b.contig.iloc[0]:
        raise SystemExit(f"locus tags {first}/{last} not uniquely found")
    return a.contig.iloc[0], min(a.start.iloc[0], b.start.iloc[0]), max(a.end.iloc[0], b.end.iloc[0])


def densest_cluster(spans, gap=20000):
    spans = sorted(spans)
    cl, cur = [], [spans[0]]
    for s in spans[1:]:
        if s[0] - max(x[1] for x in cur) > gap:
            cl.append(cur)
            cur = [s]
        else:
            cur.append(s)
    cl.append(cur)
    best = max(cl, key=len)
    return best[0][0], max(x[1] for x in best), len(best)


def lift(src_seq: str, tgt: dict, s: int, e: int):
    """Map 1-based [s, e] on src_seq to (contig, s, e) on the target assembly."""
    for c, t in tgt.items():
        if t == src_seq:
            return c, s, e, "identical"
    out = []
    for pos, anchor in ((s, src_seq[s - 1:s + 29]), (e, src_seq[e - 30:e])):
        hits = [(c, m.start()) for c, t in tgt.items() for m in re.finditer(anchor, t)]
        if len(hits) != 1:
            raise SystemExit(f"boundary {pos} anchors {len(hits)} times")
        out.append(hits[0])
    (c1, p1), (c2, p2) = out
    assert c1 == c2
    return c1, p1 + 1, p2 + 30, "anchored"




def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--src", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    src = Path(args.src)
    rows = []

    def add(acc, contig, s, e, name, cls, basis, ref, rec="", ss=None, se=None):
        rows.append({"accession": acc, "contig": contig, "start": int(s), "end": int(e), "name": name,
                     "class": cls, "basis": basis, "reference": ref, "source_record": rec,
                     "source_start": ss, "source_end": se})

    def tgt(acc):
        return fasta(glob.glob(str(src / f"fna_{acc}" / "**" / "*.fna"), recursive=True)[0])

    def gff(acc):
        return gff_genes(glob.glob(str(src / acc / "**" / "genomic.gff"), recursive=True)[0])

    # --- EDL933 ---
    feats, sq = gb_features(src / "AE005174.gb")
    T = tgt("GCF_000006665.1")
    marked, genes = {}, defaultdict(list)
    for t, s, e, q in feats:
        note = q.get("note", "")
        m = re.search(r"[Pp]rophage (CP-933\w|BP-933W)", note)
        if t == "misc_feature" and m:
            marked[m.group(1)] = (s, e)
        if t == "CDS":
            for n in set(re.findall(r"(CP-933\w|BP-933W)", " ".join(q.get(k, "") for k in ("note", "product")))):
                genes[n].append((s, e))
    for name in sorted(set(marked) | set(genes)):
        # CP-933P's author feature is its whole O-island (58 kb, non-phage genes
        # included); since 2026-09-29 (user) it takes the gene-cluster rule like
        # the other gene-annotated cryptic prophages
        if name in marked and name != "CP-933P":
            s, e = marked[name]
            basis = "author_feature"
        else:
            s, e, n = densest_cluster(genes[name])
            basis = f"author_genes(n={n})"
        c, ts, te, how = lift(sq, T, s, e)
        add("GCF_000006665.1", c, ts, te, name, "prophage" if name == "BP-933W" else "cryptic_prophage",
            f"{basis};lift={how}", REF["EDL933"], "AE005174", s, e)

    # --- SF370 ---
    feats, sq = gb_features(src / "AE004092.gb")
    T = tgt("GCF_000006785.2")
    att = defaultdict(list)
    for t, s, e, q in feats:
        m = re.search(r"bacteriophage (370\.\d) - attP", q.get("note", ""))
        if t == "misc_feature" and m:
            att[m.group(1)] += [s, e]
    for name, v in sorted(att.items()):
        c, ts, te, how = lift(sq, T, min(v), max(v))
        add("GCF_000006785.2", c, ts, te, f"SF{name}", "remnant" if name == "370.4" else "prophage",
            f"author_feature;lift={how}", REF["SF370"], "AE004092", min(v), max(v))

    # --- B. subtilis 168 ---
    feats, sq = gb_features(src / "AL009126.gb")
    T = tgt("GCF_000009045.1")
    lab = {"SPbeta": r"SP ?beta|SP-beta|SPß", "PBSX": r"PBSX", "skin": r"skin element",
           "prophage 1": r"prophage (region )?1\b", "prophage 3": r"prophage (region )?3\b",
           "prophage 6": r"prophage (region )?6\b"}
    for name, rx in lab.items():
        spans = [(s, e) for t, s, e, q in feats if t in ("CDS", "gene", "misc_feature")
                 and re.search(rx, " ".join(q.get(k, "") for k in ("note", "product", "function")), re.I)]
        if name == "SPbeta":   # a second cluster near 1.9 Mb only cites similarity to SPbeta
            spans = [x for x in spans if x[0] > 2_100_000]
        s, e, n = densest_cluster(spans)
        c, ts, te, how = lift(sq, T, s, e)
        add("GCF_000009045.1", c, ts, te, name, "prophage" if name == "SPbeta" else "remnant",
            f"author_genes(n={n});lift={how}", REF["Bsub"], "AL009126", s, e)

    # --- Sakai: author misc_features carried in the RefSeq record ---
    for line in open(glob.glob(str(src / "GCF_000008865.2" / "**" / "genomic.gff"), recursive=True)[0]):
        f = line.rstrip("\n").split("\t")
        if len(f) > 8 and f[2] == "sequence_feature":
            m = re.search(r"product=(prophage (Sp\d+)|phage-like element (SpLE\d))", f[8])
            if m:
                name = m.group(2) or m.group(3)
                add("GCF_000008865.2", f[0], f[3], f[4], name, "prophage" if name.startswith("Sp") and
                    not name.startswith("SpLE") else "remnant", "author_feature;lift=same_record",
                    REF["Sakai"], f[0], f[3], f[4])

    # --- PAO1 ---
    G = gff("GCF_000006765.1")
    for name, a, b, cls, ref in [("Pf4", "PA0715", "PA0729", "prophage", REF["PAO1_Pf4"]),
                                 ("R2/F2 pyocin", "PA0610", "PA0648", "phage_derived", REF["PAO1_pyo"])]:
        c, s, e = tag_span(G, a, b)
        add("GCF_000006765.1", c, s, e, name, cls, f"locus_tags {a}-{b};lift=same_record", ref)

    # --- N16961 ---
    G = gff("GCF_000006745.1")
    for name, a, b, cls in [("RS1-CTXphi array", "VC1452", "VC1465", "prophage"),
                            ("TLC", "VC1466", "VC1470", "phage_derived")]:
        c, s, e = tag_span(G, a, b)
        add("GCF_000006745.1", c, s, e, name, cls, f"locus_tags {a}-{b};lift=same_record", REF["N16961"])

    # --- S. aureus 8325: published phage genomes anchored by exact 31-mers ---
    T = tgt("GCF_000013425.1")
    comp = str.maketrans("ACGT", "TGCA")
    for name, acc in [("phi11", "AF424781"), ("phi12", "AF424782"), ("phi13", "AF424783")]:
        P = next(iter(fasta(src / f"{acc}.fa").values()))
        ks = {P[i:i + 31] for i in range(len(P) - 30)}
        ks |= {k.translate(comp)[::-1] for k in list(ks)}
        best = None
        for c, t in T.items():
            hits = [i for i in range(len(t) - 30) if t[i:i + 31] in ks]
            if hits:
                s, e, n = densest_cluster([(h, h + 31) for h in hits], gap=2000)
                if best is None or n > best[3]:
                    best = (c, s + 1, e, n)
        c, s, e, n = best
        frac = n / (e - s + 1 - 30)
        add("GCF_000013425.1", c, s, e, name, "prophage", f"phage_genome {acc} (31-mer cover {frac:.3f})",
            REF["8325"], acc)

    df = pd.DataFrame(rows)
    # boundary flanks from the corpus assembly, for verification against the evaluation's FASTAs
    flank = []
    for r in df.itertuples():
        t = tgt(r.accession)[r.contig]
        flank.append((t[r.start - 1:r.start + 29], t[r.end - 30:r.end]))
    df["left_flank"], df["right_flank"] = zip(*flank)
    df.insert(0, "genome_id", df.accession.map(lambda a: GENOMES[a][0]))
    df.insert(1, "split", df.accession.map(lambda a: GENOMES[a][1]))
    df.insert(2, "organism", df.accession.map(lambda a: GENOMES[a][2]))
    df["len_bp"] = df.end - df.start + 1
    df = df.sort_values(["split", "genome_id", "contig", "start"])
    df.to_csv(args.out, sep="\t", index=False)
    print(df.groupby(["split", "organism", "class"]).size().to_string())
    print(f"[done] {len(df)} elements -> {args.out}")


if __name__ == "__main__":
    main()
