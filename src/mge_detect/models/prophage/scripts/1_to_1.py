#!/usr/bin/env python3
import os, sys, random, hashlib, tempfile, shutil, subprocess

pos_bed = "windows/bed/pos_win.cleaned.bed"
neg_bed = "windows/bed/neg_win.cleaned.bed"
sizes_fp = "labels/all.contig.sizes"
all_fa = "raw/all_genomes.fa"

out_neg = "windows/bed/neg_win.balanced.bed"
out_all_bed = "windows/bed/all_win.contig.balanced.bed"
out_fa = "windows/fasta/all_512k.balanced.fa"
out_lab = "windows/labels/all_512k.balanced.labels.tsv"
out_train = "splits/train.contig.balanced.bed"
out_val = "splits/val.contig.balanced.bed"
out_test = "splits/test.contig.balanced.bed"
manifest = "windows/manifest/all_512k_manifest_balanced.csv"

for f in [pos_bed, neg_bed, all_fa]:
    if not os.path.exists(f) or os.path.getsize(f) == 0:
        print(f"MISSING {f}", file=sys.stderr)
        sys.exit(1)

if not os.path.exists(sizes_fp) or os.path.getsize(sizes_fp) == 0:
    if not os.path.exists(all_fa + ".fai"):
        subprocess.check_call(["samtools", "faidx", all_fa])
    with open(sizes_fp, "w") as w, open(all_fa + ".fai") as r:
        for line in r:
            c, L = line.split("\t", 2)[:2]
            w.write(f"{c}\t{L.strip()}\n")


def read_bed(path):
    out = []
    with open(path) as r:
        for ln in r:
            f = ln.rstrip("\n").split("\t")
            if len(f) >= 3:
                c, s, e = f[0], int(f[1]), int(f[2])
                name = f[3] if len(f) > 3 else ""
                out.append((c, s, e, name))
    return out


def write_bed(path, rows):
    with open(path, "w") as w:
        for c, s, e, n in rows:
            w.write(f"{c}\t{s}\t{e}\t{n}\n")


sizes = {}
with open(sizes_fp) as r:
    for ln in r:
        f = ln.split()
        if len(f) >= 2:
            sizes[f[0]] = int(f[1])

pos = read_bed(pos_bed)
neg = read_bed(neg_bed)


def merge_intervals(ivs):
    if not ivs:
        return []
    ivs = sorted(ivs, key=lambda x: (x[0], x[1], x[2]))
    out = []
    cur_c, cs, ce = ivs[0][0], ivs[0][1], ivs[0][2]
    for c, s, e, _ in ivs[1:]:
        if c != cur_c or s > ce:
            out.append((cur_c, cs, ce))
            cur_c, cs, ce = c, s, e
        else:
            if e > ce:
                ce = e
    out.append((cur_c, cs, ce))
    return out


forbidden = merge_intervals([(c, s, e, "") for (c, s, e, _) in pos + neg])

by_contig = {}
for c, L in sizes.items():
    by_contig[c] = []
for c, s, e in forbidden:
    by_contig.setdefault(c, []).append((s, e))
for c in by_contig:
    by_contig[c].sort()


def subtract_forbidden(c):
    L = sizes[c]
    out = []
    cur = 0
    for s, e in by_contig[c]:
        if s > cur:
            out.append((cur, s))
        if e > cur:
            cur = e
        if cur >= L:
            break
    if cur < L:
        out.append((cur, L))
    return out


allowed = {c: subtract_forbidden(c) for c in sizes}

random.seed(1337)

need = max(0, len(pos) - len(neg))
new_neg = []
if need > 0:
    pos_lengths = [e - s for (_, s, e, _) in pos]
    contigs_with_space = [
        c for c in sizes if any(b - a >= min(pos_lengths) for a, b in allowed[c])
    ]
    if not contigs_with_space:
        print("INSUFFICIENT allowed negatives: no room anywhere", file=sys.stderr)
        sys.exit(2)

    def place_window(length):
        candidates = []
        for c in contigs_with_space:
            for a, b in allowed[c]:
                span = b - a
                if span >= length:
                    candidates.append((c, a, b))
        if not candidates:
            return None
        c, a, b = random.choice(candidates)
        max_start = b - length
        if max_start < a:
            return None
        start = random.randint(a, max_start)
        end = start + length
        return (c, start, end)

    def insert_forbidden(c, s, e):
        ivs = by_contig[c]
        new = []
        placed = False
        for a, b in ivs:
            if e < a and not placed:
                new.append((s, e))
                placed = True
            new.append((a, b))
        if not placed:
            new.append((s, e))
        merged = []
        for a, b in new:
            if not merged or a > merged[-1][1]:
                merged.append((a, b))
            else:
                if b > merged[-1][1]:
                    merged[-1] = (merged[-1][0], b)
        by_contig[c] = merged
        allowed[c] = subtract_forbidden(c)

    placed = 0
    next_id = 1
    tries = 0
    cap = need * 200
    while placed < need and tries < cap:
        L = pos_lengths[placed % len(pos_lengths)]
        hit = place_window(L)
        tries += 1
        if hit is None:
            continue
        c, s, e = hit
        new_neg.append((c, s, e, f"negB|{next_id}"))
        next_id += 1
        insert_forbidden(c, s, e)
        placed += 1
    if placed < need:
        print(
            f"INSUFFICIENT allowed negatives: need {need} got {placed}", file=sys.stderr
        )
        sys.exit(2)

balanced_neg = sorted(neg + new_neg, key=lambda x: (x[0], x[1], x[2]))
write_bed(out_neg, balanced_neg)

all_rows = sorted(pos + balanced_neg, key=lambda x: (x[0], x[1], x[2]))
write_bed(out_all_bed, all_rows)


def self_overlap_count(rows):
    cprev, sp, ep, _ = rows[0]
    count = 0
    for c, s, e, _ in rows[1:]:
        if c == cprev and s < ep:
            count += 1
            ep = max(ep, e)
        else:
            cprev, sp, ep = c, s, e
    return count


ov = self_overlap_count(all_rows)
if ov != 0:
    print("OVERLAPS_PRESENT", file=sys.stderr)
    sys.exit(3)

subprocess.check_call(
    ["bedtools", "getfasta", "-fi", all_fa, "-bed", out_all_bed, "-name", "-fo", out_fa]
)

with open(out_lab, "w") as w:
    for c, s, e, n in all_rows:
        w.write(f"{n}\t{1 if n.startswith('pos|') else 0}\n")


def split_of_contig(c):
    h = hashlib.md5(c.encode()).hexdigest()
    x = int(h[:8], 16) / 0xFFFFFFFF
    return "train" if x < 0.803 else ("val" if x < 0.909 else "test")


tw, vw, ew = open(out_train, "w"), open(out_val, "w"), open(out_test, "w")
for row in all_rows:
    c, s, e, n = row
    sp = split_of_contig(c)
    (tw if sp == "train" else vw if sp == "val" else ew).write(f"{c}\t{s}\t{e}\t{n}\n")
tw.close()
vw.close()
ew.close()

with open(manifest, "w") as w:
    w.write("id,chrom,start,end,split,label\n")
    for c, s, e, n in all_rows:
        sp = split_of_contig(c)
        labv = 1 if n.startswith("pos|") else 0
        w.write(f"{n},{c},{s},{e},{sp},{labv}\n")

n_all = len(all_rows)
n_pos = sum(1 for _ in filter(lambda r: r[3].startswith("pos|"), all_rows))
print(f"balanced_total={n_all} pos={n_pos} neg={n_all-n_pos}")

