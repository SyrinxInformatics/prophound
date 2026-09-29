#!/usr/bin/env python3
import sys, os

sizes_path = "labels/all.contig.sizes"
pos_bed = "windows/bed/pos_win.cleaned.bed"
out_bed = "windows/bed/neg_win.hard.strict.bed"

hard_min = 1000
hard_max = 10000
buffer_bp = 200


def read_sizes(p):
    d = {}
    with open(p) as f:
        for line in f:
            if not line.strip():
                continue
            c, L = line.split()[:2]
            d[c] = int(L)
    return d


def read_bed(p):
    d = {}
    order = []
    with open(p) as f:
        for line in f:
            if not line.strip():
                continue
            c, s, e, id4 = line.rstrip("\n").split("\t")
            s = int(s)
            e = int(e)
            d.setdefault(c, []).append((s, e, id4))
            order.append((c, s, e, id4))
    for c in d:
        d[c].sort()
    return d, order


def overlaps(any_list, s, e):
    import bisect

    a = any_list
    i = bisect.bisect_left(a, (s, -1, ""))
    if i > 0:
        ps, pe, _ = a[i - 1]
        if not (pe <= s or e <= ps):
            return True
    if i < len(a):
        ns, ne, _ = a[i]
        if not (ne <= s or e <= ns):
            return True
    return False


contig_sizes = read_sizes(sizes_path)
pos_by_chr, pos_order = read_bed(pos_bed)
neg_by_chr = {c: [] for c in pos_by_chr.keys()}

out = []
for c, s, e, _ in pos_order:
    L = e - s
    if c not in contig_sizes:
        continue
    clen = contig_sizes[c]
    left_min = max(0, s - hard_max - L)
    left_max = max(-1, s - hard_min - L)
    right_min = min(clen - L, e + hard_min)
    right_max = min(clen - L, e + hard_max)

    placed = False
    if left_max >= left_min:
        start = left_max
        while start >= left_min:
            ss = start
            ee = ss + L
            if ee <= s - buffer_bp:
                if not overlaps(pos_by_chr.get(c, []), ss, ee) and not overlaps(
                    neg_by_chr.get(c, []), ss, ee
                ):
                    neg_by_chr[c].append((ss, ee, f"negH|{len(out)+1}"))
                    neg_by_chr[c].sort()
                    out.append((c, ss, ee, f"negH|{len(out)+1}"))
                    placed = True
                    break
            start -= max(1, L // 256)
    if not placed and right_max >= right_min:
        start = right_min
        while start <= right_max:
            ss = start
            ee = ss + L
            if ss >= e + buffer_bp:
                if not overlaps(pos_by_chr.get(c, []), ss, ee) and not overlaps(
                    neg_by_chr.get(c, []), ss, ee
                ):
                    neg_by_chr[c].append((ss, ee, f"negH|{len(out)+1}"))
                    neg_by_chr[c].sort()
                    out.append((c, ss, ee, f"negH|{len(out)+1}"))
                    placed = True
                    break
            start += max(1, L // 256)

with open(out_bed, "w") as w:
    for c, s, e, id4 in out:
        w.write(f"{c}\t{s}\t{e}\t{id4}\n")

print(f"requested={len(pos_order)} placed_hard={len(out)} out={out_bed}")

