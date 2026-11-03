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
"""Base-exact call boundaries from attL/attR direct repeats, prototyped offline
(docs/att-boundaries.md).

Integration by site-specific recombination leaves the element between two
copies of the att core, so a direct repeat with one copy near each edge of a
call marks both boundaries to the base.

Subcommands:
  run    calls rebuilt as `predict` ships them (evalkit.shipped_calls) for each
         mode, plus decoys: the same geometry shifted into host sequence away from
         every call and truth element. For each, every exact direct repeat >= 16 bp
         between the widest zones is recorded with its copy count and entropy, so
         `score` can sweep the selection rules without searching again. Writes
         <out>.items.tsv (one row per call or decoy) and <out>.cands.tsv.
  score  applies each selection rule in a grid and reports, per mode: att rate on
         calls and decoys, covered50/covered/matched recall, IoU, bp called and
         off-element, and edge errors against the truth (bp). Writes <out>.csv.

A selection rule: both copies inside zones [s - w_out, s + w_in] and
[e - w_in, e + w_out] (stopped at the call's midpoint); core length <= max_len;
<= max_copies copies in the genome (both strands); base entropy >= min_entropy;
the att interval must contain the call's high-scoring core (centres of windows
at or above the mode's seed threshold) up to core_slack bp; then the candidate
with the smallest chance expectation E = m * n * (1 - p) * p**L, accepted if
E <= e_max. m, n are the zone lengths and p = sum of squared base frequencies.

Example:
  python scripts/att_boundaries.py run --tracks analysis/ngap_dose/impute__clean__0 \\
      --genomes_dir ../../../arc_data/raw/genomes --elements scripts/eval_sizestrat_elements.csv \\
      --out analysis/att/tune42
  python scripts/att_boundaries.py score --run analysis/att/tune42 \\
      --elements scripts/eval_sizestrat_elements.csv --out analysis/att/tune42_sweep
"""

from __future__ import annotations

import argparse
import itertools
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import pyfaidx

CODE = np.full(256, 255, dtype=np.uint8)
for _i, _b in enumerate(b"ACGT"):
    CODE[_b] = _i
    CODE[_b + 32] = _i
COMP = str.maketrans("ACGTacgt", "TGCAtgca")

SEED_K = 12          # seed length; every exact repeat >= SEED_K is found
MAX_SEED_COPIES = 8  # a seed this frequent within one zone is low complexity / tandem
MIN_RECORD_LEN = 16  # shorter repeats are never significant in kb-sized zones
MIN_ELEMENT = 2000   # the att pair must enclose at least this much sequence
W_IN_MAX, W_OUT_MAX = 25000, 10000   # widest zones searched; `score` narrows them


def encode(seq: str) -> np.ndarray:
    return CODE[np.frombuffer(seq.encode(), dtype=np.uint8)]


def seed_table(a: np.ndarray, k: int = SEED_K) -> pd.DataFrame:
    """Every valid k-mer in `a` (codes 0..3, 255 = other) with its offset."""
    n = len(a) - k + 1
    if n <= 0:
        return pd.DataFrame({"code": np.zeros(0, np.uint64), "pos": np.zeros(0, np.int64)})
    bad = np.concatenate([[0], np.cumsum(a == 255)])
    ok = (bad[k:k + n] - bad[:n]) == 0
    code = np.zeros(n, dtype=np.uint64)
    a64 = np.where(a == 255, 0, a).astype(np.uint64)
    for i in range(k):
        code = code * np.uint64(4) + a64[i:i + n]
    t = pd.DataFrame({"code": code[ok], "pos": np.nonzero(ok)[0]})
    counts = t.code.map(t.code.value_counts())
    return t[counts <= MAX_SEED_COPIES]


def direct_repeats(left: str, right: str, k: int = SEED_K) -> pd.DataFrame:
    """Maximal exact direct repeats between two strings, length >= k.

    Columns: li, ri (start offsets in left/right) and length.
    """
    a, b = seed_table(encode(left), k), seed_table(encode(right), k)
    hits = a.merge(b, on="code", suffixes=("_l", "_r"))
    if hits.empty:
        return pd.DataFrame(columns=["li", "ri", "length"])
    hits["diag"] = hits.pos_r - hits.pos_l
    hits = hits.sort_values(["diag", "pos_l"])
    brk = (hits.diag.diff() != 0) | (hits.pos_l.diff() != 1)
    hits["run"] = brk.cumsum()
    g = hits.groupby("run").agg(li=("pos_l", "min"), ri=("pos_r", "min"), n=("pos_l", "size"))
    return pd.DataFrame({"li": g.li.to_numpy(), "ri": g.ri.to_numpy(), "length": (g.n + k - 1).to_numpy()})


def entropy_bits(s: str) -> float:
    c = np.array([s.count(x) for x in "ACGT"], dtype=float)
    p = c[c > 0] / c.sum()
    return float(-(p * np.log2(p)).sum())


def zones(s: int, e: int, contig_len: int, w_in: int, w_out: int) -> tuple[int, int, int, int]:
    """Left and right search zones for a call [s, e), each stopped at its midpoint."""
    mid = (s + e) // 2
    return max(0, s - w_out), min(mid, s + w_in), max(mid, e - w_in), min(contig_len, e + w_out)


def find_candidates(contig: str, genome_text: str, s: int, e: int,
                    w_in: int = W_IN_MAX, w_out: int = W_OUT_MAX) -> list:
    """Every direct repeat >= MIN_RECORD_LEN with one copy in each zone of call [s, e).

    `genome_text` is every contig joined, for counting copies on both strands.
    """
    l0, l1, r0, r1 = zones(s, e, len(contig), w_in, w_out)
    if l1 - l0 < SEED_K or r1 - r0 < SEED_K:
        return []
    reps = direct_repeats(contig[l0:l1], contig[r0:r1])
    reps = reps[reps.length >= MIN_RECORD_LEN]
    out = []
    for r in reps.itertuples():
        start, end = l0 + int(r.li), r0 + int(r.ri) + int(r.length)
        if end - start < MIN_ELEMENT:
            continue
        core = contig[start:start + int(r.length)].upper()
        rc = core.translate(COMP)[::-1]
        copies = genome_text.count(core) + (genome_text.count(rc) if rc != core else 0)
        out.append({"att_start": start, "att_end": end, "att_len": int(r.length), "att_copies": int(copies),
                    "entropy": entropy_bits(core), "att_core": core if len(core) <= 300 else core[:300] + "..."})
    return out


# --- run


def load_genome(path: Path) -> tuple[dict, str, float]:
    fa = pyfaidx.Fasta(str(path), sequence_always_upper=True, rebuild=False)
    seqs = {k: str(fa[k][:]) for k in fa.keys()}
    text = "|".join(seqs.values())
    f = np.array([text.count(x) for x in "ACGT"], dtype=float)
    f /= f.sum()
    return seqs, text, float((f ** 2).sum())


def decoy_offsets(rng, s: int, e: int, length: int, blocked: list, pad: int, n: int) -> list:
    """Up to `n` shifts d for which [s+d - pad, e+d + pad] fits the contig and avoids `blocked`."""
    span0, span1 = s - pad, e + pad
    out = []
    for _ in range(50 * n):
        if len(out) == n:
            break
        a = int(rng.integers(0, max(1, length - (span1 - span0))))
        d = a - span0
        lo, hi = span0 + d, span1 + d
        if lo < 0 or hi > length or any(lo < b1 and hi > b0 for b0, b1 in blocked):
            continue
        out.append(d)
    return out


def call_core(win: pd.DataFrame, contig: str, s: int, e: int, thr: float) -> tuple[int, int]:
    """Span of window centres inside [s, e) scoring >= thr (the call's midpoint if none)."""
    w = win[win.contig_id == contig]
    c = (w.start + w.end) // 2
    c = c[(c >= s) & (c < e) & (w.prob >= thr)]
    return (int(c.min()), int(c.max())) if len(c) else ((s + e) // 2, (s + e) // 2)


def genome_job(job: dict) -> tuple[list, list]:
    g, calls, win, truth = job["genome"], job["calls"], job["windows"], job["truth"]
    seqs, text, p = load_genome(Path(job["genomes_dir"]) / f"{g}.fna")
    rng = np.random.default_rng([job["seed"], int.from_bytes(g.encode()[-8:], "little")])
    items, cands = [], []
    for c in calls.itertuples():
        contig = seqs[c.contig_id]
        edge = getattr(c, "contig_edge", "none")
        searchable = not bool(getattr(c, "short", False)) and (pd.isna(edge) or edge == "none")
        thr = job["seed_thr"][c.mode]
        lo, hi = call_core(win, c.contig_id, int(c.start), int(c.end), thr)
        base = {"mode": c.mode, "genome_id": g, "contig_id": c.contig_id, "contig_len": len(contig),
                "call_start": int(c.start), "call_end": int(c.end), "p_match": p,
                "contig_edge": edge, "searchable": searchable}
        blocked = [(int(x.start) - W_OUT_MAX, int(x.end) + W_OUT_MAX)
                   for x in calls[calls.contig_id == c.contig_id].itertuples()]
        blocked += [(int(x.start) - W_OUT_MAX, int(x.end) + W_OUT_MAX)
                    for x in truth[truth.contig == c.contig_id].itertuples()]
        shifts = [0] + (decoy_offsets(rng, int(c.start), int(c.end), len(contig), blocked, W_OUT_MAX, job["n_decoys"])
                        if searchable else [])
        for k, d in enumerate(shifts):
            item_id = f"{c.mode}|{g}|{c.contig_id}|{c.start}|{k}"
            s, e = int(c.start) + d, int(c.end) + d
            items.append({**base, "item_id": item_id, "kind": "call" if k == 0 else "decoy", "shift": d,
                          "s": s, "e": e, "core_lo": lo + d, "core_hi": hi + d})
            if searchable:
                cands += [{"item_id": item_id, **x} for x in find_candidates(contig, text, s, e)]
    return items, cands


def cmd_run(args) -> None:
    import evalkit as ek

    ens = ek.load_tracks(args.tracks)
    seeds = ([ek.load_tracks(d) for d in args.seed_tracks] if args.seed_tracks
             else [ek.load_tracks(args.tracks, f"prob_ckpt{k}") for k in (1, 2, 3)])
    truth = pd.read_csv(args.elements) if args.elements else pd.DataFrame(columns=["accession", "contig", "start", "end"])
    per_mode = []
    for mode in args.modes:
        reg = ek.shipped_calls(ens, seeds, mode, args.genomes_dir)
        per_mode.append(reg.assign(mode=mode))
        print(f"[att] {mode}: {len(reg)} calls", flush=True)
    reg = pd.concat(per_mode, ignore_index=True)
    seed_thr = {m: ek.MODE_THRESHOLDS[m] for m in args.modes}
    jobs = [{"genome": g, "calls": sub, "windows": ens[g][0], "truth": truth[truth.accession == g],
             "genomes_dir": args.genomes_dir, "seed": args.seed, "n_decoys": args.n_decoys, "seed_thr": seed_thr}
            for g, sub in reg.groupby("genome_id", sort=True)]
    items, cands = [], []
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        for i, (it, ca) in enumerate(ex.map(genome_job, jobs)):
            items += it
            cands += ca
            print(f"[att] {i + 1}/{len(jobs)} genomes", end="\r", flush=True)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(items).to_csv(f"{args.out}.items.tsv", sep="\t", index=False)
    pd.DataFrame(cands).to_csv(f"{args.out}.cands.tsv", sep="\t", index=False)
    print(f"\n[att] wrote {args.out}.items.tsv ({len(items)}) and {args.out}.cands.tsv ({len(cands)})")


# --- selection (the rule a shipped version would implement)


def select(items: pd.DataFrame, cands: pd.DataFrame, w_in: int, w_out: int, max_len: int, max_copies: int,
           min_entropy: float, core_slack: int | None, e_max: float) -> pd.DataFrame:
    """Best accepted candidate per item (item_id, att_start, att_end, att_len, evalue)."""
    it = items[items.searchable].copy()
    mid = (it.s + it.e) // 2
    it["l0"], it["l1"] = np.maximum(0, it.s - w_out), np.minimum(mid, it.s + w_in)
    it["r0"], it["r1"] = np.maximum(mid, it.e - w_in), np.minimum(it.contig_len, it.e + w_out)
    c = cands.merge(it[["item_id", "l0", "l1", "r0", "r1", "core_lo", "core_hi", "p_match"]], on="item_id")
    ok = ((c.att_start >= c.l0) & (c.att_start + c.att_len <= c.l1)
          & (c.att_end - c.att_len >= c.r0) & (c.att_end <= c.r1)
          & (c.att_len <= max_len) & (c.att_copies <= max_copies) & (c.entropy >= min_entropy))
    if core_slack is not None:
        ok &= (c.att_start <= c.core_lo + core_slack) & (c.att_end >= c.core_hi - core_slack)
    c = c[ok].copy()
    m, n = (c.l1 - c.l0).clip(lower=1), (c.r1 - c.r0).clip(lower=1)
    c["evalue"] = m * n * (1 - c.p_match) * c.p_match ** c.att_len
    c = c[c.evalue <= e_max]
    return c.sort_values(["item_id", "evalue"]).drop_duplicates("item_id")[
        ["item_id", "att_start", "att_end", "att_len", "evalue"]]


# --- score


def edge_errors(calls: pd.DataFrame, truth: pd.DataFrame) -> np.ndarray:
    """|edge errors| (bp) of calls overlapping exactly one element, both edges pooled."""
    out = []
    for c in calls.itertuples():
        t = truth[(truth.accession == c.genome_id) & (truth.contig == c.contig_id)
                  & (truth.end > c.start) & (truth.start < c.end)]
        if len(t) == 1:
            out += [abs(int(c.start) - int(t.start.iloc[0])), abs(int(c.end) - int(t.end.iloc[0]))]
    return np.asarray(out)


def _cost(ev: pd.DataFrame, reg: pd.DataFrame) -> dict:
    """bp called and bp called off every truth element."""
    bp, on = int((reg.end - reg.start).sum()), 0
    for c in reg.itertuples():
        on += sum(max(0, min(c.end, t.end) - max(c.start, t.start))
                  for t in ev[(ev.accession == c.genome_id) & (ev.contig == c.contig_id)].itertuples())
    return {"Mb_called": bp / 1e6, "Mb_off": (bp - on) / 1e6}


def metrics(ev: pd.DataFrame, reg: pd.DataFrame, ek) -> dict:
    v = ek.score_variant(ev, reg)
    err = edge_errors(reg, ev)
    return {"covered50": ek.covered50(ev, reg), "covered": v["covered_recall"], "matched": v["matched_recall"],
            "median_iou": v["median_iou"], "absorbed": v["absorbed_elements"], **_cost(ev, reg),
            "edge_med_bp": float(np.median(err)) if len(err) else np.nan,
            "edge_le100": float(np.mean(err <= 100)) if len(err) else np.nan,
            "edge_le1k": float(np.mean(err <= 1000)) if len(err) else np.nan}


def cmd_score(args) -> None:
    import evalkit as ek

    items = pd.read_csv(f"{args.run}.items.tsv", sep="\t")
    cands = pd.read_csv(f"{args.run}.cands.tsv", sep="\t")
    truth = pd.read_csv(args.elements)
    calls_all = items[items.kind == "call"]
    rows = []
    base_by_mode = {}
    for mode, cm in calls_all.groupby("mode", sort=False):
        ev = truth[truth.accession.isin(cm.genome_id.unique())]
        reg = cm.rename(columns={"call_start": "start", "call_end": "end"})[["genome_id", "contig_id", "start", "end"]]
        base_by_mode[mode] = (ev, cm)
        rows.append({"mode": mode, "rule": "shipped", "att_calls": np.nan, "att_decoys": np.nan, **metrics(ev, reg, ek)})
    grid = itertools.product(args.w_in, args.w_out, args.max_len, args.max_copies, args.core_slack, args.e_max)
    for w_in, w_out, max_len, max_copies, slack, e_max in grid:
        sl = None if slack < 0 else slack
        sel = select(items, cands, w_in, w_out, max_len, max_copies, args.min_entropy, sl, e_max)
        hit = items.merge(sel, on="item_id", how="left")
        rule = f"win{w_in}_wout{w_out}_len{max_len}_cp{max_copies}_slack{slack}_E{e_max:g}"
        for mode, (ev, cm) in base_by_mode.items():
            h = hit[hit["mode"] == mode]
            hc, hd = h[h.kind == "call"], h[(h.kind == "decoy") & h.searchable]
            new = hc.assign(start=np.where(hc.att_start.notna(), hc.att_start, hc.call_start).astype(int),
                            end=np.where(hc.att_end.notna(), hc.att_end, hc.call_end).astype(int))
            reg = new[["genome_id", "contig_id", "start", "end"]]
            att_only = new[new.att_start.notna()][["genome_id", "contig_id", "start", "end"]]
            err_att = edge_errors(att_only, ev)
            rows.append({"mode": mode, "rule": rule, "att_calls": float(hc.att_start.notna().mean()),
                         "att_decoys": float(hd.att_start.notna().mean()), **metrics(ev, reg, ek),
                         "att_n": int(len(att_only)),
                         "att_edge_med_bp": float(np.median(err_att)) if len(err_att) else np.nan,
                         "att_edge_le100": float(np.mean(err_att <= 100)) if len(err_att) else np.nan,
                         "att_edge_le1k": float(np.mean(err_att <= 1000)) if len(err_att) else np.nan})
        print(f"[att] {rule}", flush=True)
    out = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(f"{args.out}.csv", index=False)
    print(f"[att] wrote {args.out}.csv ({len(out)} rows)")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    r = sp.add_parser("run")
    r.add_argument("--tracks", required=True, help="run dir with <genome>/windows.tsv holding prob_ckpt1..3")
    r.add_argument("--genomes_dir", required=True)
    r.add_argument("--seed_tracks", nargs="+", help="per-checkpoint run dirs, when --tracks has no prob_ckpt<k>")
    r.add_argument("--elements", help="truth table, so decoys avoid true elements")
    r.add_argument("--modes", nargs="+", default=["exploratory", "default", "conservative"])
    r.add_argument("--n_decoys", type=int, default=5)
    r.add_argument("--seed", type=int, default=0)
    r.add_argument("--workers", type=int, default=8)
    r.add_argument("--out", required=True, help="output prefix")
    r.set_defaults(fn=cmd_run)
    s = sp.add_parser("score")
    s.add_argument("--run", required=True, help="prefix given to `run --out`")
    s.add_argument("--elements", required=True)
    s.add_argument("--w_in", type=int, nargs="+", default=[10000, 15000, 25000])
    s.add_argument("--w_out", type=int, nargs="+", default=[5000, 10000])
    s.add_argument("--max_len", type=int, nargs="+", default=[100, 200, 300, 100000])
    s.add_argument("--max_copies", type=int, nargs="+", default=[2])
    s.add_argument("--min_entropy", type=float, default=1.5)
    s.add_argument("--core_slack", type=int, nargs="+", default=[-1, 0, 2000],
                   help="bp the att interval may fall inside the high-scoring core; -1 = no core rule")
    s.add_argument("--e_max", type=float, nargs="+", default=[1e-2, 1e-4, 1e-6])
    s.add_argument("--out", required=True)
    s.set_defaults(fn=cmd_score)
    args = ap.parse_args()
    args.fn(args)


if __name__ == "__main__":
    main()
