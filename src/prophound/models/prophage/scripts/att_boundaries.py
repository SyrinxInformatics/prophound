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
call marks both boundaries to the base. For each shipped call this searches a
zone around each edge for exact direct repeats (k-mer seeds extended along
their diagonal), scores each by the number of repeats that long expected by
chance between zones of that size and composition, and rejects low-complexity
and dispersed (multi-copy) repeats. The same search runs on decoy zones of the
same geometry in host sequence away from any call or truth element, which
measures the false-att rate directly instead of trusting the formula.

Subcommands:
  run    calls rebuilt as `predict` ships them (evalkit.shipped_calls) for each
         mode -> one row per call with the best repeat found and its decoys
  score  per-call table + truth -> edge errors (bp), recall and cost metrics for
         the shipped coordinates and for att-replaced ones, per threshold

Example:
  python scripts/att_boundaries.py run --tracks analysis/ngap_dose/impute__clean__0 \\
      --genomes_dir ../../../../arc_data/raw/genomes --out analysis/att/tune42_calls.tsv
  python scripts/att_boundaries.py score --calls analysis/att/tune42_calls.tsv \\
      --elements scripts/eval_sizestrat_elements.csv --out analysis/att/tune42
"""

from __future__ import annotations

import argparse
import json
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


def find_att(contig: str, genome_text: str, s: int, e: int, w_in: int, w_out: int, p_match: float,
             min_element: int = 2000, max_copies: int = 2, min_entropy: float = 1.5) -> dict | None:
    """Best att-like direct repeat for a call [s, e) on `contig` (the contig's sequence).

    The left zone is [s - w_out, s + w_in], the right [e - w_in, e + w_out], each
    stopping at the call's midpoint so the two never overlap. Returns the repeat with
    the smallest chance expectation (`evalue`), or None. `genome_text` is every contig
    joined, for counting the repeat's copies on both strands.
    """
    mid = (s + e) // 2
    l0, l1 = max(0, s - w_out), min(mid, s + w_in)
    r0, r1 = max(mid, e - w_in), min(len(contig), e + w_out)
    if l1 - l0 < SEED_K or r1 - r0 < SEED_K:
        return None
    reps = direct_repeats(contig[l0:l1], contig[r0:r1])
    if reps.empty:
        return None
    reps["start"] = l0 + reps.li
    reps["end"] = r0 + reps.ri + reps.length
    reps = reps[reps.end - reps.start >= min_element]
    m, n = l1 - l0, r1 - r0
    best = None
    for r in reps.sort_values("length", ascending=False).itertuples():
        evalue = m * n * (1 - p_match) * p_match ** r.length
        if best is not None and evalue >= best["evalue"]:
            break
        core = contig[r.start:r.start + r.length].upper()
        if entropy_bits(core) < min_entropy:
            continue
        rc = core.translate(COMP)[::-1]
        copies = genome_text.count(core) + (genome_text.count(rc) if rc != core else 0)
        if copies > max_copies:
            continue
        best = {"att_start": int(r.start), "att_end": int(r.end), "att_len": int(r.length),
                "att_core": core, "att_copies": int(copies), "evalue": float(evalue),
                "zone_left": m, "zone_right": n}
    return best


# --- run


def load_genome(path: Path) -> tuple[dict, str, float]:
    fa = pyfaidx.Fasta(str(path), sequence_always_upper=True, rebuild=False)
    seqs = {k: str(fa[k][:]) for k in fa.keys()}
    text = "|".join(seqs.values())
    f = np.array([text.count(x) for x in "ACGT"], dtype=float)
    f /= f.sum()
    return seqs, text, float((f ** 2).sum())


def decoy_offsets(rng, s: int, e: int, length: int, blocked: list, w_out: int, n: int) -> list:
    """Up to `n` shifts d for which [s+d - w_out, e+d + w_out] fits the contig and avoids `blocked`."""
    span0, span1 = s - w_out, e + w_out
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


def cmd_run(args) -> None:
    import evalkit as ek

    ens = ek.load_tracks(args.tracks)
    seeds = [ek.load_tracks(args.tracks, f"prob_ckpt{k}") for k in (1, 2, 3)]
    truth = pd.read_csv(args.elements) if args.elements else None
    rng = np.random.default_rng(args.seed)
    rows = []
    for mode in args.modes:
        reg = ek.shipped_calls(ens, seeds, mode, args.genomes_dir)
        for g, sub in reg.groupby("genome_id", sort=True):
            seqs, text, p = load_genome(Path(args.genomes_dir) / f"{g}.fna")
            for c in sub.itertuples():
                contig = seqs[c.contig_id]
                base = {"mode": mode, "genome_id": g, "contig_id": c.contig_id, "start": int(c.start),
                        "end": int(c.end), "contig_edge": getattr(c, "contig_edge", "none"),
                        "seed_support": getattr(c, "seed_support", np.nan), "p_match": p}
                hit = None
                if base["contig_edge"] in ("none", None) or pd.isna(base["contig_edge"]):
                    hit = find_att(contig, text, c.start, c.end, args.w_in, args.w_out, p,
                                   max_copies=args.max_copies)
                rows.append({**base, "kind": "call", **(hit or {})})
                # decoys: same geometry, shifted into host sequence away from calls and truth
                blocked = [(int(x.start) - args.w_out, int(x.end) + args.w_out)
                           for x in reg[(reg.genome_id == g) & (reg.contig_id == c.contig_id)].itertuples()]
                if truth is not None:
                    t = truth[(truth.accession == g) & (truth.contig == c.contig_id)]
                    blocked += [(int(x.start) - args.w_out, int(x.end) + args.w_out) for x in t.itertuples()]
                for d in decoy_offsets(rng, c.start, c.end, len(contig), blocked, args.w_out, args.n_decoys):
                    h = find_att(contig, text, c.start + d, c.end + d, args.w_in, args.w_out, p,
                                 max_copies=args.max_copies)
                    rows.append({**base, "kind": "decoy", "decoy_shift": d, **(h or {})})
        print(f"[att] {mode}: {len(reg)} calls", flush=True)
    out = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)
    print(f"[att] wrote {args.out} ({len(out)} rows)")


# --- score


def edge_errors(calls: pd.DataFrame, truth: pd.DataFrame, s_col: str, e_col: str) -> pd.DataFrame:
    """For each call overlapping exactly one element: signed edge errors (call - element)."""
    out = []
    for c in calls.itertuples():
        t = truth[(truth.accession == c.genome_id) & (truth.contig == c.contig_id)
                  & (truth.end > getattr(c, s_col)) & (truth.start < getattr(c, e_col))]
        if len(t) != 1:
            continue
        t = t.iloc[0]
        out.append({"mode": c.mode, "genome_id": c.genome_id, "start": c.start,
                    "left": int(getattr(c, s_col)) - int(t.start), "right": int(getattr(c, e_col)) - int(t.end),
                    "stratum": t.get("stratum", "")})
    return pd.DataFrame(out)


def summarize_errors(err: pd.DataFrame) -> dict:
    if err.empty:
        return {"n_calls": 0}
    a = np.abs(np.concatenate([err.left.to_numpy(), err.right.to_numpy()]))
    return {"n_calls": int(len(err)), "median_abs_bp": float(np.median(a)),
            "within_100": float(np.mean(a <= 100)), "within_500": float(np.mean(a <= 500)),
            "within_1000": float(np.mean(a <= 1000)), "within_5000": float(np.mean(a <= 5000))}


def cmd_score(args) -> None:
    import evalkit as ek

    df = pd.read_csv(args.calls, sep="\t")
    truth = pd.read_csv(args.elements)
    calls_df, dec = df[df.kind == "call"].copy(), df[df.kind == "decoy"]
    res = {}
    for mode, cm in calls_df.groupby("mode", sort=False):
        dm = dec[dec["mode"] == mode]
        n_g = cm.genome_id.nunique()
        ev = truth[truth.accession.isin(cm.genome_id.unique())] if not args.all_genomes else truth
        r = {"n_calls": int(len(cm)), "calls_per_genome": len(cm) / max(1, n_g), "by_threshold": {}}
        base = cm[["genome_id", "contig_id", "start", "end"]]
        r["shipped"] = {"covered50": ek.covered50(ev, base), **_cost(ev, base),
                        **{k: v for k, v in ek.score_variant(ev, base).items() if not k.startswith("_")},
                        "edges_all": summarize_errors(edge_errors(cm, ev, "start", "end"))}
        for thr in args.evalue:
            hit = cm.evalue.notna() & (cm.evalue <= thr)
            dhit = dm.evalue.notna() & (dm.evalue <= thr)
            new = cm.copy()
            new.loc[hit, "start"] = new.loc[hit, "att_start"].astype(int)
            new.loc[hit, "end"] = new.loc[hit, "att_end"].astype(int)
            nb = new[["genome_id", "contig_id", "start", "end"]]
            sub_ship, sub_att = cm[hit], new[hit]
            r["by_threshold"][str(thr)] = {
                "att_rate_calls": float(hit.mean()), "att_rate_decoys": float(dhit.mean()) if len(dm) else None,
                "n_decoys": int(len(dm)), "n_att": int(hit.sum()),
                "covered50": ek.covered50(ev, nb), **_cost(ev, nb),
                **{k: v for k, v in ek.score_variant(ev, nb).items() if not k.startswith("_")},
                "edges_att_calls_shipped": summarize_errors(edge_errors(sub_ship, ev, "start", "end")),
                "edges_att_calls_att": summarize_errors(edge_errors(sub_att, ev, "start", "end")),
                "edges_all": summarize_errors(edge_errors(new, ev, "start", "end")),
            }
        res[mode] = r
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(f"{args.out}.json").write_text(json.dumps(res, indent=2))
    print(json.dumps(res, indent=1)[:args.print_chars])


def _cost(ev: pd.DataFrame, reg: pd.DataFrame) -> dict:
    """bp called and bp called off every truth element."""
    bp = int((reg.end - reg.start).sum())
    on = 0
    for c in reg.itertuples():
        on += sum(max(0, min(c.end, t.end) - max(c.start, t.start))
                  for t in ev[(ev.accession == c.genome_id) & (ev.contig == c.contig_id)].itertuples())
    return {"bp_called": bp, "bp_off_element": bp - on}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    r = sp.add_parser("run")
    r.add_argument("--tracks", required=True, help="run dir with <genome>/windows.tsv holding prob_ckpt1..3")
    r.add_argument("--genomes_dir", required=True)
    r.add_argument("--elements", help="truth table, so decoys avoid true elements")
    r.add_argument("--modes", nargs="+", default=["exploratory", "default", "conservative"])
    r.add_argument("--w_in", type=int, default=20000)
    r.add_argument("--w_out", type=int, default=5000)
    r.add_argument("--max_copies", type=int, default=2)
    r.add_argument("--n_decoys", type=int, default=3)
    r.add_argument("--seed", type=int, default=0)
    r.add_argument("--out", required=True)
    r.set_defaults(fn=cmd_run)
    s = sp.add_parser("score")
    s.add_argument("--calls", required=True)
    s.add_argument("--elements", required=True)
    s.add_argument("--evalue", type=float, nargs="+", default=[1e-1, 1e-2, 1e-3, 1e-4, 1e-6])
    s.add_argument("--all_genomes", action="store_true")
    s.add_argument("--print_chars", type=int, default=6000)
    s.add_argument("--out", required=True)
    s.set_defaults(fn=cmd_score)
    args = ap.parse_args()
    args.fn(args)


if __name__ == "__main__":
    main()
