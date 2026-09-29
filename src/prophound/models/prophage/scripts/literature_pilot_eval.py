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
"""Literature-curated pilot (docs/final-benchmark.md, release checklist §3).

Scores calls against scripts/literature_truth.tsv, which is independent of
geNomad. Literature lists are incomplete (geNomad finds elements they omit), so
only recall-type quantities are reported, never precision.

Per run (seed or ensemble) and mode, calls are rebuilt from the score tracks as
the pipeline builds them (scripts/mode_recheck.py). Per element: covered
fraction, covered (>=1 bp), covered50, best single-call IoU. Reported by split
(test / val, never pooled), by class, and per genome.

Permutation test per genome (docs/final-benchmark.md): the calls on each contig
are rotated by one uniform random offset (circularly; this keeps their number,
lengths and spacing) N times. Statistics: the number of covered50 elements and
the mean best-match IoU. p = (1 + #null >= observed) / (1 + N); genomes are
combined with Fisher's method, separately for test and val.

geNomad's calls (find_proviruses) are scored the same way for comparison.

Example:
  python scripts/literature_pilot_eval.py --truth scripts/literature_truth.tsv \\
      --genomes_dir ../raw/genomes --genomad_dir ../genomad_out \\
      --run ens=analysis/litpilot/ens --out analysis/litpilot/results.json
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mode_recheck import calls, load_tracks, shipped_calls  # noqa: E402

try:  # repo layout
    from prophound.models.prophage.src.predict_genome import MODE_THRESHOLDS
except ImportError:  # HPCC data_refactored/ layout
    from src.predict_genome import MODE_THRESHOLDS


def read_fasta(path: Path) -> dict:
    seqs, name = {}, None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            else:
                seqs[name].append(line.strip().upper())
    return {k: "".join(v) for k, v in seqs.items()}


def verify_and_lengths(truth: pd.DataFrame, genomes_dir: Path) -> dict:
    """contig lengths per genome; exits if any boundary flank disagrees with the FASTA."""
    lengths, bad = {}, []
    for g, tg in truth.groupby("genome_id"):
        seqs = read_fasta(genomes_dir / f"{g}.fna")
        lengths[g] = {c: len(s) for c, s in seqs.items()}
        for r in tg.itertuples():
            s = seqs.get(r.contig, "")
            if s[r.start - 1:r.start + 29] != r.left_flank or s[r.end - 30:r.end] != r.right_flank:
                bad.append(f"{g} {r.name}")
    if bad:
        raise SystemExit(f"FASTA does not match the curated coordinates for: {bad}")
    return lengths


def overlaps(cs, ce, es, ee):
    """(n_calls, n_elements) overlap bp; everything 0-based half-open."""
    return np.clip(np.minimum(ce[:, None], ee[None, :]) - np.maximum(cs[:, None], es[None, :]), 0, None)


def element_scores(reg: pd.DataFrame, t: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for r in t.to_dict("records"):
        x = reg[(reg.genome_id == r["genome_id"]) & (reg.contig_id == r["contig"])]
        L = r["e0"] - r["s0"]
        if len(x):
            cs, ce = x.start.to_numpy(np.int64), x.end.to_numpy(np.int64)
            ov = overlaps(cs, ce, np.array([r["s0"]]), np.array([r["e0"]]))[:, 0]
            union = (ce - cs) + L - ov
            cov, iou = min(1.0, ov.sum() / L), float((ov / union).max())
        else:
            cov, iou = 0.0, 0.0
        rows.append({"genome_id": r["genome_id"], "name": r["name"], "class": r["class"], "split": r["split"],
                     "len_bp": L, "cov_frac": cov, "best_iou": iou})
    df = pd.DataFrame(rows)
    df["covered"], df["covered50"] = df.cov_frac > 0, df.cov_frac >= 0.5
    return df


def permutation(reg: pd.DataFrame, tg: pd.DataFrame, lengths: dict, n: int, rng) -> dict:
    """One genome: observed and null (covered50 count, mean best IoU) under per-contig rotation."""
    obs_c50, obs_iou = 0, 0.0
    null_c50, null_iou = np.zeros(n), np.zeros(n)
    for c, te in tg.groupby("contig"):
        L = lengths[c]
        es, ee = te.s0.to_numpy(np.int64), te.e0.to_numpy(np.int64)
        elen = ee - es
        x = reg[reg.contig_id == c]
        cs, ce = x.start.to_numpy(np.int64), x.end.to_numpy(np.int64)
        if len(x) == 0:
            continue

        def score(cs_, ce_):
            ov = overlaps(cs_, ce_, es, ee)
            ov = ov + overlaps(cs_ - L, ce_ - L, es, ee)   # the part of a rotated call that wrapped
            cov = np.minimum(1.0, ov.sum(0) / elen)
            iou = (ov / ((ce_ - cs_)[:, None] + elen[None, :] - ov)).max(0)
            return int((cov >= 0.5).sum()), float(iou.sum())

        o = score(cs, ce)
        obs_c50 += o[0]
        obs_iou += o[1]
        offs = rng.integers(0, L, size=n)
        for i, off in enumerate(offs):
            s_ = (cs + off) % L
            e_ = s_ + (ce - cs)
            r_ = score(s_, e_)
            null_c50[i] += r_[0]
            null_iou[i] += r_[1]
    m = len(tg)
    return {"n_elements": m, "obs_covered50": obs_c50, "obs_mean_iou": obs_iou / m,
            "null_covered50_mean": float(null_c50.mean()), "null_mean_iou_mean": float(null_iou.mean() / m),
            "p_covered50": float((1 + (null_c50 >= obs_c50).sum()) / (1 + n)),
            "p_mean_iou": float((1 + (null_iou >= obs_iou).sum()) / (1 + n))}


def fisher(ps: list) -> float:
    """Fisher's combined p (chi-square with 2k df, closed form for even df)."""
    if not ps:
        return float("nan")
    x = -2 * sum(math.log(p) for p in ps)
    k = len(ps)
    return float(math.exp(-x / 2) * sum((x / 2) ** i / math.factorial(i) for i in range(k)))


def summarise(es: pd.DataFrame) -> dict:
    out = {}
    for split, sub in es.groupby("split"):
        d = {"n": int(len(sub)), "covered": float(sub.covered.mean()), "covered50": float(sub.covered50.mean()),
             "median_best_iou": float(sub.best_iou.median()),
             "by_class": {c: {"n": int(len(x)), "covered50": float(x.covered50.mean())} for c, x in sub.groupby("class")},
             "named_prophages_covered50": float(sub[sub["class"] != "remnant"].covered50.mean()),
             "missed": sub[~sub.covered50][["genome_id", "name", "class", "len_bp", "cov_frac"]].to_dict("records")}
        out[split] = d
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--genomes_dir", required=True)
    ap.add_argument("--genomad_dir", default=None)
    ap.add_argument("--run", action="append", required=True, metavar="LABEL=DIR")
    ap.add_argument("--n_perm", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=20260928)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    truth = pd.read_csv(args.truth, sep="\t")
    lengths = verify_and_lengths(truth, Path(args.genomes_dir))
    truth["s0"], truth["e0"] = truth.start - 1, truth.end   # 1-based inclusive -> 0-based half-open
    print(f"[truth] {len(truth)} elements; flanks verified against {args.genomes_dir}", flush=True)
    rng = np.random.default_rng(args.seed)
    res = {"truth": args.truth, "n_perm": args.n_perm, "runs": {}}
    per_element = []

    def evaluate(label: str, reg: pd.DataFrame) -> dict:
        es = element_scores(reg, truth)
        per_element.append(es.assign(run=label))
        d = summarise(es)
        perm = {}
        for g, tg in truth.groupby("genome_id"):
            perm[g] = permutation(reg[reg.genome_id == g], tg, lengths[g], args.n_perm, rng)
        for split in d:
            gs = truth[truth.split == split].genome_id.unique()
            d[split]["fisher_p_covered50"] = fisher([perm[g]["p_covered50"] for g in gs])
            d[split]["fisher_p_mean_iou"] = fisher([perm[g]["p_mean_iou"] for g in gs])
            d[split]["calls_per_genome"] = float(reg[reg.genome_id.isin(gs)].groupby("genome_id").size()
                                                 .reindex(gs, fill_value=0).mean())
        d["permutation"] = perm
        print(f"[{label}] " + " | ".join(f"{s}: covered50 {d[s]['covered50']:.3f} (n={d[s]['n']}), "
                                         f"prophages {d[s]['named_prophages_covered50']:.3f}, "
                                         f"Fisher p {d[s]['fisher_p_covered50']:.2g}/{d[s]['fisher_p_mean_iou']:.2g}"
                                         for s in ("test", "val") if s in d), flush=True)
        return d

    loaded = {}
    for spec in args.run:
        label, _, d = spec.partition("=")
        tracks = load_tracks(d)
        missing = set(truth.genome_id) - set(tracks)
        if missing:
            raise SystemExit(f"run {label} lacks genomes {sorted(missing)}")
        loaded[label] = tracks = {g: tracks[g] for g in truth.genome_id.unique()}
        for mode in MODE_THRESHOLDS:
            res["runs"][f"{label}/{mode}"] = evaluate(f"{label}/{mode}", calls(tracks, mode))
    if {"ens", "r1", "r2", "r3"} <= set(loaded):   # the shipped pipeline: each mode's seed vote
        seeds = [loaded[s] for s in ("r1", "r2", "r3")]
        for mode in MODE_THRESHOLDS:
            res["runs"][f"ens_shipped/{mode}"] = evaluate(f"ens_shipped/{mode}",
                                                          shipped_calls(loaded["ens"], seeds, mode))

    if args.genomad_dir:
        rows = []
        for g in truth.genome_id.unique():
            p = Path(args.genomad_dir) / g / f"{g}_find_proviruses" / f"{g}_provirus.tsv"
            if p.exists() and p.stat().st_size > 10:
                t = pd.read_csv(p, sep="\t")
                rows += [{"genome_id": g, "contig_id": r.source_seq, "start": int(r.start) - 1, "end": int(r.end)}
                         for r in t.itertuples()]
        res["runs"]["genomad"] = evaluate("genomad", pd.DataFrame(rows, columns=["genome_id", "contig_id", "start", "end"]))

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(res, indent=2, default=str) + "\n")
    pd.concat(per_element).to_csv(out.with_suffix(".elements.csv"), index=False)
    print(f"[done] wrote {out}", flush=True)


if __name__ == "__main__":
    main()
