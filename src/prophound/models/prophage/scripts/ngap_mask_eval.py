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
"""Offline test of masking long N runs (docs/v1-scope.md A.3, user decision 2026-09-29).

822054 showed that `impute` turns long N gaps into random sequence the model then
calls. Proposed fix, simulated here on 822054's tracks and gap files:
  mask  a window whose gap content exceeds --frac of its length scores 0, in the
        ensemble and in every seed's track, so it can neither seed nor extend a
        call, and a masked stretch cannot be bridged by merging
  clip  gap intervals are subtracted from the calls afterwards, dropping pieces
        under --min_piece bp, so no call includes N sequence
Per condition, frac and shipped mode: covered50, elements lost against the same
setting on clean input, calls per genome, and calls lying >=50% inside random
gaps before clipping (the calls the fix exists to remove).

Example:
  python scripts/ngap_mask_eval.py --elements scripts/eval_sizestrat_elements.csv \\
      --root analysis/ngap_dose --out analysis/ngap_dose/mask_eval.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import evalkit as ek
import numpy as np
import pandas as pd
from ambiguity_dose import gap_call_stats
from evalkit import cov_fracs


def gap_bp_per_window(w: pd.DataFrame, gaps: pd.DataFrame) -> np.ndarray:
    out = np.zeros(len(w))
    for c, gg in gaps.groupby("contig_id"):
        m = (w.contig_id == c).to_numpy()
        if not m.any():
            continue
        s, e = w.start.to_numpy()[m], w.end.to_numpy()[m]
        tot = np.zeros(m.sum())
        for a, b in zip(gg.start, gg.end):
            tot += np.clip(np.minimum(e, b) - np.maximum(s, a), 0, None)
        out[m] = tot
    return out


def mask_tracks(tracks: dict, gaps: pd.DataFrame, frac: float) -> dict:
    if frac is None or gaps is None:
        return tracks
    out = {}
    for g, (w, wl) in tracks.items():
        gb = gap_bp_per_window(w, gaps[gaps.genome_id == g])
        w = w.copy()
        w.loc[gb > frac * (w.end - w.start).to_numpy(), "prob"] = 0.0
        out[g] = (w, wl)
    return out


def clip(reg: pd.DataFrame, gaps: pd.DataFrame, min_piece: int) -> pd.DataFrame:
    if gaps is None or reg.empty:
        return reg
    rows = []
    for r in reg.itertuples(index=False):
        gg = gaps[(gaps.genome_id == r.genome_id) & (gaps.contig_id == r.contig_id)
                  & (gaps.end > r.start) & (gaps.start < r.end)].sort_values("start")
        pieces, cur = [], int(r.start)
        for a, b in zip(gg.start, gg.end):
            if a > cur:
                pieces.append((cur, int(a)))
            cur = max(cur, int(b))
        if cur < int(r.end):
            pieces.append((cur, int(r.end)))
        for s, e in pieces:
            if e - s >= min_piece:
                rows.append({**r._asdict(), "start": s, "end": e})
    return pd.DataFrame(rows, columns=reg.columns)


def drop_gap_calls(reg: pd.DataFrame, gaps: pd.DataFrame, max_gap_frac: float = 0.5) -> pd.DataFrame:
    """Drop calls with >= max_gap_frac of their length in gaps; keep the rest unchanged."""
    if gaps is None or reg.empty:
        return reg
    keep = []
    for r in reg.itertuples(index=False):
        gg = gaps[(gaps.genome_id == r.genome_id) & (gaps.contig_id == r.contig_id)
                  & (gaps.end > r.start) & (gaps.start < r.end)]
        iv = sorted((max(int(a), int(r.start)), min(int(b), int(r.end))) for a, b in zip(gg.start, gg.end))
        bp, cs, ce = 0, None, None
        for a, b in iv:                        # union: element and random gaps can overlap
            if ce is None or a > ce:
                bp += 0 if ce is None else ce - cs
                cs, ce = a, b
            else:
                ce = max(ce, b)
        bp += 0 if ce is None else ce - cs
        keep.append(bp < max_gap_frac * (int(r.end) - int(r.start)))
    return reg[np.array(keep, bool)].reset_index(drop=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--elements", required=True)
    ap.add_argument("--root", required=True, help="822054's analysis/ngap_dose (impute__<pattern>__<rate>/ + gaps/)")
    ap.add_argument("--fracs", default="none,0.1,0.25,0.5")
    ap.add_argument("--min_piece", type=int, default=1000)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    ev = pd.read_csv(a.elements)
    root = Path(a.root)
    conds = sorted(d.name for d in root.iterdir() if d.is_dir() and d.name.startswith("impute__"))
    # a frac is a window-masking threshold, or one of the no-mask variants:
    # "clip" (subtract gaps from calls) and "drop" (drop calls that are >=50% gap)
    fracs = [None if f == "none" else (f if f in ("clip", "drop") else float(f)) for f in a.fracs.split(",")]
    res, base = {}, {}
    for cond in ["impute__clean__0"] + [c for c in conds if c != "impute__clean__0"]:
        d = root / cond
        ens = ek.load_tracks(str(d))
        seeds = [ek.load_tracks(str(d), f"prob_ckpt{k}") for k in (1, 2, 3)]
        _, pat, rate = cond.split("__")
        gaps = None
        if pat == "ngap":
            gaps = pd.concat([pd.read_csv(f, sep="\t") for f in (root / "gaps" / f"{pat}__{rate}").glob("*.tsv")],
                             ignore_index=True)
        for frac in fracs:
            mf = frac if isinstance(frac, float) else None
            e2, s2 = mask_tracks(ens, gaps, mf), [mask_tracks(s, gaps, mf) for s in seeds]
            for mode in ek.MODE_THRESHOLDS:
                raw = ek.vote(ek.calls(e2, mode), [ek.calls(s, mode) for s in s2], ek.MODE_MIN_SEED_SUPPORT[mode])
                if frac == "drop":
                    reg = drop_gap_calls(raw, gaps)
                elif frac is None:
                    reg = raw
                else:                                   # masking thresholds and "clip"
                    reg = clip(raw, gaps, a.min_piece)
                c50 = cov_fracs(ev, reg) >= 0.5
                key = f"{cond}/{'none' if frac is None else frac}/{mode}"
                m = {"covered50": float(c50.mean()), "calls_per_genome": len(reg) / len(ens)}
                if gaps is not None:
                    m["calls_in_random_gaps_per_genome_before_clip"] = gap_call_stats(
                        raw, gaps[gaps.kind == "random"], len(ens))["calls_in_gaps_per_genome"]
                if pat == "clean":
                    base[(frac, mode)] = c50
                b = base.get((frac, mode), base.get((None, mode)))
                m["lost_vs_clean"] = int((b & ~c50).sum()) if b is not None else None
                res[key] = m
                print(f"[{key}] covered50 {m['covered50']:.3f} lost {m['lost_vs_clean']} calls/genome "
                      f"{m['calls_per_genome']:.1f} in-gap {m.get('calls_in_random_gaps_per_genome_before_clip', '-')}",
                      flush=True)
    Path(a.out).write_text(json.dumps(res, indent=2) + "\n")
    print(f"[done] wrote {a.out}")


if __name__ == "__main__":
    main()
