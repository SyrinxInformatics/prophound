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
"""Composition edges: move an estimated edge to where the sequence turns from host-like
to element-like (docs/att-boundaries.md "Composition edges").

Prophages usually differ from their host in base composition; Alien_Hunter's change-point
step and CheckV's adjacent-window comparison use the same idea. For each `model` edge not
on a contig end, of a call at least MIN_CALL_BP long:

  - profiles: tetranucleotide frequencies of the element (the call's middle 60%) and of
    the local host (HOST_FLANK bp beyond the edge's range, on the outer side)
  - windows: WIN bp every STEP bp across the edge's range, each scored by the
    log-likelihood ratio element vs host (multinomial, pseudocount PSEUDO)
  - change-point: the split that maximises (sum of scores on the element side) - (sum on
    the host side); its first element-side window centre is the new edge
  - gates: the profiles' Jensen-Shannon divergence >= MIN_JSD, and an inward move of at
    most min(MAX_INWARD_FRAC x the call's length, MAX_INWARD_BP). A partial detection's or a
    multi-element call's middle is not the element, and a large move then cuts into it: the
    absolute ceiling stopped a 22.8 kb move that cut a 13 kb LAMBDA element off a 156 kb call

Chosen on the tuning sample and checked on two others; shipped in Default and
Conservative mode (user decision 2026-10-21), where it was recall-neutral (Default, all
three samples) and roughly halved the median edge error.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

K = 4
WIN, STEP = 1000, 100
HOST_FLANK = 20000
MIN_HOST_BP = 5000
PSEUDO = 0.5
MIN_CALL_BP = 3000
MIN_JSD = 0.01
MAX_INWARD_FRAC = 0.15
MAX_INWARD_BP = 8000
MODES = ("default", "conservative")

_CODE = np.full(256, 255, dtype=np.uint8)
for _i, _b in enumerate(b"ACGT"):
    _CODE[_b] = _i
    _CODE[_b + 32] = _i


def kmer_codes(seq: str) -> np.ndarray:
    """Tetranucleotide code per position (-1 where the k-mer has a non-ACGT base)."""
    a = _CODE[np.frombuffer(seq.encode(), dtype=np.uint8)]
    n = len(a) - K + 1
    if n <= 0:
        return np.zeros(0, np.int64)
    bad = np.concatenate([[0], np.cumsum(a == 255)])
    ok = (bad[K:K + n] - bad[:n]) == 0
    code = np.zeros(n, np.int64)
    a64 = np.where(a == 255, 0, a).astype(np.int64)
    for i in range(K):
        code = code * 4 + a64[i:i + n]
    return np.where(ok, code, -1)


def profile(codes: np.ndarray) -> np.ndarray:
    c = np.bincount(codes[codes >= 0], minlength=4 ** K).astype(float) + PSEUDO
    return c / c.sum()


def jsd(p: np.ndarray, q: np.ndarray) -> float:
    m = (p + q) / 2
    return float(0.5 * (p * np.log2(p / m)).sum() + 0.5 * (q * np.log2(q / m)).sum())


def change_point(contig: str, lo: int, hi: int, elem: np.ndarray, host: np.ndarray, side: str,
                 edge: int) -> int:
    """The best host -> element switch in [lo, hi] (`edge` if there is too little sequence)."""
    a, b = max(0, lo - WIN // 2), min(len(contig), hi + WIN // 2)
    codes = kmer_codes(contig[a:b])
    if len(codes) < WIN:
        return edge
    w = np.log(elem) - np.log(host)
    per = np.where(codes >= 0, w[np.maximum(codes, 0)], 0.0)
    c = np.concatenate([[0.0], np.cumsum(per)])
    starts = np.arange(0, len(per) - WIN + 1, STEP)
    llr = c[starts + WIN] - c[starts]
    centres = a + starts + WIN // 2
    keep = (centres >= lo) & (centres <= hi)
    llr, centres = llr[keep], centres[keep]
    if len(llr) < 3:
        return edge
    sign = 1.0 if side == "start" else -1.0     # the element lies right of a start, left of an end
    cs = np.concatenate([[0.0], np.cumsum(sign * llr)])
    j = int(np.argmax(cs[-1] - 2 * cs))
    return int(centres[j]) if j < len(centres) else int(centres[-1])


def composition_edges(calls: pd.DataFrame, seqs: dict) -> pd.DataFrame:
    """Move the `model` edges of `calls` (with their ranges) to composition change-points.

    Returns the calls with start/end moved where the gates pass, and `start_source` /
    `end_source` set to `composition` for moved edges.
    """
    out = calls.copy()
    if out.empty:
        return out
    starts, ends = out["start"].to_numpy().copy(), out["end"].to_numpy().copy()
    src_s, src_e = out["start_source"].to_numpy().copy(), out["end_source"].to_numpy().copy()
    for i, r in enumerate(out.itertuples(index=False)):
        s, e = int(r.start), int(r.end)
        L = e - s
        if r.boundary != "model" or L < MIN_CALL_BP:
            continue
        contig = seqs[r.contig_id]
        elem = profile(kmer_codes(contig[s + int(0.2 * L):e - int(0.2 * L)]))
        for side in ("start", "end"):
            if r.contig_edge in (side, "both"):
                continue
            lo, hi = int(getattr(r, f"{side}_lo")), int(getattr(r, f"{side}_hi"))
            host_seq = contig[max(0, lo - HOST_FLANK):lo] if side == "start" else contig[hi:hi + HOST_FLANK]
            if len(host_seq) < MIN_HOST_BP:
                continue
            host = profile(kmer_codes(host_seq))
            if jsd(elem, host) < MIN_JSD:
                continue
            old = s if side == "start" else e
            new = change_point(contig, lo, hi, elem, host, side, old)
            inward = (new - old) if side == "start" else (old - new)
            cap = min(MAX_INWARD_FRAC * L, MAX_INWARD_BP)
            if inward > cap:
                new = old + cap if side == "start" else old - cap
            new = int(np.round(new))
            if side == "start":
                starts[i], src_s[i] = new, ("composition" if new != old else src_s[i])
            else:
                ends[i], src_e[i] = new, ("composition" if new != old else src_e[i])
    out["start"], out["end"] = starts, ends
    out["start_source"], out["end_source"] = src_s, src_e
    if "region_len_bp" in out:
        out["region_len_bp"] = out["end"] - out["start"]
    return out
