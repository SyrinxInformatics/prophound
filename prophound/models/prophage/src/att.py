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
"""Call boundaries: att sites where the genome records them, ranges elsewhere.

docs/att-boundaries.md. Integration by site-specific recombination leaves a
prophage between two copies of the att core (attL, attR), a direct repeat
typically tens of bp long. When a call has such a pair, one copy near each edge,
its boundaries are known to the base:

  - search: every exact direct repeat (12-mer seeds extended along their
    diagonal, so exhaustive above 12 bp) between a zone around each edge of the
    call's coarse coordinates, searched 25 kb inward and 10 kb outward;
  - accept: both copies within 10 kb inward / 5 kb outward of the coarse edges,
    core <= 300 bp, at most 2 copies in the genome on either strand (dispersed
    repeats have more), base entropy >= 1.5 bits, the pair enclosing >= 2 kb,
    and E = m * n * (1 - p) * p**L <= 1e-6, the number of repeats that long
    expected by chance between the accepted zones (lengths m, n; p = sum of
    squared base frequencies). The smallest E wins; ties go to the leftmost pair.
  - coordinates: from the first base of attL to the last base of attR, so both
    core copies are inside the call.

The rule and its thresholds were chosen on the tuning sample and checked on two
others; the search zones are wider than the acceptance zones because that is how
it was measured, and a narrower search would cut some repeats differently.

Every other edge is a model estimate. Its range is relative to the call's length
L: the true edge lies from REL_OUT * L outside to REL_IN * L inside the reported
edge (the 5th and 95th percentiles of the tuning sample's model-edge errors; about
90% coverage on both check samples, but 0.78 on the held-out second lockbox, so these
are not called 90% ranges: docs/lockbox2.md "Result"). An att edge has a zero-width range, and so
does an edge on a contig end, which is the end of the sequence, not a boundary.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

SEED_K = 12              # seed length; every exact repeat >= SEED_K is found
MAX_SEED_COPIES = 8      # a seed this frequent within one zone is low complexity / tandem
MIN_REPEAT_LEN = 16      # shorter repeats are never significant in kb-sized zones
MIN_ELEMENT_BP = 2000    # the pair must enclose at least this much sequence
SEARCH_IN_BP, SEARCH_OUT_BP = 25000, 10000
ACCEPT_IN_BP, ACCEPT_OUT_BP = 10000, 5000
MAX_CORE_LEN = 300
MAX_COPIES = 2
MIN_ENTROPY_BITS = 1.5
MAX_EVALUE = 1e-6

REL_OUT = 0.055          # edge ranges, as fractions of the call's length
REL_IN = 0.494
# the same after composition edges, calibrated the same way on Default and Conservative model
# edges (0.898 on the tuning sample; 0.891 and 0.883 on the two check samples); used in every
# mode since Exploratory gained composition edges (0.89-0.91 there, LAMBDA 0.81)
RANGE_FRACS_COMPOSITION = (0.143, 0.463)

_CODE = np.full(256, 255, dtype=np.uint8)
for _i, _b in enumerate(b"ACGT"):
    _CODE[_b] = _i
    _CODE[_b + 32] = _i
_COMP = str.maketrans("ACGTacgt", "TGCAtgca")


class Genome:
    """Upper-case contig sequences, all of them joined (for copy counts), and the
    chance that two random bases match (from the genome's own composition)."""

    def __init__(self, seqs: dict[str, str]):
        self.seqs = {k: v.upper() for k, v in seqs.items()}
        self.text = "|".join(self.seqs.values())
        f = np.array([self.text.count(b) for b in "ACGT"], dtype=float)
        self.p_match = float(((f / f.sum()) ** 2).sum()) if f.sum() else 0.25

    @classmethod
    def from_fasta(cls, path: str) -> Genome:
        from pyfaidx import Fasta

        fa = Fasta(path, sequence_always_upper=True, rebuild=False)
        try:
            return cls({k: str(fa[k][:]) for k in fa.keys()})  # noqa: SIM118 -- iterating a Fasta yields records
        finally:
            fa.close()

    def copies(self, core: str) -> int:
        rc = core.translate(_COMP)[::-1]
        return self.text.count(core) + (self.text.count(rc) if rc != core else 0)


def _seeds(seq: str) -> pd.DataFrame:
    """Every valid SEED_K-mer in `seq` with its offset, minus the over-frequent ones."""
    a = _CODE[np.frombuffer(seq.encode(), dtype=np.uint8)]
    n = len(a) - SEED_K + 1
    if n <= 0:
        return pd.DataFrame({"code": np.zeros(0, np.uint64), "pos": np.zeros(0, np.int64)})
    bad = np.concatenate([[0], np.cumsum(a == 255)])
    ok = (bad[SEED_K:SEED_K + n] - bad[:n]) == 0
    code = np.zeros(n, dtype=np.uint64)
    a64 = np.where(a == 255, 0, a).astype(np.uint64)
    for i in range(SEED_K):
        code = code * np.uint64(4) + a64[i:i + n]
    t = pd.DataFrame({"code": code[ok], "pos": np.nonzero(ok)[0]})
    return t[t.code.map(t.code.value_counts()) <= MAX_SEED_COPIES]


def direct_repeats(left: str, right: str) -> pd.DataFrame:
    """Maximal exact direct repeats >= SEED_K between two strings: offsets li, ri and length."""
    hits = _seeds(left).merge(_seeds(right), on="code", suffixes=("_l", "_r"))
    if hits.empty:
        return pd.DataFrame({"li": [], "ri": [], "length": []}, dtype=np.int64)
    hits["diag"] = hits.pos_r - hits.pos_l
    hits = hits.sort_values(["diag", "pos_l"])
    hits["run"] = ((hits.diag.diff() != 0) | (hits.pos_l.diff() != 1)).cumsum()
    g = hits.groupby("run").agg(li=("pos_l", "min"), ri=("pos_r", "min"), n=("pos_l", "size"))
    return pd.DataFrame({"li": g.li.to_numpy(), "ri": g.ri.to_numpy(), "length": (g.n + SEED_K - 1).to_numpy()})


def entropy_bits(s: str) -> float:
    c = np.array([s.count(x) for x in "ACGT"], dtype=float)
    p = c[c > 0] / c.sum()
    return float(-(p * np.log2(p)).sum())


def _zones(s: int, e: int, contig_len: int, inward: int, outward: int):
    mid = (s + e) // 2
    return max(0, s - outward), min(mid, s + inward), max(mid, e - inward), min(contig_len, e + outward)


def find_att(genome: Genome, contig_id: str, s: int, e: int) -> dict | None:
    """The accepted att pair for a call whose coarse coordinates are [s, e), or None.

    Returns att_start, att_end (attL's first base to attR's end, 0-based half-open),
    att_len, att_core and evalue.
    """
    contig = genome.seqs[contig_id]
    l0, l1, r0, r1 = _zones(s, e, len(contig), SEARCH_IN_BP, SEARCH_OUT_BP)
    a0, a1, b0, b1 = _zones(s, e, len(contig), ACCEPT_IN_BP, ACCEPT_OUT_BP)
    if l1 - l0 < SEED_K or r1 - r0 < SEED_K:
        return None
    reps = direct_repeats(contig[l0:l1], contig[r0:r1])
    m, n = max(1, a1 - a0), max(1, b1 - b0)
    p = genome.p_match
    best = None
    for r in reps.itertuples(index=False):
        length = int(r.length)
        start, end = l0 + int(r.li), r0 + int(r.ri) + length
        if (length < MIN_REPEAT_LEN or length > MAX_CORE_LEN or end - start < MIN_ELEMENT_BP
                or start < a0 or start + length > a1 or end - length < b0 or end > b1):
            continue
        evalue = m * n * (1 - p) * p ** length
        if evalue > MAX_EVALUE:
            continue
        key = (evalue, start, end)
        if best is not None and key >= best[0]:  # pylint: disable=unsubscriptable-object
            continue
        core = contig[start:start + length]
        if entropy_bits(core) < MIN_ENTROPY_BITS or genome.copies(core) > MAX_COPIES:
            continue
        best = (key, {"att_start": start, "att_end": end, "att_len": length, "att_core": core, "evalue": evalue})
    return None if best is None else best[1]


def apply_boundaries(calls: pd.DataFrame, genome: Genome | None, contig_lengths: dict[str, int]) -> pd.DataFrame:
    """Att boundaries where found (when `genome` is given), and an edge range for every call.

    A call is searched when it is not a short-contig call, has no edge on a contig
    end and was not split by refinement; the search uses its coarse coordinates
    (`coarse_start`/`coarse_end`, or start/end without refinement). Adds `boundary`
    (att or model), `model_start`/`model_end` (the model's edges), `att_len`,
    `att_core` and `start_lo`/`start_hi`/`end_lo`/`end_hi`.
    """
    out = calls.copy()
    n = len(out)
    out["model_start"], out["model_end"] = out["start"].astype(int), out["end"].astype(int)
    boundary = np.array(["model"] * n, dtype=object)
    att_len = np.full(n, np.nan)
    att_core = np.array([""] * n, dtype=object)
    if genome is not None and n:
        cs = out["coarse_start"] if "coarse_start" in out else out["start"]
        ce = out["coarse_end"] if "coarse_end" in out else out["end"]
        pieces = out["n_pieces"] if "n_pieces" in out else pd.Series(1, index=out.index)
        if "short" in out:
            short = out["short"].astype(bool)
        elif "boundary_kind" in out:   # predict marks short-contig passes in boundary_kind
            short = out["boundary_kind"].astype(str).str.contains("_short")
        else:
            short = pd.Series(False, index=out.index)
        edge = out["contig_edge"] if "contig_edge" in out else pd.Series("none", index=out.index)
        starts, ends = out["start"].to_numpy().copy(), out["end"].to_numpy().copy()
        for i in range(n):
            if bool(short.iloc[i]) or edge.iloc[i] != "none" or pd.isna(cs.iloc[i]) or int(pieces.iloc[i]) != 1:
                continue
            hit = find_att(genome, out["contig_id"].iloc[i], int(cs.iloc[i]), int(ce.iloc[i]))
            if hit:
                starts[i], ends[i] = hit["att_start"], hit["att_end"]
                boundary[i], att_len[i], att_core[i] = "att", hit["att_len"], hit["att_core"]
        out["start"], out["end"] = starts.astype(int), ends.astype(int)
        if "region_len_bp" in out:
            out["region_len_bp"] = out["end"] - out["start"]
    out["boundary"] = boundary
    out["att_len"] = pd.array(att_len, dtype="Int64") if n else pd.array([], dtype="Int64")
    out["att_core"] = att_core
    out["start_source"], out["end_source"] = _sources(out, "start"), _sources(out, "end")
    return edge_ranges(out, contig_lengths)


def _sources(calls: pd.DataFrame, side: str) -> np.ndarray:
    """What placed each edge: att, contig_end, fine (moved by refinement) or coarse."""
    n = len(calls)
    edge = calls["contig_edge"].to_numpy() if "contig_edge" in calls else np.array(["none"] * n)
    att = (calls["boundary"] == "att").to_numpy()
    on_end = np.isin(edge, [side, "both"])
    moved = np.zeros(n, bool)
    if f"coarse_{side}" in calls and "refine_status" in calls:
        refined = calls["refine_status"].isin(["refined", "split"]).to_numpy()
        moved = refined & (calls[side].to_numpy() != calls[f"coarse_{side}"].to_numpy())
    return np.select([att, on_end, moved], ["att", "contig_end", "fine"], "coarse").astype(object)


def uses_composition(mode: str, setting: str = "auto") -> bool:
    """--composition-edges: auto is on in every mode (src/composition.py MODES)."""
    from . import composition

    return setting == "on" or (setting == "auto" and mode in composition.MODES)


def finish_boundaries(calls: pd.DataFrame, genome: Genome | None, contig_lengths: dict[str, int], mode: str,
                      att_search: bool = True, composition_edges: str = "auto") -> pd.DataFrame:
    """Every boundary stage after finalize, as predict runs it: att boundaries, then composition
    edges (src/composition.py; by default in every mode, with its mode's gate), then the edge ranges.
    `genome` is needed for either stage; without it only the ranges are added."""
    from . import composition

    out = apply_boundaries(calls, genome if att_search else None, contig_lengths)
    if genome is not None and uses_composition(mode, composition_edges):
        out = composition.composition_edges(out, genome.seqs, composition.min_jsd(mode))
        out = edge_ranges(out, contig_lengths, *RANGE_FRACS_COMPOSITION)
    return out


def edge_ranges(calls: pd.DataFrame, contig_lengths: dict[str, int], rel_out: float = REL_OUT,
                rel_in: float = REL_IN) -> pd.DataFrame:
    """start_lo..start_hi and end_lo..end_hi: where each true edge lies (for a model edge, fitted
    as a 90% range; 0.78 held out)."""
    out = calls.copy()
    s, e = out["start"].to_numpy(dtype=np.int64), out["end"].to_numpy(dtype=np.int64)
    span = e - s
    exact = (out["boundary"] == "att").to_numpy() if "boundary" in out else np.zeros(len(out), bool)
    edge = out["contig_edge"].to_numpy() if "contig_edge" in out else np.array(["none"] * len(out))
    lens = out["contig_id"].map(contig_lengths).to_numpy(dtype=np.int64) if len(out) else np.zeros(0, np.int64)
    out_bp, in_bp = np.round(rel_out * span).astype(np.int64), np.round(rel_in * span).astype(np.int64)
    fix_s = exact | np.isin(edge, ["start", "both"])
    fix_e = exact | np.isin(edge, ["end", "both"])
    out["start_lo"] = np.where(fix_s, s, np.maximum(0, s - out_bp))
    out["start_hi"] = np.where(fix_s, s, s + in_bp)
    out["end_lo"] = np.where(fix_e, e, e - in_bp)
    out["end_hi"] = np.where(fix_e, e, np.minimum(lens, e + out_bp))
    return out
