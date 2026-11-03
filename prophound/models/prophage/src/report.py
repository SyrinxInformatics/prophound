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
"""What a researcher reads first: the ranked, plain calls table and the terminal list.

docs/interpreting-results.md explains them; docs/output-files.md is the contract.
`regions.tsv` is this simple table, in 1-based inclusive coordinates as GenBank, GFF3
and genome browsers show them. `regions.detailed.tsv` keeps every column, 0-based
half-open like BED.

The rank orders calls by `max_prob`, ties broken by `mean_prob`: on the tuning
sample this put calls on known elements first better than `mean_prob` alone, in
every mode, and it held on the check sample (docs/att-boundaries.md "Ranking").
Ranking by resemblance to known prophages means a novel candidate can rank low.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

GUIDE_URL = "https://github.com/SyrinxInformatics/mge-detect/blob/main/docs/interpreting-results.md"
SIMPLE_COLUMNS = ["rank", "call_id", "contig", "start", "end", "length_kb", "edges", "start_range", "end_range",
                  "support", "note"]
N_FRAC_NOTE = 0.1    # mention N content from this fraction up


def rank_calls(calls: pd.DataFrame) -> pd.Series:
    """Rank 1..N within each genome: max_prob (as written, 4 dp) then mean_prob, both high first."""
    if calls.empty:
        return pd.Series([], dtype="int64", index=calls.index)
    key = calls.assign(_mx=calls["max_prob"].astype(float).round(4), _mn=calls["mean_prob"].astype(float),
                       _st=calls["start"].astype(int))
    order = key.sort_values(["_mx", "_mn", "_st"], ascending=[False, False, True], kind="stable")
    group = order["genome_id"] if "genome_id" in order else pd.Series(0, index=order.index)
    return (order.groupby(group, sort=False).cumcount() + 1).reindex(calls.index).astype("int64")


def _range(lo: int, hi: int) -> str:
    return str(lo) if lo == hi else f"{lo}-{hi}"


def _note(r) -> str:
    notes = []
    if getattr(r, "boundary", "model") == "att":
        notes.append(f"att site: {int(r.att_len)} bp repeat at both edges")
    edge = getattr(r, "contig_edge", "none")
    if getattr(r, "boundary", "model") == "model" and getattr(r, "refine_status", "") in ("refined", "split"):
        # an edge neither the fine pass nor composition moved off the coarse 1 kb edge
        stuck = [side for side in ("start", "end") if getattr(r, f"{side}_source", "") == "coarse"]
        if stuck:
            what = "start and end" if len(stuck) == 2 else stuck[0]
            notes.append(f"{what} not refined: the score stays high past {'them' if len(stuck) == 2 else 'it'}, "
                         "so the true edge is less certain")
    if edge in ("start", "end", "both"):
        where = {"start": "start", "end": "end", "both": "both ends"}[edge]
        notes.append(f"runs off the contig {where}; the element may continue")
    if float(getattr(r, "n_frac", 0) or 0) >= N_FRAC_NOTE:
        notes.append(f"{100 * float(r.n_frac):.0f}% N")
    if getattr(r, "refine_status", "") == "short_contig":
        notes.append("on a short contig")
    return "; ".join(notes)


def simple_table(detailed: pd.DataFrame, n_checkpoints: int, with_genome: bool = False) -> pd.DataFrame:
    """The researcher's table from the detailed one (0-based) -> 1-based inclusive, sorted by rank."""
    cols = (["genome"] if with_genome else []) + SIMPLE_COLUMNS
    if detailed.empty:
        return pd.DataFrame(columns=cols)
    d = detailed
    s, e = d["start"].astype(int), d["end"].astype(int)
    att = (d["boundary"] == "att") if "boundary" in d else pd.Series(False, index=d.index)
    lo_s = d["start_lo"].astype(int) if "start_lo" in d else s
    hi_s = d["start_hi"].astype(int) if "start_hi" in d else s
    lo_e = d["end_lo"].astype(int) if "end_lo" in d else e
    hi_e = d["end_hi"].astype(int) if "end_hi" in d else e
    support = (d["seed_support"].astype(int).astype(str) + f"/{n_checkpoints}") if "seed_support" in d \
        else pd.Series("1/1", index=d.index)
    out = pd.DataFrame({
        "rank": d["rank"] if "rank" in d else rank_calls(d),
        "call_id": d["call_id"],
        "contig": d["contig_id"],
        "start": s + 1,                    # 1-based first base
        "end": e,                          # 1-based last base (inclusive)
        "length_kb": np.round((e - s) / 1000, 1),
        "edges": np.where(att, "att site", "estimated"),
        "start_range": [_range(a + 1, b + 1) for a, b in zip(lo_s, hi_s)],
        "end_range": [_range(a, b) for a, b in zip(lo_e, hi_e)],
        "support": support,
        "note": [_note(r) for r in d.itertuples(index=False)],
    })
    if with_genome:
        out.insert(0, "genome", d["genome_id"].to_numpy())
    keys = (["genome"] if with_genome else []) + ["rank"]
    return out.sort_values(keys, kind="stable").reset_index(drop=True)[cols]


def terminal_lines(simple: pd.DataFrame, mode: str, limit: int = 10) -> list:
    """The ranked calls as a few plain lines for the terminal."""
    n = len(simple)
    if n == 0:
        return [f"No candidate prophage regions in {mode} mode.",
                "Exploratory mode (--mode exploratory) makes more calls, at more cost to follow up."]
    rows = [(str(r.rank), str(r.contig), f"{r.start:,}-{r.end:,}", f"{r.length_kb:.1f} kb",
             "att site" if r.edges == "att site" else "estimated",
             f"{r.support} models") for r in simple.head(limit).itertuples(index=False)]
    head = ("rank", "contig", "start-end (1-based)", "length", "edges", "support")
    w = [max(len(h), *(len(x[i]) for x in rows)) for i, h in enumerate(head)]
    fmt = "  ".join(f"{{:<{x}}}" for x in w)
    lines = [f"{n} candidate prophage region{'s' if n != 1 else ''} ({mode} mode), best first:",
             "  " + fmt.format(*head)]
    lines += ["  " + fmt.format(*x) for x in rows]
    if n > limit:
        lines.append(f"  ... and {n - limit} more in regions.tsv")
    lines.append("edges: 'att site' = exact to the base if the repeat is a real att site; "
                 "'estimated' = can be several kb off (regions.tsv gives a range for each edge).")
    lines.append(f"These are candidates to check, not confirmed prophages. How to read them: {GUIDE_URL}")
    return lines
