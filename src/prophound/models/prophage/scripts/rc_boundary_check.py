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
"""Reverse-complement symmetry of the whole call pipeline, boundary stages included.

Release text §7 and docs/scientific-validity.md "Input-robustness standing check": RC-1
(window scores) holds by construction, but whole-contig calls (RC-2) were last checked
before refinement, att sites and composition edges existed. This runs the shipped
`predict` on each genome and on its reverse complement and compares the calls.

Every contig is first trimmed to a multiple of 1 kb, so the 1 kb coarse grid and the fine
grid mirror exactly (otherwise the ragged last window moves to the start under RC;
rc_n_resilience.py RC-2). Contig ids and order are kept. RC calls are mapped back
(start = L - end, end = L - start) and compared stage by stage: the coarse call
(`coarse_start`/`coarse_end`, or the call itself when refinement left it whole), the
model edges after refinement (`model_start`/`model_end`), and the final call, with its
boundary, edge sources (mirrored: forward start = RC end) and att core (reverse
complemented). Window and fine-window scores are compared mirrored as well.

Outputs in --out: rc_boundary_check.calls.tsv (one row per forward or RC call) and
rc_boundary_check.summary.tsv (per genome and mode).

  python scripts/rc_boundary_check.py --fasta a.fna b.fna --modes default exploratory --out DIR
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

GRID = 1000
_COMP = str.maketrans("ACGTRYKMBVDHNacgtrykmbvdhn", "TGCAYRMKVBHDNtgcayrmkvbhdn")


def rc(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def read_fasta(path: Path) -> list[tuple[str, str]]:
    recs, name, buf = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(buf)))
                name, buf = line[1:], []
            else:
                buf.append(line.strip())
    if name is not None:
        recs.append((name, "".join(buf)))
    return recs


def write_fasta(path: Path, recs: list[tuple[str, str]]) -> None:
    with open(path, "w") as fh:
        for name, seq in recs:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")


def run_predict(fasta: Path, outdir: Path, mode: str) -> None:
    if (outdir / "run_metadata.json").exists():
        return
    cmd = [sys.executable, "-c", "from prophound.cli import main; main()", "predict",
           "--input", str(fasta), "--outdir", str(outdir), "--mode", mode, "-q", "--force"]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)


def mirror(df: pd.DataFrame, lengths: dict[str, int], pairs: list[tuple[str, str]]) -> pd.DataFrame:
    """Map RC coordinates back onto the forward strand: [s, e) -> [L - e, L - s)."""
    out = df.copy()
    L = out["contig_id"].map(lengths).astype(np.int64)
    for a, b in pairs:
        if a not in out.columns:
            continue
        va, vb = df[a], df[b]
        out[a] = L - vb
        out[b] = L - va
    return out


def stage_extent(df: pd.DataFrame, a: str, b: str) -> pd.Series:
    has = df[a].notna() if a in df.columns else pd.Series(False, index=df.index)
    s = np.where(has, df.get(a, df["start"]), df["start"]).astype(np.int64)
    e = np.where(has, df.get(b, df["end"]), df["end"]).astype(np.int64)
    return pd.Series(list(zip(s, e)), index=df.index)


def compare_tracks(fwd: Path, rev: Path, lengths: dict[str, int]) -> float:
    if not (fwd.exists() and rev.exists()):
        return float("nan")
    f = pd.read_csv(fwd, sep="\t")
    r = mirror(pd.read_csv(rev, sep="\t"), lengths, [("start", "end")])
    m = f.merge(r, on=["contig_id", "start", "end"], how="outer", suffixes=("_f", "_r"))
    if m["prob_f"].isna().any() or m["prob_r"].isna().any():
        return float("inf")          # the grids do not mirror
    return float((m["prob_f"] - m["prob_r"]).abs().max())


def compare_calls(f: pd.DataFrame, r: pd.DataFrame, genome: str, mode: str) -> pd.DataFrame:
    """Pair each forward call with the mirrored RC call of greatest overlap."""
    rows = []
    used = set()
    for i, a in f.iterrows():
        cand = r[(r["contig_id"] == a["contig_id"]) & (r["end"] > a["start"]) & (r["start"] < a["end"])]
        cand = cand[~cand.index.isin(used)]
        base = dict(genome=genome, mode=mode, contig_id=a["contig_id"], f_start=a["start"], f_end=a["end"],
                    f_boundary=a["boundary"])
        if cand.empty:
            rows.append(base | dict(status="forward_only"))
            continue
        ov = np.minimum(cand["end"], a["end"]) - np.maximum(cand["start"], a["start"])
        j = ov.idxmax()
        used.add(j)
        b = r.loc[j]
        same = lambda k: (a[k] == b[k]) or (pd.isna(a[k]) and pd.isna(b[k]))
        d = dict(r_start=b["start"], r_end=b["end"], r_boundary=b["boundary"],
                 d_start=int(b["start"] - a["start"]), d_end=int(b["end"] - a["end"]),
                 coarse_same=a["_coarse"] == b["_coarse"], model_same=a["_model"] == b["_model"],
                 final_same=(a["start"], a["end"]) == (b["start"], b["end"]),
                 boundary_same=same("boundary"),
                 sources_same=same("start_source") and same("end_source"),
                 att_same=same("att_core"), support_same=same("seed_support"))
        first = ("same" if d["final_same"] and d["boundary_same"] and d["sources_same"] and d["att_same"]
                 else "coarse" if not d["coarse_same"]
                 else "refinement" if not d["model_same"]
                 else "att/composition")
        rows.append(base | d | dict(status="paired", first_diff=first,
                                    f_sources=f"{a['start_source']}/{a['end_source']}",
                                    r_sources=f"{b['start_source']}/{b['end_source']}"))
    for j, b in r.iterrows():
        if j not in used:
            rows.append(dict(genome=genome, mode=mode, contig_id=b["contig_id"], r_start=b["start"],
                             r_end=b["end"], r_boundary=b["boundary"], status="rc_only"))
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", nargs="+", required=True, type=Path)
    ap.add_argument("--modes", nargs="+", default=["default", "exploratory"])
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    all_calls, summary = [], []
    for fa in args.fasta:
        genome = fa.name.split(".fna")[0].split(".fa")[0]
        gdir = args.out / genome
        gdir.mkdir(exist_ok=True)
        recs = [(n, s[: len(s) - len(s) % GRID]) for n, s in read_fasta(fa)]
        lengths = {n.split()[0]: len(s) for n, s in recs}
        fwd_fa, rev_fa = gdir / f"{genome}.fna", gdir / f"{genome}.rc.fna"
        if not fwd_fa.exists():
            write_fasta(fwd_fa, recs)
            write_fasta(rev_fa, [(n, rc(s)) for n, s in recs])
        for mode in args.modes:
            fd, rd = gdir / f"{mode}_fwd", gdir / f"{mode}_rc"
            run_predict(fwd_fa, fd, mode)
            run_predict(rev_fa, rd, mode)
            f = pd.read_csv(fd / "regions.detailed.tsv", sep="\t")
            r0 = pd.read_csv(rd / "regions.detailed.tsv", sep="\t")
            r = mirror(r0, lengths, [("start", "end"), ("coarse_start", "coarse_end"),
                                     ("model_start", "model_end"), ("start_lo", "end_hi"),
                                     ("start_hi", "end_lo")])
            r["start_source"], r["end_source"] = r0["end_source"], r0["start_source"]
            if "att_core" in r.columns:
                r["att_core"] = r0["att_core"].map(lambda s: rc(s) if isinstance(s, str) else s)
            for df in (f, r):
                df["_coarse"] = stage_extent(df, "coarse_start", "coarse_end")
                df["_model"] = stage_extent(df, "model_start", "model_end")
            c = compare_calls(f, r, genome, mode)
            all_calls.append(c)
            paired = c[c["status"] == "paired"]
            summary.append(dict(
                genome=genome, mode=mode, fwd_calls=len(f), rc_calls=len(r),
                paired=len(paired), identical=int((paired["first_diff"] == "same").sum()),
                first_diff_coarse=int((paired["first_diff"] == "coarse").sum()),
                first_diff_refinement=int((paired["first_diff"] == "refinement").sum()),
                first_diff_att_composition=int((paired["first_diff"] == "att/composition").sum()),
                forward_only=int((c["status"] == "forward_only").sum()),
                rc_only=int((c["status"] == "rc_only").sum()),
                att_calls_fwd=int((f["boundary"] == "att").sum()), att_calls_rc=int((r["boundary"] == "att").sum()),
                max_window_diff=compare_tracks(fd / "windows.tsv", rd / "windows.tsv", lengths),
                max_fine_diff=compare_tracks(fd / "fine_windows.tsv", rd / "fine_windows.tsv", lengths)))
            print(pd.DataFrame([summary[-1]]).to_string(index=False, header=len(summary) == 1), flush=True)

    calls = pd.concat(all_calls, ignore_index=True)
    calls.to_csv(args.out / "rc_boundary_check.calls.tsv", sep="\t", index=False)
    s = pd.DataFrame(summary)
    s.to_csv(args.out / "rc_boundary_check.summary.tsv", sep="\t", index=False)
    print("\n" + s.to_string(index=False))


if __name__ == "__main__":
    main()
