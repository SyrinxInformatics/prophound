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
"""Ambiguity dose-response (docs/release-scope-0.1.0.md A.3).

The corpus cannot measure how recall degrades with IUPAC ambiguity codes (3.7%
of genomes carry any, at 4.7 ppm), so codes are spiked into clean genomes with
known truth. Substitution only: lengths and coordinates are unchanged, so the
element truth still applies.

  spike  write one genome's FASTA with codes at a given pattern and rate.
         `scattered`: single codes drawn from N,R,Y,K,M,S,W at `rate` ppm
         (the consensus/polishing case behind the EDL933 failure).
         `nruns`: 100 bp runs of N at `rate` per Mb (scaffold gaps).
         `ngap`: scaffold gaps of `rate` bp of N (added 2026-09-29), 2 per Mb at
         random plus, with --elements, one inside every element (capped at 25%
         of its length, placed in its middle half). Gap starts are the same draw
         at every length. --gaps_out records every gap (kind random/element).
         Positions come from one seeded draw per genome and pattern, and a
         lower rate uses a prefix of it, so rates are nested and the curve is
         not re-randomised at every step.
  eval   rebuild calls from each condition's score tracks (all modes) and report
         covered / covered50 recall overall and per stratum, calls per genome,
         and elements lost relative to the same policy's clean run. Runs made
         with --write-seed-scores are scored as shipped (each mode's seed vote);
         with --gaps_dir, `ngap` conditions also report calls lying mostly
         inside random gaps, i.e. calls on imputed sequence.

Example:
  python scripts/ambiguity_dose.py spike --fasta ../raw/genomes/G.fna --pattern scattered \\
      --rate 300 --out /tmp/G.fna
  python scripts/ambiguity_dose.py eval --elements scripts/eval_sizestrat_elements.csv \\
      --root analysis/ambiguity_dose --out analysis/ambiguity_dose/results.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

CODES = np.frombuffer(b"NRYKMSW", dtype=np.uint8)
RUN_BP = 100
MAX_SCATTERED_PPM = 10_000   # draw size cap; the job's rates stay below it
MAX_RUNS_PER_MB = 50
GAPS_PER_MB = 2
EL_GAP_MAX_FRAC = 0.25


def read_fasta(path: str) -> list:
    out, name, desc = [], None, ""
    buf = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    out.append((name, desc, "".join(buf)))
                name, desc, buf = line[1:].split()[0], line[1:].rstrip("\n"), []
            else:
                buf.append(line.strip())
    if name is not None:
        out.append((name, desc, "".join(buf)))
    return out


def rng_for(genome: str, contig: str, pattern: str) -> np.random.Generator:
    h = hashlib.sha256(f"{genome}|{contig}|{pattern}".encode()).digest()
    return np.random.default_rng(int.from_bytes(h[:8], "little"))


def spike_seq(seq: str, genome: str, contig: str, pattern: str, rate: float, elements=None) -> tuple:
    """Returns (sequence, bp changed, gaps) -- gaps as (start, end, kind), ngap only."""
    a = np.frombuffer(seq.upper().encode(), dtype=np.uint8).copy()
    L = len(a)
    rng = rng_for(genome, contig, pattern)
    gaps = []
    if pattern == "ngap":
        g = int(rate)
        n = int(round(GAPS_PER_MB * L / 1e6))
        for st in rng.integers(0, max(1, L - g), size=n):
            gaps.append((int(st), int(min(L, st + g)), "random"))
        for e in ([] if elements is None else elements.itertuples()):
            eg = int(min(g, EL_GAP_MAX_FRAC * (e.end - e.start)))
            er = rng_for(genome, contig, f"ngap_el|{e.start}")
            lo = e.start + (e.end - e.start) // 4
            hi = max(lo + 1, e.start + 3 * (e.end - e.start) // 4 - eg)
            st = int(er.integers(lo, hi))
            gaps.append((st, st + eg, "element"))
        for st, en, _ in gaps:
            a[st:en] = ord("N")
        n_bp = int((a != np.frombuffer(seq.upper().encode(), dtype=np.uint8)).sum())
        return a.tobytes().decode(), n_bp, gaps
    if pattern == "scattered":
        n_max = int(round(MAX_SCATTERED_PPM * L / 1e6))
        pos = rng.choice(L, size=min(n_max, L), replace=False)
        codes = CODES[rng.integers(0, len(CODES), size=len(pos))]
        k = int(round(rate * L / 1e6))
        a[pos[:k]] = codes[:k]
        n_bp = k
    elif pattern == "nruns":
        n_max = int(round(MAX_RUNS_PER_MB * L / 1e6)) + 1
        starts = rng.integers(0, max(1, L - RUN_BP), size=n_max)
        k = int(round(rate * L / 1e6))
        for s in starts[:k]:
            a[s:s + RUN_BP] = ord("N")
        # runs may overlap, so count what actually changed
        n_bp = int((a != np.frombuffer(seq.upper().encode(), dtype=np.uint8)).sum())
    else:
        raise SystemExit(f"unknown pattern {pattern}")
    return a.tobytes().decode(), n_bp, gaps


def cmd_spike(args) -> None:
    genome = Path(args.fasta).stem
    ev = pd.read_csv(args.elements) if args.elements else None
    total, spiked, gap_rows = 0, 0, []
    with open(args.out, "w") as fh:
        for name, desc, seq in read_fasta(args.fasta):
            els = None if ev is None else ev[(ev.accession == genome) & (ev.contig == name)]
            s, n, gaps = spike_seq(seq, genome, name, args.pattern, args.rate, els) if args.rate > 0 else (seq, 0, [])
            gap_rows += [{"genome_id": genome, "contig_id": name, "start": a, "end": b, "kind": k} for a, b, k in gaps]
            total += len(seq)
            spiked += n
            fh.write(f">{desc}\n")
            for i in range(0, len(s), 80):
                fh.write(s[i:i + 80] + "\n")
    if args.gaps_out:
        Path(args.gaps_out).parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(gap_rows, columns=["genome_id", "contig_id", "start", "end", "kind"]).to_csv(
            args.gaps_out, sep="\t", index=False)
    print(f"[spike] {genome} {args.pattern} {args.rate:g}: {spiked} bp of {total} ({1e6 * spiked / max(1, total):.1f} ppm)"
          + (f", {len(gap_rows)} gaps" if gap_rows else ""))


def gap_call_stats(reg: pd.DataFrame, gaps: pd.DataFrame, n_genomes: int) -> dict:
    """Calls with >= 50% of their bp inside random N gaps, and the fraction of
    random-gap bp under any call."""
    n_in, gap_bp, gap_called = 0, 0, 0
    for (g, c), gg in gaps.groupby(["genome_id", "contig_id"]):
        rc = reg[(reg.genome_id == g) & (reg.contig_id == c)]
        iv = sorted(zip(gg.start, gg.end))
        for r in rc.itertuples():
            ov = sum(max(0, min(int(r.end), b) - max(int(r.start), a)) for a, b in iv)
            n_in += ov >= 0.5 * (int(r.end) - int(r.start))
        for a, b in iv:
            gap_bp += b - a
            gap_called += sum(max(0, min(int(r.end), b) - max(int(r.start), a)) for r in rc.itertuples())
    return {"calls_in_gaps_per_genome": n_in / max(1, n_genomes),
            "frac_gap_bp_called": gap_called / max(1, gap_bp)}


def cmd_eval(args) -> None:
    from evalkit import calls, cov_fracs, load_tracks, shipped_calls  # here, so `spike` runs without torch

    ev = pd.read_csv(args.elements)
    root = Path(args.root)
    res, got = {}, {}
    for d in sorted(p for p in root.iterdir() if p.is_dir() and p.name.count("__") == 2):
        policy, pattern, rate = d.name.split("__")
        tracks = load_tracks(str(d))
        if len(tracks) != ev.accession.nunique():
            print(f"[warn] {d.name}: {len(tracks)} genomes scored of {ev.accession.nunique()}", flush=True)
        first = next(iter(d.glob("*/windows.tsv")))
        voted = "prob_ckpt3" in pd.read_csv(first, sep="\t", nrows=0).columns
        seeds = [load_tracks(str(d), f"prob_ckpt{k}") for k in (1, 2, 3)] if voted else None
        gaps = None
        if pattern == "ngap" and args.gaps_dir:
            gp = Path(args.gaps_dir) / f"{pattern}__{rate}"
            gaps = pd.concat([pd.read_csv(f, sep="\t") for f in gp.glob("*.tsv")], ignore_index=True)
            gaps = gaps[gaps.kind == "random"]
        for mode in ("exploratory", "default", "conservative"):
            reg = shipped_calls(tracks, seeds, mode) if voted else calls(tracks, mode)
            frac = cov_fracs(ev, reg)
            c50 = frac >= 0.5
            got[(policy, pattern, float(rate), mode)] = c50
            res[f"{d.name}/{mode}"] = {
                "policy": policy, "pattern": pattern, "rate": float(rate), "mode": mode, "seed_vote": voted,
                "covered": float((frac > 0).mean()), "covered50": float(c50.mean()),
                "calls_per_genome": float(len(reg) / len(tracks)), "n_genomes": len(tracks),
                "by_stratum": {s: float(c50[(ev.stratum == s).to_numpy()].mean()) for s in sorted(ev.stratum.unique())}}
            if gaps is not None:
                # calls lying mostly inside random gaps sit on imputed sequence
                res[f"{d.name}/{mode}"].update(gap_call_stats(reg, gaps, len(tracks)))
    for (policy, pattern, rate, mode), c50 in got.items():
        base = got.get((policy, "clean", 0.0, mode))
        r = res[f"{policy}__{pattern}__{rate:g}/{mode}"]
        if base is not None:
            r["lost_vs_clean"], r["gained_vs_clean"] = int((base & ~c50).sum()), int((~base & c50).sum())
        print(f"[{policy} {pattern} {rate:g} {mode}] covered50 {r['covered50']:.3f} covered {r['covered']:.3f} "
              f"calls/genome {r['calls_per_genome']:.1f} lost {r.get('lost_vs_clean')}", flush=True)
    out = Path(args.out)
    out.write_text(json.dumps(res, indent=2) + "\n")
    flat = [{k: v for k, v in r.items() if k != "by_stratum"} | {f"c50_{s}": v for s, v in r["by_stratum"].items()}
            for r in res.values()]
    pd.DataFrame(flat).sort_values(["mode", "pattern", "policy", "rate"]).to_csv(out.with_suffix(".csv"), index=False)
    print(f"[done] wrote {out}", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("spike")
    s.add_argument("--fasta", required=True)
    s.add_argument("--pattern", choices=("scattered", "nruns", "ngap"), required=True)
    s.add_argument("--rate", type=float, required=True,
                   help="ppm (scattered), runs per Mb (nruns) or gap length in bp (ngap)")
    s.add_argument("--elements", default=None, help="ngap: also put one gap inside every element")
    s.add_argument("--gaps_out", default=None, help="ngap: write the gap intervals here (TSV)")
    s.add_argument("--out", required=True)
    e = sub.add_parser("eval")
    e.add_argument("--elements", required=True)
    e.add_argument("--root", required=True, help="dirs named <policy>__<pattern>__<rate> with <genome>/windows.tsv")
    e.add_argument("--gaps_dir", default=None, help="ngap: <pattern>__<rate>/<genome>.tsv from spike --gaps_out")
    e.add_argument("--out", required=True)
    args = ap.parse_args()
    cmd_spike(args) if args.cmd == "spike" else cmd_eval(args)


if __name__ == "__main__":
    main()
