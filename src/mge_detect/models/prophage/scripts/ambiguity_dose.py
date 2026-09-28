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
"""Ambiguity dose-response (docs/v1-scope.md A.3).

The corpus cannot measure how recall degrades with IUPAC ambiguity codes (3.7%
of genomes carry any, at 4.7 ppm), so codes are spiked into clean genomes with
known truth. Substitution only: lengths and coordinates are unchanged, so the
element truth still applies.

  spike  write one genome's FASTA with codes at a given pattern and rate.
         `scattered`: single codes drawn from N,R,Y,K,M,S,W at `rate` ppm
         (the consensus/polishing case behind the EDL933 failure).
         `nruns`: 100 bp runs of N at `rate` per Mb (scaffold gaps).
         Positions come from one seeded draw per genome and pattern, and a
         lower rate uses a prefix of it, so rates are nested and the curve is
         not re-randomised at every step.
  eval   rebuild calls from each condition's score tracks (all modes) and report
         covered / covered50 recall overall and per stratum, calls per genome,
         and elements lost relative to the same policy's clean run.

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
import sys
from pathlib import Path

import numpy as np
import pandas as pd

CODES = np.frombuffer(b"NRYKMSW", dtype=np.uint8)
RUN_BP = 100
MAX_SCATTERED_PPM = 10_000   # draw size cap; the job's rates stay below it
MAX_RUNS_PER_MB = 50


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


def spike_seq(seq: str, genome: str, contig: str, pattern: str, rate: float) -> tuple:
    a = np.frombuffer(seq.upper().encode(), dtype=np.uint8).copy()
    L = len(a)
    rng = rng_for(genome, contig, pattern)
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
    return a.tobytes().decode(), n_bp


def cmd_spike(args) -> None:
    genome = Path(args.fasta).stem
    total, spiked = 0, 0
    with open(args.out, "w") as fh:
        for name, desc, seq in read_fasta(args.fasta):
            s, n = spike_seq(seq, genome, name, args.pattern, args.rate) if args.rate > 0 else (seq, 0)
            total += len(seq)
            spiked += n
            fh.write(f">{desc}\n")
            for i in range(0, len(s), 80):
                fh.write(s[i:i + 80] + "\n")
    print(f"[spike] {genome} {args.pattern} {args.rate:g}: {spiked} bp of {total} ({1e6 * spiked / max(1, total):.1f} ppm)")


def cov_fracs(ev: pd.DataFrame, reg: pd.DataFrame) -> np.ndarray:
    """Per element, the fraction of its bp under the union of overlapping calls."""
    out = []
    for e in ev.itertuples():
        x = reg[(reg.genome_id == e.accession) & (reg.contig_id == e.contig) & (reg.end > e.start) & (reg.start < e.end)]
        iv = sorted((max(int(a), e.start), min(int(b), e.end)) for a, b in zip(x.start, x.end))
        bp, cs, ce = 0, None, None
        for a, b in iv:
            if ce is None or a > ce:
                bp += (ce - cs) if ce is not None else 0
                cs, ce = a, b
            else:
                ce = max(ce, b)
        bp += (ce - cs) if ce is not None else 0
        out.append(bp / (e.end - e.start))
    return np.array(out)


def cmd_eval(args) -> None:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from mode_recheck import calls, load_tracks

    ev = pd.read_csv(args.elements)
    root = Path(args.root)
    res, got = {}, {}
    for d in sorted(p for p in root.iterdir() if p.is_dir() and p.name.count("__") == 2):
        policy, pattern, rate = d.name.split("__")
        tracks = load_tracks(str(d))
        if len(tracks) != ev.accession.nunique():
            print(f"[warn] {d.name}: {len(tracks)} genomes scored of {ev.accession.nunique()}", flush=True)
        for mode in ("exploratory", "default", "conservative"):
            frac = cov_fracs(ev, calls(tracks, mode))
            c50 = frac >= 0.5
            got[(policy, pattern, float(rate), mode)] = c50
            res[f"{d.name}/{mode}"] = {
                "policy": policy, "pattern": pattern, "rate": float(rate), "mode": mode,
                "covered": float((frac > 0).mean()), "covered50": float(c50.mean()),
                "calls_per_genome": float(len(calls(tracks, mode)) / len(tracks)), "n_genomes": len(tracks),
                "by_stratum": {s: float(c50[(ev.stratum == s).to_numpy()].mean()) for s in sorted(ev.stratum.unique())}}
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
    s.add_argument("--pattern", choices=("scattered", "nruns"), required=True)
    s.add_argument("--rate", type=float, required=True, help="ppm (scattered) or runs per Mb (nruns)")
    s.add_argument("--out", required=True)
    e = sub.add_parser("eval")
    e.add_argument("--elements", required=True)
    e.add_argument("--root", required=True, help="dirs named <policy>__<pattern>__<rate> with <genome>/windows.tsv")
    e.add_argument("--out", required=True)
    args = ap.parse_args()
    cmd_spike(args) if args.cmd == "spike" else cmd_eval(args)


if __name__ == "__main__":
    main()
