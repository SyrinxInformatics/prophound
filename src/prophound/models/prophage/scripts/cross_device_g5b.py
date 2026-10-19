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
"""Do the calls on this machine match the recorded CUDA calls? (docs/benchmarks.md)

Re-runs the 45 G5b runs (9 genomes: complete plus four shredded conditions) with the
current code on --device, with the arguments `g5_prophound.slurm` used on ARC's V100
(`--mode default --write-seed-scores`). Then it compares, per run:
  - the calls (contig, start, end) and each call's seed_support, which must be equal;
  - every window score, ensemble and per checkpoint, whose largest difference is reported.

The local inputs must match the recorded MD5s (analysis/g5b/input_md5.tsv).
Run from src/prophound/models/prophage:
  python scripts/cross_device_g5b.py --device mps --workdir <scratch> \
      --out analysis/benchmarks/cross_device_g5b_mps.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd

G5B = Path("analysis/g5b")
CONDITIONS = ("complete", "mag_s0", "mag_s1", "mag_s2", "draft_s0")
CHECKPOINTS = [f"runs/prophage_v2_pool_r{k}/checkpoints/best.pt" for k in (1, 2, 3)]
CONFIG = "runs/prophage_v2_pool_r1/config_snapshot.yaml"
CALL_COLUMNS = ["contig_id", "start", "end", "seed_support"]


def local_input(cond: str, genome: str) -> Path:
    return (G5B / "fasta" if cond == "complete" else G5B / "shreds" / cond) / f"{genome}.fna"


def check_inputs(runs: list[tuple[str, str]]) -> None:
    """Refuse to run unless every local input has the MD5 recorded for the ARC run."""
    recorded = {}
    with open(G5B / "input_md5.tsv") as fh:
        for line in fh:
            if "\t" in line:
                path, md5 = line.rstrip("\n").split("\t")[:2]
                recorded[path.replace("fasta/", "complete/", 1)] = md5   # complete genomes live in fasta/
    bad = []
    for cond, genome in runs:
        key = next((k for k in recorded if k.endswith(f"{cond}/{genome}.fna")), None)
        if key is None or hashlib.md5(local_input(cond, genome).read_bytes()).hexdigest() != recorded[key]:
            bad.append(f"{cond}/{genome}")
    if bad:
        raise SystemExit(f"inputs differ from the ARC runs: {bad[:5]}")


def compare(local: Path, cuda: Path) -> dict:
    la, ca = (pd.read_csv(d / "regions.tsv", sep="\t") for d in (local, cuda))
    lw, cw = (pd.read_csv(d / "windows.tsv", sep="\t") for d in (local, cuda))
    prob_cols = ["prob"] + [c for c in cw.columns if c.startswith("prob_ckpt")]
    same_windows = lw[["contig_id", "start", "end"]].equals(cw[["contig_id", "start", "end"]])
    return {
        "n_calls_local": len(la), "n_calls_cuda": len(ca),
        "calls_equal": la[CALL_COLUMNS].astype(str).equals(ca[CALL_COLUMNS].astype(str)),
        "windows_equal": bool(same_windows),
        "max_score_diff": float(max((lw[c] - cw[c]).abs().max() for c in prob_cols)) if same_windows else None,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--device", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    runs = [(c, g.name) for c in CONDITIONS for g in sorted((G5B / "prophound" / c).iterdir()) if g.is_dir()]
    check_inputs(runs)
    work = Path(args.workdir)
    results = {}
    for cond, genome in runs:
        outdir = work / cond / genome
        cmd = [sys.executable, "-c", "from prophound.cli import main; main()", "predict",
               "--input", str(local_input(cond, genome).resolve()), "--outdir", str(outdir.resolve()),
               "--config", str(Path(CONFIG).resolve()), "--checkpoint", *(str(Path(c).resolve()) for c in CHECKPOINTS),
               "--mode", "default", "--device", args.device, "--write-seed-scores"]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode:
            raise SystemExit(f"{cond}/{genome} failed:\n{proc.stderr[-1500:]}")
        results[f"{cond}/{genome}"] = r = compare(outdir, G5B / "prophound" / cond / genome)
        print(f"[cross-device] {cond}/{genome}: calls {'equal' if r['calls_equal'] else 'DIFFER'} "
              f"({r['n_calls_local']} vs {r['n_calls_cuda']}), max score diff {r['max_score_diff']:.2e}", flush=True)
    summary = {"device": args.device, "reference": "CUDA (V100S), g5_prophound.slurm job 826741",
               "n_runs": len(results),
               "runs_with_equal_calls": sum(r["calls_equal"] for r in results.values()),
               "calls_local": sum(r["n_calls_local"] for r in results.values()),
               "calls_cuda": sum(r["n_calls_cuda"] for r in results.values()),
               "max_score_diff": max(r["max_score_diff"] for r in results.values()), "runs": results}
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(summary, indent=2) + "\n")
    print(f"[done] {summary['runs_with_equal_calls']} of {summary['n_runs']} runs with equal calls; "
          f"max score difference {summary['max_score_diff']:.2e}; wrote {args.out}")


if __name__ == "__main__":
    main()
