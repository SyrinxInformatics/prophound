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
(`--mode default --write-seed-scores`), and compares the window scores with the recorded
CUDA ones. Only windows.tsv was pulled from ARC for G5b (calls were rebuilt offline), so
the test is on the decisions the pipeline makes from scores:
  - every window score, ensemble and per checkpoint, against each mode's seed/extend
    thresholds (0.90, 0.95, 0.99, 0.995);
  - every 1 kb bin's mean score against the 0.1 bridge gate.
Everything after those comparisons is integer logic, so if no decision flips between
devices, every mode's calls are identical. Also reports the largest score difference.

The local inputs must match the recorded MD5s (analysis/g5b/input_md5.tsv).
Run from prophound/models/prophage:
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

import numpy as np
import pandas as pd

G5B = Path("analysis/g5b")
CONDITIONS = ("complete", "mag_s0", "mag_s1", "mag_s2", "draft_s0")
CHECKPOINTS = [f"runs/prophage_v2_pool_r{k}/checkpoints/best.pt" for k in (1, 2, 3)]
CONFIG = "runs/prophage_v2_pool_r1/config_snapshot.yaml"
THRESHOLDS = (0.90, 0.95, 0.99, 0.995)   # seed/extend thresholds of all modes and short-contig scales
BRIDGE_PROB = 0.1
BIN_BP = 1000


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


def bin_means(windows: pd.DataFrame) -> np.ndarray:
    """Mean ensemble score of the windows covering each 1 kb bin, all contigs concatenated."""
    out = []
    for _, c in windows.groupby("contig_id", sort=True):
        n = int(np.ceil(c.end.max() / BIN_BP))
        total, count = np.zeros(n), np.zeros(n)
        for start, end, prob in zip(c.start.to_numpy(), c.end.to_numpy(), c.prob.to_numpy()):
            total[start // BIN_BP: -(-end // BIN_BP)] += prob
            count[start // BIN_BP: -(-end // BIN_BP)] += 1
        out.append(total[count > 0] / count[count > 0])
    return np.concatenate(out) if out else np.zeros(0)


def compare(local: Path, cuda: Path) -> dict:
    lw, cw = (pd.read_csv(d / "windows.tsv", sep="\t") for d in (local, cuda))
    coords = ["contig_id", "start", "end"]
    if not lw[coords].astype(str).equals(cw[coords].astype(str)):
        raise SystemExit(f"{local}: the windows differ from the CUDA run's")
    prob_cols = ["prob"] + [c for c in cw.columns if c.startswith("prob_ckpt")]
    a, b = lw[prob_cols].to_numpy(), cw[prob_cols].to_numpy()
    score_flips = int(sum(((a >= t) != (b >= t)).sum() for t in THRESHOLDS))
    bin_flips = int(((bin_means(lw) >= BRIDGE_PROB) != (bin_means(cw) >= BRIDGE_PROB)).sum())
    return {"n_windows": len(lw), "max_score_diff": float(np.abs(a - b).max()),
            "threshold_flips": score_flips, "bridge_flips": bin_flips,
            "decisions_equal": score_flips == 0 and bin_flips == 0}


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
        print(f"[cross-device] {cond}/{genome}: {r['n_windows']} windows, max score diff {r['max_score_diff']:.2e}, "
              f"decision flips {r['threshold_flips']} (thresholds) + {r['bridge_flips']} (bridge gate)", flush=True)
    summary = {"device": args.device, "reference": "CUDA (V100S), g5_prophound.slurm job 826741",
               "n_runs": len(results),
               "runs_with_no_decision_flip": sum(r["decisions_equal"] for r in results.values()),
               "windows": sum(r["n_windows"] for r in results.values()),
               "threshold_flips": sum(r["threshold_flips"] for r in results.values()),
               "bridge_flips": sum(r["bridge_flips"] for r in results.values()),
               "max_score_diff": max(r["max_score_diff"] for r in results.values()), "runs": results}
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(summary, indent=2) + "\n")
    print(f"[done] {summary['runs_with_no_decision_flip']} of {summary['n_runs']} runs with no decision flip; "
          f"max score difference {summary['max_score_diff']:.2e}; wrote {args.out}")


if __name__ == "__main__":
    main()
