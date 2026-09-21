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



"""Downstream quick checks on already-scored windows.tsv files (no re-inference).

Written 2026-09-17 for the error analysis in docs/model-quality-baseline.md and
the design in docs/multi-model-architecture.md. Compares, on the 18-genome
eval set, the frozen `build_region_table` (intersection semantics: a 1 kb bin
is in only if the MEAN of the ~10 windows covering it >= extend_threshold)
against a coverage-fraction builder (bin is in if >= `frac` of its covering
windows individually score >= extend_threshold; frac->0 is the union), with
optional replicate ensembling (mean / max of window probs) and a minimum
region length. The coverage builder is `predict_genome.build_region_table_coverage`
(`mge-detect predict --region-builder coverage --coverage-frac F`); this script
just sweeps it. Each variant's regions.tsv is written under --out_dir/<variant>/
<genome>/ and scored with eval_region_metrics.py.

Usage:
    python coverage_builder_checks.py --predict_dirs A=/path/predict_A [B=/path/predict_B ...] \
        --out_dir /scratch/variants [--fracs 0.01,0.1,0.3] [--modes default,exploratory] [--min_len 0]

Each --predict_dirs entry must contain <genome>/windows.tsv per genome in
eval_genome_list.txt (i.e. `mge-detect predict` output dirs). With >=2 entries,
mean- and max-ensembles across them are added as extra inputs.
"""
from __future__ import annotations

import argparse
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from src.predict_genome import build_region_table, build_region_table_coverage  # noqa: E402

MODES = {"conservative": (0.999, 0.99), "default": (0.995, 0.99), "exploratory": (0.95, 0.95)}


def coverage_builder(
    w: pd.DataFrame, frac: float, seed: float, ext: float, gap: int = 0, bridge: float = 0.0
) -> pd.DataFrame:
    _, regions = build_region_table_coverage(
        w, seed_threshold=seed, extend_threshold=ext, min_coverage_frac=frac,
        max_gap_bp=gap, bridge_prob=bridge,
    )
    return regions


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--predict_dirs", nargs="+", required=True, help="NAME=/path/to/predict_dir entries")
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--genome_list", default=str(HERE / "eval_genome_list.txt"))
    ap.add_argument("--genomad_dir", default=str(HERE / "eval_genomad_truth"))
    ap.add_argument("--fracs", default="0.01,0.1,0.3,0.5")
    ap.add_argument("--modes", default="default")
    ap.add_argument("--gaps", default="", help="Fragment-merge sweep: comma-separated max gap sizes in bp (e.g. '5000,20000,50000'). Empty = no merge variants (docs/gap-filling.md).")
    ap.add_argument("--bridges", default="0", help="With --gaps: comma-separated bridge_prob gates (mean covering-window prob required across the gap). 0 = distance-only merge.")
    ap.add_argument("--min_len", type=int, default=0, help="Also emit each variant filtered to regions >= this many bp (0 = off)")
    ap.add_argument("--strata", default="", help="Passed through to eval_region_metrics.py: kb breakpoints for per-size-class recall (e.g. '5,8,10,15,30,60').")
    args = ap.parse_args()

    genomes = [l.strip() for l in open(args.genome_list) if l.strip()]
    inputs = dict(e.split("=", 1) for e in args.predict_dirs)
    fracs = [float(f) for f in args.fracs.split(",")]
    modes = args.modes.split(",")
    gaps = [int(g) for g in args.gaps.split(",") if g.strip()]
    bridges = [float(b) for b in args.bridges.split(",") if b.strip()]
    out = Path(args.out_dir)
    variants: dict = {}

    for g in genomes:
        tracks = {n: pd.read_csv(Path(p) / g / "windows.tsv", sep="\t") for n, p in inputs.items()}
        if len(tracks) >= 2:
            base = next(iter(tracks.values()))
            stack = np.stack([t.prob.values for t in tracks.values()])
            for name, arr in (("ensmean", stack.mean(0)), ("ensmax", stack.max(0))):
                t = base.copy()
                t["prob"] = arr
                tracks[name] = t
        tables = {}
        for name, w in tracks.items():
            for mode in modes:
                seed, ext = MODES[mode]
                _, frozen = build_region_table(w, seed_threshold=seed, extend_threshold=ext)
                tables[f"{name}_frozen_{mode}"] = frozen
                for f in fracs:
                    tables[f"{name}_cov{f}_{mode}"] = coverage_builder(w, f, seed, ext)
                    for gap in gaps:
                        for br in bridges:
                            key = f"{name}_cov{f}_{mode}_gap{gap // 1000}k"
                            if br > 0:
                                key += f"br{br:g}"
                            tables[key] = coverage_builder(w, f, seed, ext, gap=gap, bridge=br)
        if args.min_len > 0:
            for k, t in list(tables.items()):
                tables[f"{k}_min{args.min_len}"] = t[t.region_len_bp >= args.min_len] if len(t) else t
        for k, t in tables.items():
            d = out / k / g
            d.mkdir(parents=True, exist_ok=True)
            t.to_csv(d / "regions.tsv", sep="\t", index=False)
            variants[k] = True

    strata_rows: dict = {}
    print(f"{'variant':40s} {'recall':>14s} {'precision':>15s} {'medIoU':>7s} {'start -/+':>10s} {'end -/+':>9s}")
    for k in variants:
        cmd = [sys.executable, str(HERE / "eval_region_metrics.py"), "--genome_list", args.genome_list,
               "--genomad_dir", args.genomad_dir, "--predict_dir", str(out / k), "--label", k]
        if args.strata:
            cmd += ["--strata", args.strata]
        res = subprocess.run(cmd, capture_output=True, text=True).stdout.splitlines()
        strata_rows[k] = [l for l in res if l.split() and l.split()[0].endswith("kb") and "%" in l]
        first = lambda p: next(l for l in res if l.startswith(p))  # noqa: E731
        ss = first("start error sign").split(":")[1].split("/")
        es = first("end error sign").split(":")[1].split("/")
        print(f"{k:40s} {first('region-level recall').split(':')[1].strip():>14s} "
              f"{first('region-level precision').split(':')[1].strip():>15s} "
              f"{first('IoU:').split('median=')[1].split()[0]:>7s} "
              f"{ss[0].split()[0]}/{ss[1].split()[0]:>6s} {es[0].split()[0]}/{es[1].split()[0]:>5s}")


    if args.strata:
        print("\n=== recall by element size class, per variant ===")
        for k, rows in strata_rows.items():
            if not rows:
                continue
            compact = "  ".join(
                f"{r.split()[0]}={r.split()[3]}/IoU{r.split()[4]}" for r in rows
            )
            print(f"{k:40s} {compact}")


if __name__ == "__main__":
    main()
