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
"""Scores the lockbox once (docs/lockbox.md), from the shipped CLI's regions.tsv
per mode (scripts/lockbox_run.slurm). Metric definitions are eval_region_metrics.py's
own functions (covered bp, conflation), so the preregistered definitions are the
code that produced every earlier figure.

Endpoints, as preregistered:
  primary    Default covered50 over all elements, 95% Wilson and genome-bootstrap intervals
  secondary  the same for Exploratory and Conservative; calls/genome; covered (>=1 bp);
             conflated calls and absorbed elements; covered50 per stratum with n;
             the gap to the 42- and 18-genome estimates
Descriptive, not preregistered: calls/genome on the genomes with no element.

Before scoring, every run is checked: one regions.tsv and run_metadata.json per genome
and mode, the mode as named, edge snap at 0.5, impute, and one checkpoint-hash set.
Nothing is dropped after scoring. An --out that already exists is never overwritten.

Run from src/prophound/models/prophage, once:
  python scripts/lockbox_score.py --elements analysis/lockbox/lockbox_elements.csv \\
      --genome_list analysis/lockbox/lockbox_genome_list.txt \\
      --run default=analysis/lockbox/prophound/default \\
      --run exploratory=analysis/lockbox/prophound/exploratory \\
      --run conservative=analysis/lockbox/prophound/conservative \\
      --out analysis/lockbox/results.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from evalkit import covered_bp, load_genome_list, load_pred_regions, overlap, wilson

MODES = ("default", "exploratory", "conservative")
# Final-mode ensemble covered50 on the reused samples (docs/release-checklist.md,
# "Operating modes in plain language"; 822223 / 821361), for the optimism gap.
REFERENCE = {
    "42": {"exploratory": 0.994, "default": 0.972, "conservative": 0.944},
    "18": {"exploratory": 1.000, "default": 1.000, "conservative": 0.962},
}
STRATA = ["<5kb", "5-8kb", "8-10kb", "10-15kb", "15-30kb", "30-60kb", "60kb+"]


def check_runs(genomes: list, run_dir: Path, mode: str, check_metadata: bool) -> dict:
    missing, bad = [], []
    shas = set()
    for g in genomes:
        d = run_dir / g
        if not (d / "regions.tsv").exists() or not (d / "run_metadata.json").exists():
            missing.append(g)
            continue
        if not check_metadata:
            continue
        m = json.loads((d / "run_metadata.json").read_text())
        if (m.get("mode") != mode or (m.get("edge_snap") or {}).get("frac") != 0.5
                or m.get("ambiguity_policy") != "impute"):
            bad.append(g)
        shas.add(json.dumps(m.get("checkpoint_sha256"), sort_keys=True))
    if missing or bad or len(shas) > 1:
        raise SystemExit(f"{run_dir}: {len(missing)} genome(s) without outputs (e.g. {missing[:2]}), "
                         f"{len(bad)} with wrong mode/snap/policy (e.g. {bad[:2]}), "
                         f"{len(shas)} checkpoint-hash sets. Fix or exclude before scoring.")
    return {"n_runs": len(genomes), "checkpoint_sha256": json.loads(shas.pop()) if shas else None}


def score_mode(ev: pd.DataFrame, genomes: list, run_dir: Path) -> tuple[pd.DataFrame, dict]:
    rows, calls = [], {}
    n_conflated = n_absorbed = 0
    for g in genomes:
        pred = load_pred_regions(run_dir, g)
        calls[g] = len(pred)
        tr = ev[ev.accession == g]
        for r in tr.itertuples():
            on = pred[pred.contig_id == r.contig]
            iv = list(zip(on.start.astype(int), on.end.astype(int)))
            rows.append({"accession": g, "contig": r.contig, "start": int(r.start), "end": int(r.end),
                         "stratum": r.stratum,
                         "covered_frac": covered_bp(int(r.start), int(r.end), iv) / max(1, int(r.end - r.start)),
                         "n_calls": sum(1 for s, e in iv if overlap(r.start, r.end, s, e) > 0)})
        # conflation, as in eval_region_metrics.evaluate: a call overlapping >=2 elements
        for p in pred.itertuples():
            k = sum(1 for t in tr.itertuples() if t.contig == p.contig_id and overlap(t.start, t.end, p.start, p.end) > 0)
            if k >= 2:
                n_conflated += 1
                n_absorbed += k
    el = pd.DataFrame(rows)
    el["covered50"] = el.covered_frac >= 0.5
    el["covered"] = el.n_calls > 0
    return el, {"calls": calls, "conflated_calls": n_conflated, "absorbed_elements": n_absorbed}


def genome_bootstrap(el: pd.DataFrame, genomes: list, n_boot: int, seed: int) -> list:
    per = el.groupby("accession").covered50.agg(["sum", "count"]).reindex(genomes, fill_value=0)
    k, n = per["sum"].to_numpy(), per["count"].to_numpy()
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(genomes), size=(n_boot, len(genomes)))
    den = n[idx].sum(1)
    est = k[idx].sum(1)[den > 0] / den[den > 0]
    return [float(np.percentile(est, 2.5)), float(np.percentile(est, 97.5))]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--elements", required=True, help="accession,contig,start,end,...,stratum")
    ap.add_argument("--genome_list", required=True)
    ap.add_argument("--run", action="append", required=True, metavar="MODE=DIR", help="DIR/<genome>/regions.tsv")
    ap.add_argument("--out", required=True)
    ap.add_argument("--n_boot", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=20260929)
    ap.add_argument("--no_metadata_check", action="store_true",
                    help="validation on older runs only (they predate edge snap); never for the lockbox")
    a = ap.parse_args()

    out = Path(a.out)
    if out.exists():
        raise SystemExit(f"{out} exists: the lockbox is scored once, and its result is never overwritten")
    genomes = load_genome_list(Path(a.genome_list))
    ev = pd.read_csv(a.elements)
    stray = sorted(set(ev.accession) - set(genomes))
    if stray:
        raise SystemExit(f"{len(stray)} element genome(s) not in the genome list, e.g. {stray[:2]}")
    runs = dict(s.split("=", 1) for s in a.run)
    if set(runs) - set(MODES):
        raise SystemExit(f"--run modes must be among {MODES}")
    phage_free = [g for g in genomes if g not in set(ev.accession)]

    res = {"n_genomes": len(genomes), "n_elements": int(len(ev)), "n_genomes_with_elements": int(ev.accession.nunique()),
           "n_genomes_without_elements": len(phage_free), "modes": {}}
    per_element = []
    for mode in [m for m in MODES if m in runs]:
        run_dir = Path(runs[mode])
        res["modes"][mode] = r = {"runs": check_runs(genomes, run_dir, mode, not a.no_metadata_check)}
        el, extra = score_mode(ev, genomes, run_dir)
        k, n = int(el.covered50.sum()), len(el)
        n_calls = sum(extra["calls"].values())
        r.update({
            "covered50": k / n, "k": k, "n": n,
            "covered50_wilson95": wilson(k, n),
            "covered50_genome_boot95": genome_bootstrap(el, genomes, a.n_boot, a.seed),
            "covered": float(el.covered.mean()),
            "calls_per_genome": n_calls / len(genomes),
            "conflated_calls": extra["conflated_calls"], "absorbed_elements": extra["absorbed_elements"],
            "by_stratum": {s: {"n": int((el.stratum == s).sum()), "k": int(el[el.stratum == s].covered50.sum()),
                               "covered50": float(el[el.stratum == s].covered50.mean()),
                               "wilson95": wilson(int(el[el.stratum == s].covered50.sum()), int((el.stratum == s).sum()))}
                           for s in STRATA if (el.stratum == s).any()},
            "gap_to_reused": {lab: ref[mode] - k / n for lab, ref in REFERENCE.items()},
            "descriptive_phage_free_calls_per_genome":
                float(np.mean([extra["calls"][g] for g in phage_free])) if phage_free else float("nan"),
        })
        per_element.append(el.assign(mode=mode))
        print(f"{mode:12s} covered50 {k / n:.3f} ({k}/{n}) Wilson [{r['covered50_wilson95'][0]:.3f}, "
              f"{r['covered50_wilson95'][1]:.3f}] genome-boot [{r['covered50_genome_boot95'][0]:.3f}, "
              f"{r['covered50_genome_boot95'][1]:.3f}] | covered {r['covered']:.3f} | calls/genome "
              f"{r['calls_per_genome']:.1f} (no-element genomes {r['descriptive_phage_free_calls_per_genome']:.1f}) | "
              f"conflated {r['conflated_calls']} absorbing {r['absorbed_elements']} | gap to 42 "
              f"{r['gap_to_reused']['42']:+.3f}, to 18 {r['gap_to_reused']['18']:+.3f}")
        for s, v in r["by_stratum"].items():
            print(f"    {s:>8s} n={v['n']:3d} covered50 {v['covered50']:.3f} [{v['wilson95'][0]:.3f}, {v['wilson95'][1]:.3f}]")

    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(res, indent=1))
    pd.concat(per_element, ignore_index=True).to_csv(out.with_suffix(".elements.csv"), index=False)
    print(f"[done] wrote {out} and {out.with_suffix('.elements.csv')}")


if __name__ == "__main__":
    main()
