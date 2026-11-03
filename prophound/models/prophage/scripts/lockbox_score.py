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
"""Scores the lockbox once (analysis/lockbox/, analysis/lockbox2/), from the shipped CLI's regions.tsv
per mode (scripts/lockbox_run.slurm). Metric definitions are eval_region_metrics.py's
own functions (covered bp, conflation), so the preregistered definitions are the
code that produced every earlier figure.

Endpoints, as preregistered:
  primary    Default covered50 over all elements, 95% Wilson and genome-bootstrap intervals
  secondary  the same for Exploratory and Conservative; calls/genome; covered (>=1 bp);
             conflated calls and absorbed elements; covered50 per stratum with n;
             the gap to the 42- and 18-genome estimates
Descriptive, not preregistered: calls/genome on the genomes with no element.

With --boundaries (the second lockbox, docs/lockbox2.md), also its boundary endpoints:
median IoU per element; the share of calls with att boundaries; coverage of the model
edges' ranges (start_lo..start_hi, end_lo..end_hi) on calls overlapping one element;
att edges' distance to geNomad's; bp called and bp off every element. The run check
then also requires the 0.1.0 edge stages (refinement, att search, composition).

Before scoring, every run is checked: one regions.tsv and run_metadata.json per genome
and mode, the mode as named, edge snap at 0.5, impute, and one checkpoint-hash set.
Nothing is dropped after scoring. An --out that already exists is never overwritten.

Run from prophound/models/prophage, once:
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
from evalkit import calls_file, covered_bp, load_genome_list, load_pred_regions, overlap, wilson

MODES = ("default", "exploratory", "conservative")
# Final-mode ensemble covered50 on the reused samples (docs/release-text-0.1.0.md
# "Operating modes"; 822223 / 821361), for the optimism gap.
REFERENCE = {
    "42": {"exploratory": 0.994, "default": 0.972, "conservative": 0.944},
    "18": {"exploratory": 1.000, "default": 1.000, "conservative": 0.962},
}
# the 0.1.0 edge stages per mode (run_metadata.json boundary_refinement; docs/lockbox2.md)
# (Exploratory gained composition edges on 2026-10-23, amendment in docs/lockbox2.md)
BOUNDARY_STAGES = dict.fromkeys(MODES, "fine_pass_then_att_then_composition")
STRATA = ["<5kb", "5-8kb", "8-10kb", "10-15kb", "15-30kb", "30-60kb", "60kb+"]


def check_runs(genomes: list, run_dir: Path, mode: str, check_metadata: bool, boundaries: bool = False) -> dict:
    missing, bad = [], []
    shas = set()
    for g in genomes:
        d = run_dir / g
        if not calls_file(d).exists() or not (d / "run_metadata.json").exists():
            missing.append(g)
            continue
        if not check_metadata:
            continue
        m = json.loads((d / "run_metadata.json").read_text())
        if (m.get("mode") != mode or (m.get("edge_snap") or {}).get("frac") != 0.5
                or m.get("ambiguity_policy") != "impute"
                or (boundaries and m.get("boundary_refinement") != BOUNDARY_STAGES[mode])):
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


def boundary_metrics(ev: pd.DataFrame, genomes: list, run_dir: Path) -> dict:
    """The second lockbox's boundary endpoints (docs/lockbox2.md), from full regions.tsv files."""
    ious, cover, att_err = [], [], []
    n_calls = n_att = 0
    bp = off = 0
    for g in genomes:
        tsv = calls_file(run_dir / g)
        pred = pd.read_csv(tsv, sep="\t") if tsv.exists() and tsv.stat().st_size > 1 else pd.DataFrame()
        tr = ev[ev.accession == g]
        for r in tr.itertuples():
            on = pred[(pred.contig_id == r.contig) & (pred.end > r.start) & (pred.start < r.end)] if len(pred) else pred
            ious.append(max(((min(c.end, r.end) - max(c.start, r.start)) / (max(c.end, r.end) - min(c.start, r.start))
                             for c in on.itertuples()), default=0.0))
        if pred.empty:
            continue
        n_calls += len(pred)
        n_att += int((pred.boundary == "att").sum())
        for c in pred.itertuples():
            t = tr[(tr.contig == c.contig_id) & (tr.end > c.start) & (tr.start < c.end)]
            bp += int(c.end - c.start)
            off += int(c.end - c.start) - sum(overlap(x.start, x.end, c.start, c.end) for x in t.itertuples())
            if len(t) != 1:
                continue
            x = t.iloc[0]
            edges = [("start", int(x.start), c.contig_edge not in ("start", "both")),
                     ("end", int(x.end), c.contig_edge not in ("end", "both"))]
            for side, true, real in edges:
                if not real:
                    continue
                if c.boundary == "att":
                    att_err.append(abs(int(getattr(c, side)) - true))
                else:
                    cover.append(int(getattr(c, f"{side}_lo")) <= true <= int(getattr(c, f"{side}_hi")))
    att_err = np.asarray(att_err)
    return {"median_iou": float(np.median(ious)) if ious else None,
            "att_frac_calls": n_att / n_calls if n_calls else None, "n_calls": n_calls, "n_att_calls": n_att,
            "model_edge_range_coverage": float(np.mean(cover)) if cover else None, "n_model_edges": len(cover),
            "att_edge_vs_genomad": {"n": int(len(att_err)),
                                    "median_bp": float(np.median(att_err)) if len(att_err) else None,
                                    "le100": float(np.mean(att_err <= 100)) if len(att_err) else None,
                                    "le1000": float(np.mean(att_err <= 1000)) if len(att_err) else None},
            "Mb_called": bp / 1e6, "Mb_off_element": off / 1e6}


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
    ap.add_argument("--boundaries", action="store_true", help="also the boundary endpoints (docs/lockbox2.md)")
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
        res["modes"][mode] = r = {"runs": check_runs(genomes, run_dir, mode, not a.no_metadata_check, a.boundaries)}
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
        if a.boundaries:
            r["boundaries"] = b = boundary_metrics(ev, genomes, run_dir)
        per_element.append(el.assign(mode=mode))
        print(f"{mode:12s} covered50 {k / n:.3f} ({k}/{n}) Wilson [{r['covered50_wilson95'][0]:.3f}, "
              f"{r['covered50_wilson95'][1]:.3f}] genome-boot [{r['covered50_genome_boot95'][0]:.3f}, "
              f"{r['covered50_genome_boot95'][1]:.3f}] | covered {r['covered']:.3f} | calls/genome "
              f"{r['calls_per_genome']:.1f} (no-element genomes {r['descriptive_phage_free_calls_per_genome']:.1f}) | "
              f"conflated {r['conflated_calls']} absorbing {r['absorbed_elements']} | gap to 42 "
              f"{r['gap_to_reused']['42']:+.3f}, to 18 {r['gap_to_reused']['18']:+.3f}")
        if a.boundaries:
            ae = b["att_edge_vs_genomad"]
            print(f"    boundaries: median IoU {b['median_iou']:.3f} | att on {b['att_frac_calls']:.3f} of calls | "
                  f"model-edge range coverage {b['model_edge_range_coverage']:.3f} (n={b['n_model_edges']}) | "
                  f"att edges vs geNomad median {ae['median_bp']} bp (n={ae['n']}) | "
                  f"{b['Mb_called']:.2f} Mb called, {b['Mb_off_element']:.2f} off")
        for s, v in r["by_stratum"].items():
            print(f"    {s:>8s} n={v['n']:3d} covered50 {v['covered50']:.3f} [{v['wilson95'][0]:.3f}, {v['wilson95'][1]:.3f}]")

    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(res, indent=1))
    pd.concat(per_element, ignore_index=True).to_csv(out.with_suffix(".elements.csv"), index=False)
    print(f"[done] wrote {out} and {out.with_suffix('.elements.csv')}")


if __name__ == "__main__":
    main()
