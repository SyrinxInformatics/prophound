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
"""Corrective analyses C2, C3, C4 and C8 of docs/research-integrity.md.

  recall       C2 + C3. Shipped modes (mode_recheck.shipped_calls) on a sample
               with r1..r3 + ens run dirs. Per mode: recall at >=1 bp and at
               25 / 50 / 75% of each element's bp under calls, each with a
               binomial Wilson 95% interval (elements independent) and a
               genome-level cluster bootstrap 95% interval (genomes resampled
               with replacement); calls per genome with a bootstrap interval.
               Without --elements (a geNomad-negative sample): calls per genome.
  lit          C4. Literature-pilot sensitivity from <results>.elements.csv and
               the truth file: covered50 per run on all rows, without the
               coordinate classes that rest on a curation judgement, and
               without remnants; genome-bootstrap intervals throughout.
  errortables  C8. Re-derives docs/error-analysis.md's scratchpad-only tables
               from the tracks: seed support of ensemble calls split by
               on/off-target, and the type II miss list with each missed
               element's maximum ensemble window score. Pre-vote Default
               (0.995, no vote) by default, as the doc measured it.

Example:
  python scripts/integrity_checks.py recall --elements scripts/eval_sizestrat_elements.csv \\
      --root analysis/multiscale --suffix _3000 --out analysis/integrity/recall42.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mode_recheck as mr  # noqa: E402
from ambiguity_dose import cov_fracs  # noqa: E402

SEEDS = ("r1", "r2", "r3")
CUTS = {"covered": 0.0, "covered25": 0.25, "covered50": 0.5, "covered75": 0.75}
B = 10_000


def wilson(k: int, n: int, z: float = 1.959964) -> list:
    if n == 0:
        return [float("nan"), float("nan")]
    p = k / n
    c = (p + z * z / (2 * n)) / (1 + z * z / n)
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return [float(c - h), float(c + h)]


def boot_ratio(num: np.ndarray, den: np.ndarray, rng) -> list:
    """95% percentile interval of sum(num)/sum(den) with genomes resampled."""
    idx = rng.integers(0, len(num), size=(B, len(num)))
    r = num[idx].sum(1) / np.maximum(1, den[idx].sum(1))
    return [float(np.percentile(r, 2.5)), float(np.percentile(r, 97.5))]


def hit(frac: np.ndarray, cut: float) -> np.ndarray:
    return frac > 0 if cut == 0 else frac >= cut


def load(root: str, suffix: str) -> tuple:
    tr = {r: mr.load_tracks(str(Path(root) / f"{r}{suffix}")) for r in (*SEEDS, "ens")}
    return tr["ens"], [tr[s] for s in SEEDS]


def cmd_recall(a) -> None:
    rng = np.random.default_rng(a.seed)
    ens, seeds = load(a.root, a.suffix)
    genomes = sorted(ens)
    ev = pd.read_csv(a.elements) if a.elements else None
    res = {"root": a.root, "n_genomes": len(genomes), "bootstrap_resamples": B, "modes": {}}
    rows = []
    for mode in mr.MODE_THRESHOLDS:
        reg = mr.shipped_calls(ens, seeds, mode)
        per_g_calls = reg.groupby("genome_id").size().reindex(genomes, fill_value=0).to_numpy()
        m = {"threshold": mr.MODE_THRESHOLDS[mode], "min_seed_support": mr.MODE_MIN_SEED_SUPPORT[mode],
             "calls_per_genome": float(per_g_calls.mean()),
             "calls_per_genome_ci": boot_ratio(per_g_calls.astype(float), np.ones(len(genomes)), rng)}
        if ev is not None:
            frac = cov_fracs(ev, reg)
            g_of = ev.accession.to_numpy()
            m["n_elements"] = int(len(ev))
            for name, cut in CUTS.items():
                h = hit(frac, cut)
                num = np.array([h[g_of == g].sum() for g in genomes], dtype=float)
                den = np.array([(g_of == g).sum() for g in genomes], dtype=float)
                m[name] = {"value": float(h.mean()), "k": int(h.sum()), "wilson95": wilson(int(h.sum()), len(h)),
                           "genome_bootstrap95": boot_ratio(num, den, rng)}
            m["covered50_by_stratum"] = {s: {"n": int(len(sub)), "value": float(hit(frac[sub.index], 0.5).mean()),
                                             "wilson95": wilson(int(hit(frac[sub.index], 0.5).sum()), len(sub))}
                                         for s, sub in ev.reset_index(drop=True).groupby("stratum")}
        res["modes"][mode] = m
        line = {"mode": mode, "calls_per_genome": m["calls_per_genome"],
                "calls_ci_lo": m["calls_per_genome_ci"][0], "calls_ci_hi": m["calls_per_genome_ci"][1]}
        if ev is not None:
            for name in CUTS:
                line.update({name: m[name]["value"], f"{name}_wilson_lo": m[name]["wilson95"][0],
                             f"{name}_wilson_hi": m[name]["wilson95"][1], f"{name}_boot_lo": m[name]["genome_bootstrap95"][0],
                             f"{name}_boot_hi": m[name]["genome_bootstrap95"][1]})
        rows.append(line)
        print(json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in line.items()}), flush=True)
    write(a.out, res, rows)


def cmd_lit(a) -> None:
    rng = np.random.default_rng(a.seed)
    truth = pd.read_csv(a.truth, sep="\t")
    basis = truth.basis.str.split(";").str[0]
    truth["coord_class"] = np.select(
        [basis.str.startswith("author_feature;") | (basis == "author_feature") | basis.str.startswith("phage_genome"),
         basis.str.startswith("locus_tags")],
        ["author_coordinates", "locus_tag_span"], default="gene_span_judgement")
    el = pd.read_csv(a.elements_csv)
    key = ["genome_id", "name"]
    # class comes from the truth file, so a class change applies without re-scoring
    el = el.drop(columns=["class"]).merge(truth[key + ["coord_class", "class"]], on=key, how="left")
    if el.coord_class.isna().any():
        raise SystemExit("elements.csv rows without a truth match")
    # headline = class `prophage` only; the other classes are reported apart (user, 2026-09-29)
    subsets = {"headline_prophage": lambda d: d[d["class"] == "prophage"],
               "cryptic_prophage": lambda d: d[d["class"] == "cryptic_prophage"],
               "remnant": lambda d: d[d["class"] == "remnant"],
               "phage_derived": lambda d: d[d["class"] == "phage_derived"],
               "all": lambda d: d,
               "author_coordinates_only": lambda d: d[d.coord_class == "author_coordinates"],
               "no_gene_span_judgement": lambda d: d[d.coord_class != "gene_span_judgement"],
               "no_remnants": lambda d: d[~d["class"].isin(["remnant", "phage_derived"])],
               "named_prophages_author_coordinates": lambda d: d[(d["class"] == "prophage")
                                                                 & (d.coord_class == "author_coordinates")]}
    res = {"coord_class_counts": truth.coord_class.value_counts().to_dict(), "runs": {}}
    rows = []
    for run in a.runs.split(","):
        d0 = el[el.run == run]
        if d0.empty:
            raise SystemExit(f"no rows for run {run}")
        res["runs"][run] = {}
        for sname, f in subsets.items():
            d = f(d0)
            gs = sorted(d.genome_id.unique())
            num = np.array([d[d.genome_id == g].covered50.sum() for g in gs], dtype=float)
            den = np.array([(d.genome_id == g).sum() for g in gs], dtype=float)
            v = {"n": int(len(d)), "n_genomes": len(gs), "covered50": float(d.covered50.mean()),
                 "wilson95": wilson(int(d.covered50.sum()), len(d)), "genome_bootstrap95": boot_ratio(num, den, rng)}
            res["runs"][run][sname] = v
            rows.append({"run": run, "subset": sname, **{k: x for k, x in v.items() if not isinstance(x, list)},
                         "boot_lo": v["genome_bootstrap95"][0], "boot_hi": v["genome_bootstrap95"][1]})
            print(f"[{run} {sname}] n={v['n']} covered50 {v['covered50']:.3f} "
                  f"boot [{v['genome_bootstrap95'][0]:.3f}, {v['genome_bootstrap95'][1]:.3f}]", flush=True)
    # which elements each run covers that the reference run does not, with their coordinate class
    if a.reference:
        ref = el[el.run == a.reference].set_index(key)
        for run in a.runs.split(","):
            if run == a.reference:
                continue
            x = el[el.run == run].set_index(key)
            only = x[x.covered50 & ~ref.reindex(x.index).covered50.fillna(False).astype(bool)]
            res["runs"][run]["covered_beyond_reference"] = [
                {"genome_id": g, "name": n, "class": r["class"], "coord_class": r.coord_class}
                for (g, n), r in only.iterrows()]
            # paired difference (run minus reference) per subset, genomes resampled together
            res["runs"][run]["paired_difference_vs_reference"] = {}
            for sname, f in subsets.items():
                d, dr = f(x.reset_index()), f(ref.reset_index())
                m = d.merge(dr[key + ["covered50"]], on=key, suffixes=("", "_ref"))
                gs = sorted(m.genome_id.unique())
                diff = np.array([(m[m.genome_id == g].covered50.astype(int) - m[m.genome_id == g].covered50_ref.astype(int)).sum()
                                 for g in gs], dtype=float)
                den = np.array([(m.genome_id == g).sum() for g in gs], dtype=float)
                v = {"n": int(len(m)), "elements_gained": int((m.covered50 & ~m.covered50_ref).sum()),
                     "elements_lost": int((~m.covered50 & m.covered50_ref).sum()),
                     "difference": float(diff.sum() / den.sum()), "genome_bootstrap95": boot_ratio(diff, den, rng)}
                res["runs"][run]["paired_difference_vs_reference"][sname] = v
                print(f"[{run} - {a.reference} {sname}] +{v['elements_gained']} / -{v['elements_lost']} "
                      f"diff {v['difference']:+.3f} boot [{v['genome_bootstrap95'][0]:+.3f}, "
                      f"{v['genome_bootstrap95'][1]:+.3f}]", flush=True)
    write(a.out, res, rows)


def cmd_errortables(a) -> None:
    ens, seeds = load(a.root, a.suffix)
    ev = pd.read_csv(a.elements)
    saved = dict(mr.MODE_THRESHOLDS)
    mr.MODE_THRESHOLDS["default"] = a.threshold
    try:
        reg = mr.calls(ens, "default")
        sreg = [mr.calls(s, "default") for s in seeds]
    finally:
        mr.MODE_THRESHOLDS.update(saved)
    try:
        from mge_detect.models.prophage.src.predict_genome import seed_support
    except ImportError:
        from src.predict_genome import seed_support
    support = []
    for g, sub in reg.groupby("genome_id", sort=False):
        support += seed_support(sub, [s[s.genome_id == g] for s in sreg])
    reg = reg.assign(seed_support=support)
    reg["on_target"] = [mr.element_overlap_bp(ev, r) > 0 for r in reg.itertuples()]
    table = reg.groupby(["seed_support", "on_target"]).size().unstack(fill_value=0)
    table.columns = ["off_target" if not c else "on_target" for c in table.columns]
    frac = cov_fracs(ev, reg)
    misses = []
    for e, f in zip(ev.itertuples(), frac):
        if f >= 0.5:
            continue
        w = ens[e.accession][0]
        w = w[(w.contig_id == e.contig) & (w.end > e.start) & (w.start < e.end)]
        misses.append({"accession": e.accession, "contig": e.contig, "start": int(e.start), "end": int(e.end),
                       "stratum": e.stratum, "cov_frac": float(f), "max_window_prob": float(w.prob.max()) if len(w) else None})
    res = {"threshold": a.threshold, "vote": "none", "n_calls": int(len(reg)),
           "seed_support_by_target": table.reset_index().to_dict(orient="records"),
           "covered50_misses": misses}
    print(table.to_string(), flush=True)
    print(f"[misses] {len(misses)}: max window prob " + ", ".join(f"{m['max_window_prob']:.3f}" for m in misses), flush=True)
    write(a.out, res, None)


def write(out: str, res: dict, rows) -> None:
    p = Path(out)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(res, indent=2, default=str) + "\n")
    if rows:
        pd.DataFrame(rows).to_csv(p.with_suffix(".csv"), index=False)
    print(f"[done] wrote {out}", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("recall")
    r.add_argument("--elements", default=None)
    r.add_argument("--root", required=True, help="holds r1 r2 r3 ens run dirs (+ --suffix)")
    r.add_argument("--suffix", default="")
    r.add_argument("--seed", type=int, default=20260929)
    r.add_argument("--out", required=True)
    l = sub.add_parser("lit")
    l.add_argument("--truth", required=True)
    l.add_argument("--elements_csv", required=True, help="literature_pilot_eval <out>.elements.csv")
    l.add_argument("--runs", required=True, help="comma list of run labels, e.g. ens/default,genomad")
    l.add_argument("--reference", default=None, help="run to compare against, e.g. genomad")
    l.add_argument("--seed", type=int, default=20260929)
    l.add_argument("--out", required=True)
    e = sub.add_parser("errortables")
    e.add_argument("--elements", required=True)
    e.add_argument("--root", required=True)
    e.add_argument("--suffix", default="")
    e.add_argument("--threshold", type=float, default=0.995, help="pre-vote Default, as error-analysis measured")
    e.add_argument("--out", required=True)
    a = ap.parse_args()
    {"recall": cmd_recall, "lit": cmd_lit, "errortables": cmd_errortables}[a.cmd](a)


if __name__ == "__main__":
    main()
