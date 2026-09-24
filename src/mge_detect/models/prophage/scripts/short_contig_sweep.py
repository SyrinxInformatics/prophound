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
"""Calibrate --short-contigs scan thresholds on shredded genomes.

Runs are `predict_genome --short-contigs scan` over shred_genomes.py output, so
every contig is 3-10 kb and scored at 5 kb or 3 kb windows. windows.tsv keeps
every score, so each threshold is applied offline here. Per scale and threshold:

  - element recall: covered (>=1 bp) and covered50 (>=50% of bp), in genome
    coordinates, overall and per size stratum
  - fragment detection: of fragments that are >=50% element, the fraction with
    a call (what a MAG user sees); of fragments with no element, the fraction
    with a call
  - false-positive load: calls on element-free sequence per Mb, and the fraction
    of element-free bp inside calls

The same false-positive load is computed for --reference runs (the intact-genome
10 kb pipeline at its shipped settings), so a scale's threshold can be chosen to
match the load users already accept at 10 kb, not borrowed from it.

Example:
  python scripts/short_contig_sweep.py --elements scripts/eval_sizestrat_elements.csv \\
      --fasta_dir analysis/short_contig_calib/fasta \\
      --run ens=analysis/short_contig_calib/ens --reference ens=analysis/multiscale/ens_3000 \\
      --out analysis/short_contig_calib/sweep.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from pyfaidx import Fasta

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from compare_variants import score_variant  # noqa: E402
from refine_sweep import covered50, load_run  # noqa: E402
from shred_genomes import SEP  # noqa: E402

try:  # repo layout
    from mge_detect.models.prophage.src.predict_genome import build_region_table_coverage
except ImportError:  # HPCC data_refactored/ layout
    from src.predict_genome import build_region_table_coverage


def element_bp(ev: pd.DataFrame, g: str, c: str, s: int, e: int) -> int:
    x = ev[(ev.accession == g) & (ev.contig == c) & (ev.end > s) & (ev.start < e)]
    return int((np.minimum(x.end, e) - np.maximum(x.start, s)).clip(lower=0).sum())


def fp_load(ev: pd.DataFrame, reg: pd.DataFrame, free_bp: float) -> dict:
    """Calls overlapping no element, per Mb of element-free sequence, and the
    fraction of element-free bp inside calls."""
    n_false, bp_in = 0, 0
    for r in reg.itertuples():
        eb = element_bp(ev, r.genome_id, r.contig_id, int(r.start), int(r.end))
        n_false += eb == 0
        bp_in += (int(r.end) - int(r.start)) - eb
    return {"false_calls_per_mb": n_false / max(1.0, free_bp / 1e6),
            "frac_free_bp_called": bp_in / max(1.0, free_bp)}


def fragments(fasta_dir: str, genomes) -> pd.DataFrame:
    rows = []
    for g in genomes:
        fa = Fasta(str(Path(fasta_dir) / f"{g}.fna"), as_raw=True)
        for fid in fa.keys():
            contig, off = fid.rsplit(SEP, 1)
            rows.append({"genome_id": g, "frag": fid, "contig": contig, "start": int(off),
                         "end": int(off) + len(fa[fid])})
        fa.close()
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--elements", required=True)
    ap.add_argument("--fasta_dir", required=True, help="shred_genomes.py output")
    ap.add_argument("--run", action="append", required=True, metavar="LABEL=DIR")
    ap.add_argument("--reference", action="append", default=[], metavar="LABEL=DIR",
                    help="intact-genome 10 kb runs (windows.tsv + run_metadata.json per genome)")
    ap.add_argument("--thresholds", default="0.9,0.95,0.98,0.99,0.995,0.998,0.999,0.9995")
    ap.add_argument("--coverage_frac", type=float, default=0.1)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    ev = pd.read_csv(args.elements)
    genomes = sorted(ev.accession.unique())
    fr = fragments(args.fasta_dir, genomes)
    fr["el_bp"] = [element_bp(ev, f.genome_id, f.contig, f.start, f.end) for f in fr.itertuples()]
    fr["len"] = fr.end - fr.start
    fr["kind"] = np.where(fr.el_bp == 0, "neg", np.where(fr.el_bp >= 0.5 * fr.len, "pos", "mixed"))
    total_bp = float(fr.len.sum())
    free_bp = total_bp - float(ev.assign(n=ev.end - ev.start).n.sum())
    res = {"n_fragments": int(len(fr)), "fragment_kinds": fr.kind.value_counts().to_dict(),
           "total_bp": total_bp, "element_free_bp": free_bp, "reference": {}, "runs": {}}
    flat = []

    for spec in args.reference:
        label, _, d = spec.partition("=")
        coarse, _, meta, _, _ = load_run(d)
        m = {"covered_recall": score_variant(ev, coarse)["covered_recall"],
             "covered50_recall": covered50(ev, coarse), **fp_load(ev, coarse, free_bp)}
        res["reference"][label] = m
        print(f"[reference {label}] {m}", flush=True)
        flat.append({"run": f"reference_{label}", "scale": 10000, "threshold": meta["seed_threshold"], **m})

    for spec in args.run:
        label, _, d = spec.partition("=")
        w = pd.concat([pd.read_csv(p, sep="\t", usecols=["genome_id", "contig_id", "start", "end", "prob"])
                       for p in Path(d).glob("*/windows.tsv")], ignore_index=True)
        w["fasta_path"] = "x"
        w["scale"] = w.end - w.start
        res["runs"][label] = {}
        for thr in [float(x) for x in args.thresholds.split(",")]:
            regs = []
            for scale, sw in w.groupby("scale"):
                _, r = build_region_table_coverage(sw.drop(columns="scale"), seed_threshold=thr,
                                                   extend_threshold=thr,
                                                   min_coverage_frac=args.coverage_frac, max_gap_bp=0)
                regs.append(r.assign(scale=int(scale)))
            reg = pd.concat(regs, ignore_index=True)
            # back to genome coordinates
            called_frags = set(reg.contig_id)
            if len(reg):
                parts = reg.contig_id.str.rsplit(SEP, n=1, expand=True)
                reg = reg.assign(frag=reg.contig_id, contig_id=parts[0],
                                 start=reg.start + parts[1].astype(int), end=reg.end + parts[1].astype(int))
            fr["called"] = fr.frag.isin(called_frags)
            m = {"covered_recall": score_variant(ev, reg)["covered_recall"],
                 "covered50_recall": covered50(ev, reg), **fp_load(ev, reg, free_bp),
                 "pos_fragment_detect": float(fr[fr.kind == "pos"].called.mean()),
                 "neg_fragment_called": float(fr[fr.kind == "neg"].called.mean()),
                 "calls_per_genome": len(reg) / max(1, len(genomes)),
                 "by_stratum": {str(st): {"covered_recall": score_variant(sub, reg)["covered_recall"],
                                          "covered50_recall": covered50(sub, reg)}
                                for st, sub in ev.groupby("stratum")}}
            # per-scale fragment detection: which window size saw the fragment
            fscale = w.groupby("contig_id").scale.first()
            fr["scale"] = fr.frag.map(fscale).fillna(0).astype(int)
            m["by_scale"] = {str(int(sc)): {"pos_fragment_detect": float(sub[sub.kind == "pos"].called.mean()),
                                            "neg_fragment_called": float(sub[sub.kind == "neg"].called.mean()),
                                            "n_pos": int((sub.kind == "pos").sum()),
                                            "n_neg": int((sub.kind == "neg").sum())}
                             for sc, sub in fr.groupby("scale")}
            res["runs"][label][f"{thr:g}"] = m
            flat.append({"run": label, "scale": "5000+3000", "threshold": thr,
                         **{k: v for k, v in m.items() if not isinstance(v, dict)}})
            print(f"[{label} t={thr:g}] covered {m['covered_recall']:.3f} c50 {m['covered50_recall']:.3f} "
                  f"pos-frag {m['pos_fragment_detect']:.3f} neg-frag {m['neg_fragment_called']:.3f} "
                  f"false/Mb {m['false_calls_per_mb']:.2f} free-bp {m['frac_free_bp_called']:.4f}", flush=True)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(res, indent=2) + "\n")
    pd.DataFrame(flat).to_csv(Path(args.out).with_suffix(".csv"), index=False)
    print(f"[done] wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
