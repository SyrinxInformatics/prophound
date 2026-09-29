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
"""Re-derive SHORT_CONTIG_THRESHOLDS under the per-mode seed vote.

short_contig_sweep.py (818235) matched each mode's short-contig threshold to the
false-call load of the intact-genome 10 kb pipeline in the same mode. Since
2026-09-29 the modes carry a seed vote (predict_genome.MODE_MIN_SEED_SUPPORT)
that applies to short-contig calls too, and it lowers the load on both sides.
Same method, vote applied on both sides:

  reference  the 42 intact genomes' 10 kb tracks (r1..r3 + ens), shipped calls
             per mode (mode_recheck.shipped_calls)
  shredded   818235's 3-10 kb fragment tracks (r1..r3 + ens); per threshold t,
             calls at seed = extend = t (no merging), then the mode's vote,
             then mapped back to genome coordinates

The recommended threshold per mode is the lowest t whose false calls per Mb of
element-free sequence do not exceed the reference's. It is printed, not written
into predict_genome.py.

Example:
  python scripts/short_contig_vote_sweep.py --elements scripts/eval_sizestrat_elements.csv \\
      --fasta_dir analysis/short_contig_calib/fasta --shred_root analysis/short_contig_calib \\
      --ref_root analysis/multiscale --ref_suffix _3000 --out analysis/short_contig_calib/vote_sweep.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import mode_recheck as mr  # noqa: E402
from compare_variants import score_variant  # noqa: E402
from refine_sweep import covered50  # noqa: E402
from short_contig_sweep import element_bp, fp_load, fragments  # noqa: E402
from shred_genomes import SEP  # noqa: E402

SEEDS = ("r1", "r2", "r3")


def shred_calls(w: pd.DataFrame, thr: float, frac: float) -> pd.DataFrame:
    regs = []
    for scale, sw in w.groupby("scale"):
        _, r = mr.build_region_table_coverage(sw.drop(columns="scale"), seed_threshold=thr, extend_threshold=thr,
                                              min_coverage_frac=frac, max_gap_bp=0)
        regs.append(r)
    regs = [r for r in regs if len(r)]
    return pd.concat(regs, ignore_index=True) if regs else pd.DataFrame(
        columns=["genome_id", "contig_id", "start", "end", "max_prob"])


def to_genome(reg: pd.DataFrame) -> pd.DataFrame:
    if reg.empty:
        return reg
    parts = reg.contig_id.str.rsplit(SEP, n=1, expand=True)
    return reg.assign(contig_id=parts[0], start=reg.start + parts[1].astype(int), end=reg.end + parts[1].astype(int))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--elements", required=True)
    ap.add_argument("--fasta_dir", required=True, help="shred_genomes.py output")
    ap.add_argument("--shred_root", required=True, help="holds r1/ r2/ r3/ ens/ fragment runs")
    ap.add_argument("--ref_root", required=True, help="holds the intact-genome 10 kb runs")
    ap.add_argument("--ref_suffix", default="", help="run dir suffix under --ref_root, e.g. _3000")
    ap.add_argument("--thresholds", default="0.9,0.95,0.97,0.98,0.985,0.99,0.995,0.998,0.999")
    ap.add_argument("--coverage_frac", type=float, default=0.1)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    ev = pd.read_csv(args.elements)
    genomes = sorted(ev.accession.unique())
    fr = fragments(args.fasta_dir, genomes)
    fr["el_bp"] = [element_bp(ev, f.genome_id, f.contig, f.start, f.end) for f in fr.itertuples()]
    free_bp = float((fr.end - fr.start).sum()) - float((ev.end - ev.start).sum())

    ref = {r: mr.load_tracks(str(Path(args.ref_root) / f"{r}{args.ref_suffix}")) for r in (*SEEDS, "ens")}
    shred = {}
    for r in (*SEEDS, "ens"):
        w = pd.concat([pd.read_csv(p, sep="\t", usecols=["genome_id", "contig_id", "start", "end", "prob"])
                       for p in (Path(args.shred_root) / r).glob("*/windows.tsv")], ignore_index=True)
        shred[r] = w.assign(fasta_path="x", scale=w.end - w.start)
        if not shred[r].genome_id.isin(genomes).all() or shred[r].genome_id.nunique() != len(genomes):
            raise SystemExit(f"{r}: fragment tracks cover {shred[r].genome_id.nunique()} of {len(genomes)} genomes")

    res, flat = {"element_free_bp": free_bp, "modes": {}}, []
    for mode, k in mr.MODE_MIN_SEED_SUPPORT.items():
        rc = mr.shipped_calls(ref["ens"], [ref[s] for s in SEEDS], mode)
        ref_m = {"threshold": mr.MODE_THRESHOLDS[mode], "min_seed_support": k,
                 "covered50_recall": covered50(ev, rc), **fp_load(ev, rc, free_bp)}
        print(f"[{mode} reference 10 kb, t={ref_m['threshold']:g} k={k}] false/Mb "
              f"{ref_m['false_calls_per_mb']:.2f} covered50 {ref_m['covered50_recall']:.3f}", flush=True)
        rows = []
        for thr in [float(x) for x in args.thresholds.split(",")]:
            reg = mr.vote(shred_calls(shred["ens"], thr, args.coverage_frac),
                          [shred_calls(shred[s], thr, args.coverage_frac) for s in SEEDS], k)
            called = set(reg.contig_id)
            reg = to_genome(reg)
            pos = fr[fr.el_bp >= 0.5 * (fr.end - fr.start)]
            m = {"threshold": thr, "covered_recall": score_variant(ev, reg)["covered_recall"] if len(reg) else 0.0,
                 "covered50_recall": covered50(ev, reg), **fp_load(ev, reg, free_bp),
                 "pos_fragment_detect": float(pos.frag.isin(called).mean()),
                 "calls_per_genome": len(reg) / len(genomes)}
            rows.append(m)
            flat.append({"mode": mode, "min_seed_support": k, **m})
            print(f"[{mode} shredded t={thr:g} k={k}] false/Mb {m['false_calls_per_mb']:.2f} "
                  f"covered50 {m['covered50_recall']:.3f} pos-frag {m['pos_fragment_detect']:.3f}", flush=True)
        ok = [r for r in rows if r["false_calls_per_mb"] <= ref_m["false_calls_per_mb"]]
        rec = min(ok, key=lambda r: r["threshold"]) if ok else None
        res["modes"][mode] = {"reference": ref_m, "sweep": rows,
                              "recommended_threshold": rec["threshold"] if rec else None}
        print(f"[{mode}] recommended short-contig threshold: {rec['threshold'] if rec else 'none in grid'} "
              f"(shipped provisional {mr.SHORT_CONTIG_THRESHOLDS[mode]})", flush=True)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(res, indent=2) + "\n")
    pd.DataFrame(flat).to_csv(Path(args.out).with_suffix(".csv"), index=False)
    print(f"[done] wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
