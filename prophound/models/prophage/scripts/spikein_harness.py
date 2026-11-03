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


"""Synthetic spike-in harness: measure detection and delimitation at any element size.

The corpus cannot answer "how small an
element can the pipeline resolve" -- geNomad's floor is 3.50 kb and only 2 of
its 3,970 elements are below 5 kb -- so this builds the test cases instead:
real element sequence of a chosen size, implanted at a known coordinate in
verified host, scored at a chosen window scale.

That gives *n* by choice rather than by inheritance, at any size, and it
measures the two limits separately (they are not the same number):

  * **detection** -- does any window over the implant clear a threshold
  * **delimitation** -- how close are the called edges to the known ones

Verified host matters. An earlier hand-run of this experiment picked "host" that
overlapped a known prophage and reported a 700 bp element scoring 1.0000; with
genuine host the same case scores 0.077. So every host slice is scored first and
rejected if it is not quiet (`--max_host_prob`).

**Type-agnostic by design.** `--element_fasta` takes any sequence, so the same
harness serves IS elements, ICE or integrons later without change -- the element
type axis of the design doc. With no `--element_fasta` it carves implants from
geNomad-called prophage regions of the donor genome.

What it does NOT establish: the implant context is artificial (no integration
site, no host adaptation, composition inherited from the donor), so it answers
"what size can this pipeline resolve", not "does it find real IS elements".
That needs curated truth or an orthogonal database, held out.

Example:
  python scripts/spikein_harness.py \
      --host_fasta raw/genomes/GCF_000006665.1_ASM666v1_genomic.fna \
      --host_contig NC_002655.2 --host_start 2400000 --host_end 2440000 \
      --donor_tsv scripts/eval_genomad_truth/<acc>/<acc>_find_proviruses/<acc>_provirus.tsv \
      --config runs/prophage_v2_pool_r1/config_snapshot.yaml \
      --checkpoint runs/prophage_v2_pool_r1/checkpoints/best.pt \
      --sizes 700,1000,1500,2000,3000,5000,8000 --scales 3000,10000 \
      --replicates 10 --out_dir analysis/spikein
"""

from __future__ import annotations

import argparse
import json
import random
import sys
from pathlib import Path

import pandas as pd
import torch
from pyfaidx import Fasta

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.predict_genome import encode_window  # noqa: E402
from src.train_prophage import build_model, load_cfg  # noqa: E402


def read_bed_positives(bed_path: Path):
    out = []
    for line in open(bed_path):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 4 or parts[3] != "positive":
            continue
        out.append((parts[0], int(parts[1]), int(parts[2])))
    return out


def load_model(config: str, checkpoint: str, device: torch.device):
    cfg = load_cfg(config)
    model = build_model(cfg, device)
    ck = torch.load(checkpoint, map_location="cpu", weights_only=False)
    state = next((ck[k] for k in ("model", "model_state_dict", "state_dict")
                  if isinstance(ck, dict) and k in ck), ck)
    model.load_state_dict(state)
    model.eval()
    return model, cfg


def scan(model, seq: str, scale: int, stride: int, device: torch.device,
         policy: str, batch: int = 32) -> pd.DataFrame:
    """Score `seq` with `scale`-bp windows every `stride` bp."""
    starts = list(range(0, max(1, len(seq) - scale + 1), stride))
    probs = []
    with torch.no_grad():
        for i in range(0, len(starts), batch):
            chunk = [seq[st : st + scale] for st in starts[i : i + batch]]
            chunk = [c for c in chunk if len(c) == scale]
            if not chunk:
                continue
            x = torch.stack([encode_window(c, policy) for c in chunk]).to(device)
            y = model(x)
            if isinstance(y, tuple):
                y = y[0]
            probs.extend(torch.sigmoid(y.reshape(-1).float()).cpu().numpy().tolist())
    starts = starts[: len(probs)]
    return pd.DataFrame({"start": starts, "end": [s + scale for s in starts], "prob": probs})


def call_regions(track: pd.DataFrame, threshold: float) -> list:
    """Contiguous runs of windows at or above `threshold`, as [start, end) spans."""
    hits = track[track.prob >= threshold]
    if hits.empty:
        return []
    spans, cur = [], None
    for r in hits.itertuples():
        if cur is None:
            cur = [r.start, r.end]
        elif r.start <= cur[1]:
            cur[1] = max(cur[1], r.end)
        else:
            spans.append(tuple(cur))
            cur = [r.start, r.end]
    spans.append(tuple(cur))
    return spans


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--host_fasta", required=True)
    ap.add_argument("--host_contig", required=True)
    ap.add_argument("--host_start", type=int, required=True)
    ap.add_argument("--host_end", type=int, required=True)
    ap.add_argument("--element_fasta", help="any element sequence; omit to carve from --donor_bed")
    ap.add_argument("--donor_bed", help="BED of positives in --host_fasta (beds_pos/, HPCC-only)")
    ap.add_argument("--donor_tsv", help="geNomad *_provirus.tsv for --host_fasta; the checked-in "
                                        "eval_genomad_truth*/ dirs hold these, so this works locally")
    ap.add_argument("--config", required=True)
    ap.add_argument("--checkpoint", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--sizes", default="700,1000,1500,2000,3000,5000,8000")
    ap.add_argument("--scales", default="3000,10000")
    ap.add_argument("--stride", type=int, default=250)
    ap.add_argument("--replicates", type=int, default=10, help="distinct implant sequences per size")
    ap.add_argument("--whole_elements", action="store_true",
                    help="implant each --element_fasta record WHOLE at its own length, ignoring "
                         "--sizes. Tests whether a complete element detects better than a random "
                         "same-length fragment of a larger one -- the main caveat on the size curve.")
    ap.add_argument("--threshold", type=float, default=0.995)
    ap.add_argument("--max_host_prob", type=float, default=0.5,
                    help="reject the host slice if any window scores above this")
    ap.add_argument("--ambiguity_policy", default="n-channel")
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    device = torch.device(args.device)
    model, _ = load_model(args.config, args.checkpoint, device)
    rng = random.Random(args.seed)
    sizes = [int(x) for x in args.sizes.split(",") if x.strip()]
    scales = [int(x) for x in args.scales.split(",") if x.strip()]

    fa = Fasta(args.host_fasta, as_raw=True, sequence_always_upper=True)
    try:
        host = str(fa[args.host_contig][args.host_start : args.host_end])
        # --- host verification: a "quiet" host is the whole basis of the measurement
        ctrl = scan(model, host, min(scales), 1000, device, args.ambiguity_policy)
        if ctrl.empty:
            raise SystemExit("host slice too short to scan")
        host_max, host_med = float(ctrl.prob.max()), float(ctrl.prob.median())
        print(f"[host] {len(host)} bp, window-score median {host_med:.4f}, max {host_max:.4f}", flush=True)
        if host_max > args.max_host_prob:
            raise SystemExit(
                f"host slice is not quiet (max {host_max:.4f} > {args.max_host_prob}); "
                "it likely overlaps a real element -- pick a different interval")

        # --- implant source
        donors = []
        if args.element_fasta:
            efa = Fasta(args.element_fasta, as_raw=True, sequence_always_upper=True)
            try:
                donors = [str(efa[k][:]) for k in efa.keys()]
            finally:
                efa.close()
        elif args.donor_bed:
            for contig, s, e in read_bed_positives(Path(args.donor_bed)):
                if contig in fa and (e - s) >= max(sizes):
                    donors.append(str(fa[contig][s:e]))
        elif args.donor_tsv:
            t = pd.read_csv(args.donor_tsv, sep="\t")
            for r in t.itertuples():
                contig = str(r.source_seq)
                if contig in fa and (int(r.end) - int(r.start)) >= max(sizes):
                    donors.append(str(fa[contig][int(r.start) : int(r.end)]))
        if not donors:
            raise SystemExit("no implant source: pass --element_fasta, --donor_bed or --donor_tsv")
        print(f"[donor] {len(donors)} source sequences, longest {max(len(d) for d in donors)} bp", flush=True)

        rows = []
        if args.whole_elements:
            for i, el in enumerate(donors):
                implant = el
                size = len(implant)
                pos = len(host) // 2
                chim = host[:pos] + implant + host[pos:]
                true_s, true_e = pos, pos + size
                for scale in scales:
                    if len(chim) < scale:
                        continue
                    track = scan(model, chim, scale, args.stride, device, args.ambiguity_policy)
                    over = track[(track.end > true_s) & (track.start < true_e)]
                    spans = call_regions(track, args.threshold)
                    hit = [sp for sp in spans if sp[1] > true_s and sp[0] < true_e]
                    if hit:
                        best = max(hit, key=lambda sp: min(sp[1], true_e) - max(sp[0], true_s))
                        inter = min(best[1], true_e) - max(best[0], true_s)
                        union = max(best[1], true_e) - min(best[0], true_s)
                        iou = inter / union if union > 0 else 0.0
                        called_len = best[1] - best[0]
                    else:
                        iou, called_len = 0.0, 0
                    rows.append(dict(size_bp=size, scale_bp=scale, replicate=i,
                                     max_prob=float(over.prob.max()) if len(over) else 0.0,
                                     n_windows_over=int((over.prob >= args.threshold).sum()),
                                     detected=bool(hit), iou=iou, called_len_bp=called_len,
                                     n_calls_overlapping=len(hit), n_calls_total=len(spans)))
                if (i + 1) % 5 == 0:
                    print(f"[whole] {i + 1}/{len(donors)} elements", flush=True)
            sizes = []  # skip the size loop below
        for size in sizes:
            usable = [d for d in donors if len(d) >= size]
            if not usable:
                print(f"[skip] no donor >= {size} bp", flush=True)
                continue
            for rep in range(args.replicates):
                d = usable[rep % len(usable)]
                off = rng.randrange(0, len(d) - size + 1)
                implant = d[off : off + size]
                pos = len(host) // 2
                chim = host[:pos] + implant + host[pos:]
                true_s, true_e = pos, pos + size
                for scale in scales:
                    if len(chim) < scale:
                        continue
                    track = scan(model, chim, scale, args.stride, device, args.ambiguity_policy)
                    over = track[(track.end > true_s) & (track.start < true_e)]
                    spans = call_regions(track, args.threshold)
                    hit = [sp for sp in spans if sp[1] > true_s and sp[0] < true_e]
                    if hit:
                        # best-overlapping call
                        best = max(hit, key=lambda sp: min(sp[1], true_e) - max(sp[0], true_s))
                        inter = min(best[1], true_e) - max(best[0], true_s)
                        union = max(best[1], true_e) - min(best[0], true_s)
                        iou = inter / union if union > 0 else 0.0
                        called_len = best[1] - best[0]
                    else:
                        iou, called_len = 0.0, 0
                    rows.append(dict(size_bp=size, scale_bp=scale, replicate=rep,
                                     max_prob=float(over.prob.max()) if len(over) else 0.0,
                                     n_windows_over=int((over.prob >= args.threshold).sum()),
                                     detected=bool(hit), iou=iou, called_len_bp=called_len,
                                     n_calls_overlapping=len(hit), n_calls_total=len(spans)))
            print(f"[size] {size} bp done", flush=True)
    finally:
        fa.close()

    df = pd.DataFrame(rows)
    df.to_csv(out / "spikein_raw.csv", index=False)
    if args.whole_elements and not df.empty:
        df["size_bin"] = pd.cut(df.size_bp, [0, 6000, 9000, 12000, 10**9],
                                labels=["<6kb", "6-9kb", "9-12kb", "12kb+"])
        binned = (df.groupby(["scale_bp", "size_bin"], observed=True)
                    .agg(n=("replicate", "size"), detect_rate=("detected", "mean"),
                         median_max_prob=("max_prob", "median"), median_iou=("iou", "median"))
                    .reset_index())
        print("\nWHOLE complete elements, binned by true size:", flush=True)
        print(binned.to_string(index=False), flush=True)
    summary = (df.groupby(["scale_bp", "size_bp"])
                 .agg(n=("replicate", "size"), detect_rate=("detected", "mean"),
                      median_max_prob=("max_prob", "median"), median_iou=("iou", "median"),
                      median_called_len=("called_len_bp", "median"))
                 .reset_index())
    res = {"host": {"fasta": args.host_fasta, "contig": args.host_contig,
                    "start": args.host_start, "end": args.host_end,
                    "median_window_prob": host_med, "max_window_prob": host_max},
           "params": {"threshold": args.threshold, "stride": args.stride,
                      "ambiguity_policy": args.ambiguity_policy, "replicates": args.replicates},
           "by_scale_and_size": summary.to_dict("records"),
           "note": ("detect_rate is detection; median_iou is delimitation. They are different "
                    "limits -- a scan can fire on an element far smaller than it can delimit. "
                    "Implant context is artificial; this bounds pipeline resolution, not "
                    "real-world sensitivity to a given element class.")}
    (out / "results.json").write_text(json.dumps(res, indent=2) + "\n")
    pd.set_option("display.width", 200)
    print("\n" + summary.round(4).to_string(index=False), flush=True)
    print(f"\n[done] wrote {out/'results.json'} and {out/'spikein_raw.csv'}", flush=True)


if __name__ == "__main__":
    main()
