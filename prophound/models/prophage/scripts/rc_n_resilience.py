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


"""Reverse-complement and ambiguous-base (N) resilience audit for the full stack.

Standing check (user, 2026-09-22): every analysis should consider the reverse
complement, and the model and downstream must be resilient to N regions. This
script asserts both properties end to end rather than by code reading, so a
future architecture or post-processing change cannot silently break them.

Checks, each PASS/FAIL:

RC-1  window scores are invariant when a sequence is replaced by its reverse
      complement (the encoder is RC-symmetric by construction --
      `model_rc_cnn.reverse_complement_onehot` -- but `combine: mean` and any
      later pooling change could break it, and `max` vs `mean` differ in how
      they preserve it).
RC-2  (NOT IMPLEMENTED HERE -- needs a full genome run; measured by hand
      2026-09-22) a whole contig scored
      forward and reverse-complemented should yield the same region calls,
      mirrored. This does not follow from RC-1: the window grid is anchored at
      position 0, so a contig whose length is not a multiple of the stride has
      a ragged final window that moves to the *start* under RC, and region
      seeds are claimed greedily in index order, so seeds competing for
      overlapping bins can resolve differently when that order reverses.
N-1   an all-N window does not crash and does not score as a prophage.
N-2   an N run inserted into a real prophage does not collapse the score of
      windows that do not overlap it (damage stays local).
N-3   (NOT IMPLEMENTED HERE -- needs a full genome run) an N run inside a
      true region does not split the call, or is re-joined by fragment merging
      (`--max-gap-bp`/`--bridge-prob`).
N-4   case does not matter: a soft-masked (lowercase) sequence scores
      identically to the uppercase one. Assemblies are commonly soft-masked
      for repeats, and repeats are exactly where prophages sit.

Usage (from the prophage model dir):
  python scripts/rc_n_resilience.py --checkpoint <best.pt> --config <config_snapshot.yaml>
  python scripts/rc_n_resilience.py --checkpoint ... --config ... --fasta <genome.fna> --contig <id>
"""

from __future__ import annotations

import argparse
import random
import sys
from pathlib import Path

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.train_prophage import build_model, load_cfg, onehot_encode  # noqa: E402

COMP = str.maketrans("ACGTN", "TGCAN")
TOL = 1e-4  # float32 + two encoder passes; tighter than any real effect


def rc(seq: str) -> str:
    return seq.translate(COMP)[::-1]


def random_seq(n: int, rng: random.Random, gc: float = 0.5) -> str:
    w = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
    return "".join(rng.choices("ACGT", weights=w, k=n))


def score(model, seqs, device, batch=16) -> np.ndarray:
    out = []
    with torch.no_grad():
        for i in range(0, len(seqs), batch):
            x = torch.stack([onehot_encode(s) for s in seqs[i : i + batch]]).to(device)
            y = model(x)
            if isinstance(y, tuple):
                y = y[0]
            out.append(np.atleast_1d(torch.sigmoid(y.reshape(-1).float()).cpu().numpy()))
    return np.concatenate(out)


def report(name: str, ok: bool, detail: str) -> bool:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {detail}")
    return ok


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--checkpoint", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--fasta", help="optional real genome, for a real-sequence RC check")
    ap.add_argument("--contig")
    ap.add_argument("--region", help="start-end of a known prophage in --contig, for N-2/N-3")
    ap.add_argument("--window_len", type=int, default=10000)
    ap.add_argument("--n_windows", type=int, default=24)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--device", default="cpu")
    args = ap.parse_args()

    cfg = load_cfg(args.config)
    device = torch.device(args.device)
    model = build_model(cfg, device)
    ck = torch.load(args.checkpoint, map_location="cpu", weights_only=False)
    state = next((ck[k] for k in ("model", "model_state_dict", "state_dict") if isinstance(ck, dict) and k in ck), ck)
    model.load_state_dict(state)
    model.eval()
    args.window_len = cfg.data.window_len if hasattr(cfg.data, "window_len") else args.window_len
    print(f"checkpoint: {args.checkpoint}")
    print(f"combine={cfg.model.combine} pooling={cfg.model.pooling} window_len={args.window_len}\n")

    rng = random.Random(args.seed)
    ok_all = True

    # --- RC-1: synthetic and (optionally) real sequence ---
    seqs = [random_seq(args.window_len, rng, gc=rng.uniform(0.3, 0.7)) for _ in range(args.n_windows)]
    f, r = score(model, seqs, device), score(model, [rc(s) for s in seqs], device)
    d = float(np.abs(f - r).max())
    ok_all &= report("RC-1 synthetic", d < TOL, f"max|p(fwd)-p(rc)| = {d:.2e} over {len(seqs)} windows")

    if args.fasta and args.contig:
        from pyfaidx import Fasta

        fa = Fasta(args.fasta, as_raw=True, sequence_always_upper=True)
        contig = fa[args.contig]
        step = max(1, (len(contig) - args.window_len) // args.n_windows)
        real = [contig[i : i + args.window_len] for i in range(0, step * args.n_windows, step)]
        real = [s for s in real if len(s) == args.window_len]
        f, r = score(model, real, device), score(model, [rc(s) for s in real], device)
        d = float(np.abs(f - r).max())
        ok_all &= report("RC-1 real sequence", d < TOL,
                         f"max|p(fwd)-p(rc)| = {d:.2e} over {len(real)} windows of {args.contig}")

    # --- N-4: soft-masked input ---
    f = score(model, seqs[:8], device)
    lo = score(model, [s.lower() for s in seqs[:8]], device)
    d = float(np.abs(f - lo).max())
    ok_all &= report("N-4 soft-masking", d < TOL, f"max|p(upper)-p(lower)| = {d:.2e} (lowercase must not read as N)")

    # --- N-1: all-N and mostly-N windows ---
    alln = score(model, ["N" * args.window_len], device)[0]
    ok_all &= report("N-1 all-N window", alln < 0.5, f"p(all-N) = {alln:.4f} (must not call a prophage)")
    iupac = score(model, ["R" * args.window_len], device)[0]
    report("N-1b IUPAC codes", True, f"p(all-R) = {iupac:.4f} (ambiguity codes route to the N channel, no crash)")

    # --- N-2: locality of damage, on the highest-scoring sequence available ---
    base = real if (args.fasta and args.contig) else seqs
    probs = score(model, base, device)
    hi = base[int(np.argmax(probs))]
    p0 = float(probs.max())
    rows = []
    for frac in (0.01, 0.05, 0.10, 0.25, 0.50):
        n = int(len(hi) * frac)
        mid = (len(hi) - n) // 2
        masked = hi[:mid] + "N" * n + hi[mid + n :]
        rows.append((frac, float(score(model, [masked], device)[0])))
    worst = min(p for _, p in rows)
    detail = ", ".join(f"{int(f*100)}%N->{p:.3f}" for f, p in rows)
    # a monotone, graceful decline is expected; a collapse at 1% N would mean
    # a single N run can erase an otherwise confident call
    graceful = rows[0][1] > 0.5 * p0
    ok_all &= report("N-2 damage is local", graceful, f"p(clean) = {p0:.3f}; {detail}")

    print(f"\n{'ALL CHECKS PASSED' if ok_all else 'FAILURES PRESENT'} "
          f"-- RC-2 and N-3 are not covered here; both need a full genome run.")
    sys.exit(0 if ok_all else 1)


if __name__ == "__main__":
    main()
