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
"""Robustness check 1, the shuffled-genome null (docs/robustness-checks.md).

  prepare  check each lockbox-2 FASTA against its 851599 input_sha256, shuffle every contig with
           esl-shuffle -d, and check names, lengths and mono-/dinucleotide counts are unchanged
  run      prophound predict (Default) on each shuffled genome; skips genomes already done
  score    calls per Mb on shuffled sequence against the real-genome comparators

Run from src/prophound/models/prophage:
  python scripts/shuffle_null.py prepare --fasta_dir ../../../../arc_data/lockbox2_fasta --out ../../../../arc_data/lockbox2_shuffled
  python scripts/shuffle_null.py run --shuffled ../../../../arc_data/lockbox2_shuffled --runs analysis/robustness/shuffle_null
  python scripts/shuffle_null.py score --runs analysis/robustness/shuffle_null --out analysis/robustness
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from collections import Counter
from pathlib import Path

import pandas as pd
from scipy.stats import chi2

L2 = Path("analysis/lockbox2")
SEED = 20261028
ESL = ["pixi", "exec", "-c", "conda-forge", "-c", "bioconda", "-s", "hmmer", "--", "esl-shuffle"]


def genomes() -> list[str]:
    return [l.strip() for l in open(L2 / "lockbox_genome_list.txt") if l.strip()]


def read_fasta(p: Path) -> list[tuple[str, str]]:
    out, name, chunks = [], None, []
    for line in open(p):
        if line.startswith(">"):
            if name is not None:
                out.append((name, "".join(chunks)))
            name, chunks = line[1:].split(None, 1)[0], []
        else:
            chunks.append(line.strip().upper())
    if name is not None:
        out.append((name, "".join(chunks)))
    return out


def composition(s: str) -> tuple[Counter, Counter]:
    return Counter(s), Counter(s[i:i + 2] for i in range(len(s) - 1))


def cmd_prepare(a) -> None:
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    report = []
    for g in genomes():
        src = Path(a.fasta_dir) / f"{g}.fna"
        want = json.load(open(L2 / "prophound" / "default" / g / "run_metadata.json"))["input_sha256"]
        got = hashlib.sha256(src.read_bytes()).hexdigest() if src.exists() else None
        row = dict(genome=g, sha256_ok=got == want)
        if row["sha256_ok"]:
            dst, tmp = out / f"{g}.fna", out / f"{g}.esl.tmp"
            subprocess.run(ESL + ["-d", "--seed", str(SEED), "-o", str(tmp), str(src)], check=True, capture_output=True)
            # esl-shuffle names each contig "<name>-shuffled" and drops the description; put the original header back
            headers = [l.rstrip("\n") for l in open(src) if l.startswith(">")]
            with open(dst, "w") as f:
                for h, (_, s) in zip(headers, read_fasta(tmp)):
                    f.write(h + "\n" + "\n".join(s[i:i + 80] for i in range(0, len(s), 80)) + "\n")
            tmp.unlink()
            orig, shuf = read_fasta(src), read_fasta(dst)
            row["contigs"] = len(orig)
            row["composition_ok"] = len(orig) == len(shuf) and all(
                n1 == n2 and len(s1) == len(s2) and composition(s1) == composition(s2) and s1 != s2
                for (n1, s1), (n2, s2) in zip(orig, shuf))
            row["non_acgt"] = sum(len(s) - sum(s.count(b) for b in "ACGT") for _, s in orig)
        report.append(row)
        print(row, flush=True)
    t = pd.DataFrame(report); t.to_csv(out / "prepare_report.tsv", sep="\t", index=False)
    print(f"\nsha256 ok {int(t.sha256_ok.sum())}/{len(t)}; composition ok {int(t.composition_ok.fillna(False).sum())}")


def cmd_run(a) -> None:
    runs = Path(a.runs); runs.mkdir(parents=True, exist_ok=True)
    report = pd.read_csv(Path(a.shuffled) / "prepare_report.tsv", sep="\t")
    ok = report[report.sha256_ok & report.composition_ok.fillna(False).astype(bool)].genome
    for g in ok:
        if (runs / g / "run_metadata.json").exists():
            continue
        with open(runs / f"{g}.log", "w") as log:
            subprocess.run(["prophound", "predict", "--input", str(Path(a.shuffled) / f"{g}.fna"),
                            "--outdir", str(runs / g), "--force", "-q"], check=True, stdout=log, stderr=subprocess.STDOUT)
        print(g, flush=True)


def poisson_ci(k: int, t: float) -> tuple[float, float]:
    lo = chi2.ppf(0.025, 2 * k) / 2 / t if k else 0.0
    return lo, chi2.ppf(0.975, 2 * k + 2) / 2 / t


def cmd_score(a) -> None:
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    truth = pd.read_csv(L2 / "lockbox_elements.csv")
    xc = pd.read_csv("analysis/extra_calls/regions.tsv", sep="\t")
    xc = xc[xc.kind == "candidate_default"]
    rows = []
    for d in sorted(Path(a.runs).iterdir()):
        if not (d / "run_metadata.json").exists():
            continue
        g = d.name; m = json.load(open(d / "run_metadata.json")); bp = m["input"]["total_bp"]
        c = pd.read_csv(d / "regions.detailed.tsv", sep="\t")
        real = pd.read_csv(L2 / "prophound" / "default" / g / "regions.detailed.tsv", sep="\t")
        seq = "".join(s for _, s in read_fasta(Path(a.shuffled) / f"{g}.fna")) if a.shuffled else ""
        rows.append(dict(genome=g, mb=bp / 1e6, gc=(seq.count("G") + seq.count("C")) / len(seq) if seq else None,
                         element_free=g not in set(truth.accession), calls=len(c), bp_called=int((c.end - c.start).sum()),
                         real_calls=len(real), real_bp_called=int((real.end - real.start).sum()), real_offgenomad=int((xc.genome == g).sum()),
                         real_nosignal=int(((xc.genome == g) & (xc.bin == 4)).sum()),
                         call_kb=";".join(f"{(e - s) / 1000:.1f}" for s, e in zip(c.start, c.end)),
                         support=";".join(map(str, c.seed_support)) if "seed_support" in c else ""))
    t = pd.DataFrame(rows); t.to_csv(out / "shuffle_null.tsv", sep="\t", index=False)
    mb = t.mb.sum(); k = int(t.calls.sum())
    rate, (lo, hi) = k / mb, poisson_ci(k, mb)
    R = t.real_nosignal.sum() / mb
    free = t[t.element_free]
    res = dict(genomes=len(t), mb=round(mb, 1), calls=k, calls_per_mb=rate, ci95=[lo, hi],
               calls_per_genome_median=float(t.calls.median()), calls_per_genome_max=int(t.calls.max()),
               genomes_with_a_call=int((t.calls > 0).sum()), share_bp_called=float(t.bp_called.sum() / (mb * 1e6)),
               # added 2026-10-28 after the interim look (post hoc): the bp comparator the preregistered share lacked
               real_share_bp_called=float(t.real_bp_called.sum() / (mb * 1e6)),
               genomes_over_half_called=int((t.bp_called / (t.mb * 1e6) >= 0.5).sum()),
               real_calls_per_mb=float(t.real_calls.sum() / mb), real_offgenomad_per_mb=float(t.real_offgenomad.sum() / mb),
               real_nosignal_per_mb_R=float(R), ratio_to_R=rate / R if R else None,
               element_free_real_per_mb=float(free.real_calls.sum() / free.mb.sum()) if len(free) else None,
               reading="under 0.1 R" if rate < 0.1 * R else "at least 0.5 R" if rate >= 0.5 * R else "in between")
    if t.gc.notna().any():
        res["spearman_calls_per_mb_vs_gc"] = float((t.calls / t.mb).corr(t.gc, method="spearman"))
    json.dump(res, open(out / "shuffle_null.json", "w"), indent=1)
    print(json.dumps(res, indent=1))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("prepare"); p.add_argument("--fasta_dir", required=True); p.add_argument("--out", required=True)
    r = sub.add_parser("run"); r.add_argument("--shuffled", required=True); r.add_argument("--runs", required=True)
    s = sub.add_parser("score"); s.add_argument("--runs", required=True); s.add_argument("--out", required=True)
    s.add_argument("--shuffled", default=None)
    a = ap.parse_args()
    dict(prepare=cmd_prepare, run=cmd_run, score=cmd_score)[a.cmd](a)


if __name__ == "__main__":
    main()
