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
"""Reference outputs for behaviour-preserving changes to the shipped CLI.

docs/tech-debt-audit.md Phase 0. A refactor of the hash-checked modules keeps the
lockbox figure only if every output stays the same, so this records what the
current code writes for a matrix of inputs and options, and compares a later
tree against it.

  inputs   build the reference inputs from a local genome (deterministic)
  run      run every case through `prophound predict` from a given source tree
  compare  compare two run directories: output files byte for byte, run_metadata.json
           as JSON, the user-facing stderr lines ([warn], [summary], error:), exit codes.
           With --across-devices, the CPU/MPS/CUDA agreement check instead: the same calls
           and messages, scores within --tolerance

Run from the repo root. CPU is the binding check; MPS is allowed for quick iteration only if two
MPS runs of the same tree are identical. Runs on different devices are never compared:
  python src/prophound/models/prophage/scripts/reference_outputs.py inputs --out analysis_dir/inputs
  python .../reference_outputs.py run --code /path/to/tree --inputs .../inputs --out .../ref_a
  python .../reference_outputs.py compare .../ref_a .../ref_b
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from pyfaidx import Fasta

MODEL_DIR = Path(__file__).resolve().parents[1]
CHECKPOINTS = [MODEL_DIR / f"runs/prophage_v2_pool_r{k}/checkpoints/best.pt" for k in (1, 2, 3)]
CONFIG = MODEL_DIR / "runs/prophage_v2_pool_r1/config_snapshot.yaml"
SOURCE = MODEL_DIR / "analysis/g5/fasta/GCF_000027305.1.fna"   # H. influenzae Rd KW20, NC_000907
SOURCE_CONTIG = "NC_000907"
# Both geNomad elements (1,498,836-1,513,073 and 1,559,961-1,594,275) with flanks.
REGION = (1_440_000, 1_660_000)
BACKGROUND = (0, 150_000)
SHRED_SEED = 20261019
# Shred piece lengths, cycled: under 3 kb (never scanned), 3-10 kb (short-contig scan), 10 kb+.
SHRED_LENGTHS = (2_500, 6_000, 9_000, 4_000, 25_000, 15_000, 3_500, 40_000, 7_500, 12_000)
STDERR_PREFIXES = ("[warn]", "[summary]", "error:")
# volatile parts of a kept line, replaced before comparing
VOLATILE = [(re.compile(r"^\[summary\] ran on \S+ in .*$"), "[summary] ran on <device> in <time>")]


@dataclass(frozen=True)
class Case:
    name: str
    input: str
    args: tuple[str, ...] = ()
    checkpoints: tuple[int, ...] = (1, 2, 3)


CASES = [
    *(Case(f"complete_{m}", "complete", ("--mode", m)) for m in ("default", "exploratory", "conservative")),
    *(Case(f"shred_{m}", "shred", ("--mode", m)) for m in ("default", "exploratory", "conservative")),
    *(Case(f"gappy_{m}", "gappy", ("--mode", m)) for m in ("default", "exploratory", "conservative")),
    Case("complete_outputs", "complete", ("--write-seed-scores", "--write-region-fasta")),
    Case("short_only_default", "short_only"),
    Case("too_short_default", "too_short"),
    Case("complete_refine", "complete",
         ("--refine-window-len", "5000", "--edge-alpha", "0.6", "--edge-alpha-ramp", "15000", "60000")),
    Case("conflated_default", "conflated"),
    Case("conflated_refine", "conflated", ("--refine-window-len", "5000", "--edge-alpha", "0.6")),
    Case("conflated_refine_split", "conflated",
         ("--refine-window-len", "3000", "--split-min-gap-bp", "2000", "--split-max-prob", "0.3")),
    Case("complete_refine_unsupported", "complete", ("--refine-window-len", "5000", "--refine-threshold", "1.0")),
    Case("complete_mean_ensemble", "complete", ("--ensemble", "mean")),
    Case("complete_vote2", "complete", ("--min-seed-support", "2")),
    Case("complete_single_checkpoint", "complete", checkpoints=(1,)),
    Case("complete_thresholds", "complete",
         ("--threshold", "0.9", "--seed-threshold", "0.97", "--extend-threshold", "0.8")),
    Case("complete_builder_opts", "complete", ("--coverage-frac", "0.3", "--max-gap-bp", "0")),
    Case("complete_bridge0", "complete", ("--max-gap-bp", "20000", "--bridge-prob", "0")),
    Case("complete_grid_offset", "complete", ("--grid-offset", "777")),
    Case("complete_window8k", "complete", ("--window-len", "8000", "--stride", "500")),
    Case("gappy_nchannel", "gappy", ("--ambiguity-policy", "n-channel")),
    Case("gappy_uniform", "gappy", ("--ambiguity-policy", "uniform")),
    Case("shred_skip_short", "shred", ("--short-contigs", "skip")),
    Case("shred_short_threshold", "shred", ("--short-contig-threshold", "0.9")),
    Case("shred_no_snap", "shred", ("--edge-snap-frac", "0")),
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")


def gappy(region: str) -> str:
    """The region with two 20 kb N runs and IUPAC codes inside the second element.

    One N run sits in flank, one ends inside the first element, so the N filter
    and imputation over a gap are both exercised. The scattered codes in the
    second element are the EDL933 case (docs/scientific-validity.md "1a"), and the
    NN pair in it is the two-contiguous-codes case.
    """
    seq = list(region)
    for start in (165_000, 50_000):
        seq[start:start + 20_000] = "N" * 20_000
    codes = "RYKMSW"
    for k, pos in enumerate(range(125_000, 150_000, 1_500)):
        seq[pos] = codes[k % len(codes)]
    seq[132_200:132_202] = "NN"
    return "".join(seq)


def shred(region: str) -> list[tuple[str, str]]:
    rng = np.random.default_rng(SHRED_SEED)
    lengths = list(rng.permutation(SHRED_LENGTHS))
    pieces, pos, k = [], 0, 0
    while pos < len(region):
        n = int(lengths[k % len(lengths)])
        pieces.append((f"shred_{k + 1:03d}", region[pos:pos + n]))
        pos += n
        k += 1
    return pieces


def cmd_inputs(args: argparse.Namespace) -> None:
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    fa = Fasta(str(SOURCE), as_raw=True, sequence_always_upper=True)
    region = fa[SOURCE_CONTIG][REGION[0]:REGION[1]]
    background = fa[SOURCE_CONTIG][BACKGROUND[0]:BACKGROUND[1]]
    fa.close()
    inputs = {
        "complete": [("hin_region", region), ("hin_background", background)],
        "gappy": [("hin_region_gappy", gappy(region))],
        "shred": shred(region),
        "short_only": [("short_9500", region[118_000:127_500]), ("short_7000", region[60_000:67_000]),
                       ("short_5200", region[140_000:145_200]), ("short_3100", region[10_000:13_100]),
                       ("short_2800", region[20_000:22_800])],
        # both elements 12 kb apart with background between: one coarse call, split by refinement
        "conflated": [("hin_conflated", region[53_000:73_100] + background[10_000:22_000] + region[119_900:159_000])],
        "too_short": [("tiny_2900", region[118_000:120_900]), ("tiny_1500", region[60_000:61_500])],
    }
    sums = {}
    for name, records in inputs.items():
        path = out / f"{name}.fna"
        write_fasta(path, records)
        sums[path.name] = sha256(path)
        print(f"[inputs] {path.name}: {len(records)} contig(s), {sum(len(s) for _, s in records):,} bp")
    (out / "SHA256SUMS.json").write_text(json.dumps({"source": str(SOURCE), "source_sha256": sha256(SOURCE),
                                                     "files": sums}, indent=2) + "\n")


def user_messages(stderr: str, outdir: Path, exit_code: int) -> list[str]:
    """The lines a user acts on; on failure, also the last line (the exception message)."""
    lines = [ln.replace(str(outdir), "<OUT>") for ln in stderr.splitlines()]
    kept = [ln for ln in lines if ln.startswith(STDERR_PREFIXES)]
    for pattern, repl in VOLATILE:
        kept = [pattern.sub(repl, ln) for ln in kept]
    if exit_code and lines and lines[-1] not in kept:
        kept.append(lines[-1])
    return kept


def cmd_run(args: argparse.Namespace) -> None:
    code, inputs, out = Path(args.code).resolve(), Path(args.inputs).resolve(), Path(args.out).resolve()
    if out.exists() and any(out.iterdir()):
        raise SystemExit(f"{out} is not empty")
    missing = [p for p in [*CHECKPOINTS, CONFIG] if not p.exists()]
    if missing:
        raise SystemExit(f"missing: {', '.join(map(str, missing))}")
    env = {**os.environ, "PYTHONPATH": str(code / "src"), "OMP_NUM_THREADS": str(args.threads)}
    cases = [c for c in CASES if not args.only or c.name in args.only]
    manifest = {"code": str(code), "inputs": json.loads((inputs / "SHA256SUMS.json").read_text()),
                "device": args.device, "threads": args.threads, "cases": {}}
    for case in cases:
        outdir = out / case.name
        cmd = [sys.executable, "-c", "from prophound.cli import main; main()", "predict",
               "--input", str(inputs / f"{case.input}.fna"), "--outdir", str(outdir),
               "--config", str(CONFIG), "--checkpoint", *(str(CHECKPOINTS[k - 1]) for k in case.checkpoints),
               "--device", args.device, *case.args]
        proc = subprocess.run(cmd, env=env, capture_output=True, text=True)
        outdir.mkdir(parents=True, exist_ok=True)
        (outdir / "stderr.txt").write_text(proc.stderr)
        manifest["cases"][case.name] = {"args": list(case.args), "input": case.input,
                                        "checkpoints": list(case.checkpoints), "exit_code": proc.returncode,
                                        "messages": user_messages(proc.stderr, outdir, proc.returncode)}
        print(f"[run] {case.name}: exit {proc.returncode}", flush=True)
    probe = subprocess.run([sys.executable, "-c", "import prophound.cli, torch; print(prophound.cli.__file__); "
                            "print(torch.__version__)"], env=env, capture_output=True, text=True)
    manifest["imported_from"], manifest["torch"] = probe.stdout.split()
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


def compare_case(a: Path, b: Path, ignore_keys: set[str]) -> list[str]:
    diffs = []
    files_a = {p.name for p in a.iterdir() if p.name != "stderr.txt"}
    files_b = {p.name for p in b.iterdir() if p.name != "stderr.txt"}
    if files_a != files_b:
        diffs.append(f"files differ: only in a {sorted(files_a - files_b)}, only in b {sorted(files_b - files_a)}")
    for name in sorted(files_a & files_b):
        if name == "run_metadata.json":
            ma, mb = (json.loads((d / name).read_text()) for d in (a, b))
            for key in sorted((set(ma) | set(mb)) - ignore_keys):
                if ma.get(key) != mb.get(key):
                    diffs.append(f"run_metadata.json[{key}]: {ma.get(key)!r} != {mb.get(key)!r}")
        elif (a / name).read_bytes() != (b / name).read_bytes():
            diffs.append(f"{name} differs")
    return diffs


def compare_case_across_devices(a: Path, b: Path, tol: float) -> list[str]:
    """Same calls and metadata (device aside); scores within `tol`."""
    import pandas as pd

    diffs = []
    for name in ("regions.tsv", "windows.tsv"):
        ta, tb = (pd.read_csv(d / name, sep="\t") for d in (a, b))
        if list(ta.columns) != list(tb.columns) or len(ta) != len(tb):
            diffs.append(f"{name}: shape or columns differ ({ta.shape} vs {tb.shape})")
            continue
        prob_cols = [c for c in ta.columns if "prob" in c]
        exact_cols = [c for c in ta.columns if c not in prob_cols]
        if not ta[exact_cols].equals(tb[exact_cols]):
            bad = [c for c in exact_cols if not ta[c].equals(tb[c])]
            diffs.append(f"{name}: columns differ beyond tolerance: {bad}")
        worst = max(((ta[c] - tb[c]).abs().max() for c in prob_cols), default=0.0)
        # regions.tsv rounds probabilities to 4 dp, so a tiny score difference can round 1e-4 apart
        if worst > tol + (1e-4 if name == "regions.tsv" else 0.0):
            diffs.append(f"{name}: max score difference {worst:.2e} > {tol:g}")
    ma, mb = (json.loads((d / "run_metadata.json").read_text()) for d in (a, b))
    for key in sorted((set(ma) | set(mb)) - {"device"}):
        if ma.get(key) != mb.get(key):
            diffs.append(f"run_metadata.json[{key}] differs")
    return diffs


def cmd_compare(args: argparse.Namespace) -> None:
    a, b = Path(args.a), Path(args.b)
    ma, mb = (json.loads((d / "manifest.json").read_text()) for d in (a, b))
    if ma["inputs"]["files"] != mb["inputs"]["files"]:
        raise SystemExit("the two runs used different inputs")
    if not args.across_devices and (ma["device"], ma["threads"]) != (mb["device"], mb["threads"]):
        raise SystemExit(f"device/threads differ: {ma['device']}/{ma['threads']} vs {mb['device']}/{mb['threads']}")
    ignore = set(args.ignore_metadata_key)
    names = set(ma["cases"]) & set(mb["cases"]) if args.common_only else set(ma["cases"]) | set(mb["cases"])
    n_bad = 0
    for name in sorted(names):
        ca, cb = ma["cases"].get(name), mb["cases"].get(name)
        if ca is None or cb is None:
            print(f"[compare] {name}: only in {'b' if ca is None else 'a'}")
            n_bad += 1
            continue
        diffs = []
        if ca["exit_code"] != cb["exit_code"]:
            diffs.append(f"exit code {ca['exit_code']} != {cb['exit_code']}")
        if ca["messages"] != cb["messages"] and not args.files_only:
            diffs.append(f"messages differ:\n    a: {ca['messages']}\n    b: {cb['messages']}")
        if ca["exit_code"] == 0 and args.across_devices:
            diffs += compare_case_across_devices(a / name, b / name, args.tolerance)
        elif ca["exit_code"] == 0:
            diffs += compare_case(a / name, b / name, ignore)
        n_bad += bool(diffs)
        print(f"[compare] {name}: " + ("identical" if not diffs else "DIFFERS\n  " + "\n  ".join(diffs)))
    print(f"[compare] {len(names) - n_bad} identical, {n_bad} differ"
          + (f" (ignoring metadata keys {sorted(ignore)})" if ignore else ""))
    if n_bad:
        raise SystemExit(1)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("inputs")
    s.add_argument("--out", required=True)
    s = sub.add_parser("run")
    s.add_argument("--code", required=True, help="source tree (repo root) whose src/ is imported")
    s.add_argument("--inputs", required=True)
    s.add_argument("--out", required=True)
    s.add_argument("--threads", type=int, default=4)
    s.add_argument("--device", choices=("cpu", "mps"), default="cpu")
    s.add_argument("--only", nargs="*", default=None, help="case names to run (default: all)")
    s = sub.add_parser("compare")
    s.add_argument("a")
    s.add_argument("b")
    s.add_argument("--ignore-metadata-key", action="append", default=[],
                   help="a run_metadata.json key expected to change (state why in the amendment)")
    s.add_argument("--common-only", action="store_true", help="compare only the cases both runs have")
    s.add_argument("--across-devices", action="store_true",
                   help="CPU/MPS/CUDA agreement: same calls and messages, scores within --tolerance")
    s.add_argument("--tolerance", type=float, default=1e-4)
    s.add_argument("--files-only", action="store_true",
                   help="skip the stderr messages (for a change that rewords them on purpose)")
    args = ap.parse_args()
    {"inputs": cmd_inputs, "run": cmd_run, "compare": cmd_compare}[args.cmd](args)


if __name__ == "__main__":
    main()
