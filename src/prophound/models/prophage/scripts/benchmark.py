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
"""Wall-clock and memory benchmark of `prophound predict` (docs/benchmarks.md).

  timing   time the shipped CLI end to end: per input, one warm-up run, then
           --repeats timed runs; wall seconds, peak RSS, and the regions.tsv hash
           (repeat runs must agree)
  profile  split one input's scoring into window reading + encoding (CPU) and the
           model forward pass (device), in process, with the code from --code

The machine, OS, power source, torch build and code commit are recorded with every
result. Run on a quiet machine, on mains power, from the repo root:
  python src/prophound/models/prophage/scripts/benchmark.py timing --code . --device mps \\
      --input typical=path.fna --out analysis/benchmarks/mps_timing.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import re
import statistics
import subprocess
import sys
import time
from pathlib import Path

MODEL_DIR = Path(__file__).resolve().parents[1]
CHECKPOINTS = [MODEL_DIR / f"runs/prophage_v2_pool_r{k}/checkpoints/best.pt" for k in (1, 2, 3)]
CONFIG = MODEL_DIR / "runs/prophage_v2_pool_r1/config_snapshot.yaml"


def run_text(cmd: list[str]) -> str:
    try:
        return subprocess.run(cmd, capture_output=True, text=True, check=False).stdout.strip()
    except OSError:
        return ""


def machine_info(code: Path, env: dict) -> dict:
    probe = subprocess.run(
        [sys.executable, "-c", "import torch, json; print(json.dumps({'torch': torch.__version__, "
         "'threads': torch.get_num_threads(), 'mps': torch.backends.mps.is_available(), "
         "'cuda': torch.cuda.is_available(), 'mkldnn': torch.backends.mkldnn.is_available()}))"],
        env=env, capture_output=True, text=True, check=True)
    power = run_text(["pmset", "-g", "batt"]).splitlines()
    return {
        "machine": platform.machine(), "python": platform.python_version(),
        "os": f"{run_text(['sw_vers', '-productName'])} {run_text(['sw_vers', '-productVersion'])}".strip(),
        "model": run_text(["sysctl", "-n", "hw.model"]), "chip": run_text(["sysctl", "-n", "machdep.cpu.brand_string"]),
        "cores": run_text(["sysctl", "-n", "hw.ncpu"]),
        "performance_cores": run_text(["sysctl", "-n", "hw.perflevel0.physicalcpu"]),
        "efficiency_cores": run_text(["sysctl", "-n", "hw.perflevel1.physicalcpu"]),
        "memory_gb": round(int(run_text(["sysctl", "-n", "hw.memsize"]) or 0) / 2**30, 1),
        "power": power[0] if power else "", "battery": power[1].strip() if len(power) > 1 else "",
        "omp_num_threads": env.get("OMP_NUM_THREADS"),
        "code_commit": run_text(["git", "-C", str(code), "rev-parse", "HEAD"]),
        "code_dirty": bool(run_text(["git", "-C", str(code), "status", "--porcelain", "--", "src"])),
        **json.loads(probe.stdout),
    }


def fasta_bp(path: Path) -> tuple[int, int]:
    n_contigs = bp = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                n_contigs += 1
            else:
                bp += len(line.strip())
    return n_contigs, bp


def timed_run(code_env: dict, fasta: Path, outdir: Path, device: str, mode: str) -> dict:
    cmd = ["/usr/bin/time", "-l", sys.executable, "-c", "from prophound.cli import main; main()", "predict",
           "--input", str(fasta), "--outdir", str(outdir), "--config", str(CONFIG),
           "--checkpoint", *map(str, CHECKPOINTS), "--device", device, "--mode", mode]
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, env=code_env, capture_output=True, text=True)
    wall = time.perf_counter() - t0
    if proc.returncode:
        raise SystemExit(f"run failed ({fasta}):\n{proc.stderr[-2000:]}")
    rss = re.search(r"(\d+)\s+maximum resident set size", proc.stderr)
    meta = json.loads((outdir / "run_metadata.json").read_text())
    return {"wall_s": round(wall, 2), "peak_rss_mb": round(int(rss.group(1)) / 2**20, 1) if rss else None,
            "n_windows": meta["n_windows"], "n_regions": meta["n_regions"],
            "regions_sha256": hashlib.sha256((outdir / "regions.tsv").read_bytes()).hexdigest()}


def cmd_timing(args: argparse.Namespace) -> None:
    code = Path(args.code).resolve()
    env = {**os.environ, "PYTHONPATH": str(code / "src")}
    if args.threads:
        env["OMP_NUM_THREADS"] = str(args.threads)
    work = Path(args.workdir).resolve()
    result = {"machine": machine_info(code, env), "device": args.device, "mode": args.mode,
              "repeats": args.repeats, "warmup": not args.no_warmup, "inputs": {}}
    for spec in args.input:
        label, _, path = spec.partition("=")
        fasta = Path(path).resolve()
        n_contigs, bp = fasta_bp(fasta)
        runs = []
        for k in range(args.repeats + (0 if args.no_warmup else 1)):
            r = timed_run(env, fasta, work / f"{label}_{k}", args.device, args.mode)
            r["warmup"] = k == 0 and not args.no_warmup
            runs.append(r)
            print(f"[timing] {label} run {k}{' (warm-up)' if r['warmup'] else ''}: {r['wall_s']} s, "
                  f"{r['peak_rss_mb']} MB", flush=True)
        timed = [r for r in runs if not r["warmup"]]
        walls = [r["wall_s"] for r in timed]
        result["inputs"][label] = {
            "path": str(fasta), "n_contigs": n_contigs, "bp": bp, "runs": runs,
            "median_s": statistics.median(walls), "min_s": min(walls), "max_s": max(walls),
            "median_s_per_mb": round(statistics.median(walls) / (bp / 1e6), 2),
            "median_peak_rss_mb": statistics.median(r["peak_rss_mb"] for r in timed if r["peak_rss_mb"]),
            "regions_identical": len({r["regions_sha256"] for r in runs}) == 1,
        }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2) + "\n")
    print(f"[done] wrote {out}")


def sync(device) -> None:
    import torch
    if device.type == "mps":
        torch.mps.synchronize()
    elif device.type == "cuda":
        torch.cuda.synchronize()


def cmd_profile(args: argparse.Namespace) -> None:
    """Time one input's coarse pass stage by stage, with the code from --code."""
    code = Path(args.code).resolve()
    sys.path.insert(0, str(code / "src"))
    import torch
    from torch.utils.data import DataLoader

    from prophound.models.prophage.src import predict_genome as pg

    cfg = pg.load_cfg(str(CONFIG))
    device = pg.resolve_device(args.device)
    fasta = str(Path(args.input).resolve())
    windows = pg.build_windows(fasta, Path(fasta).stem, cfg.data)
    dataset = pg.GenomeWindowDataset(windows, window_len=cfg.data.window_len)
    loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=False, collate_fn=pg.collate_batch)
    model = pg.build_model(cfg, device).eval()
    t0 = time.perf_counter()
    batches = [x for x, _ in loader]
    t_encode = time.perf_counter() - t0
    t_load = t_forward = t_transfer = 0.0
    with torch.no_grad():
        for ckpt in CHECKPOINTS:
            t0 = time.perf_counter()
            model.load_state_dict(torch.load(ckpt, map_location=device, weights_only=False)["model"])
            sync(device)
            t_load += time.perf_counter() - t0
            for x in batches:
                t0 = time.perf_counter()
                xd = x.to(device)
                sync(device)
                t_transfer += time.perf_counter() - t0
                t0 = time.perf_counter()
                model(xd)
                sync(device)
                t_forward += time.perf_counter() - t0
    dataset.close()
    n = len(windows)
    # the CLI encodes once per checkpoint pass, so its encoding cost is 3x the one pass timed here
    result = {"machine": machine_info(code, os.environ.copy()), "device": str(device), "input": fasta,
              "n_windows": n, "batch_size": args.batch_size,
              "encode_one_pass_s": round(t_encode, 2), "encode_as_cli_s": round(3 * t_encode, 2),
              "checkpoint_load_s": round(t_load, 2), "transfer_s": round(t_transfer, 2),
              "forward_s": round(t_forward, 2),
              "ms_per_window": {"encode": round(1e3 * t_encode / n, 3),
                                "forward_per_checkpoint": round(1e3 * t_forward / n / len(CHECKPOINTS), 3)}}
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "machine"}, indent=2))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("timing")
    s.add_argument("--code", default=".", help="source tree (repo root) whose src/ is imported")
    s.add_argument("--device", default="mps")
    s.add_argument("--mode", default="default")
    s.add_argument("--input", action="append", required=True, metavar="LABEL=FASTA")
    s.add_argument("--repeats", type=int, default=3)
    s.add_argument("--no-warmup", action="store_true")
    s.add_argument("--threads", type=int, default=0, help="OMP_NUM_THREADS (0: torch's default)")
    s.add_argument("--workdir", required=True, help="scratch directory for the runs' outputs")
    s.add_argument("--out", required=True)
    s = sub.add_parser("profile")
    s.add_argument("--code", default=".")
    s.add_argument("--device", default="mps")
    s.add_argument("--input", required=True)
    s.add_argument("--batch-size", type=int, default=128)
    s.add_argument("--out", required=True)
    args = ap.parse_args()
    {"timing": cmd_timing, "profile": cmd_profile}[args.cmd](args)


if __name__ == "__main__":
    main()
