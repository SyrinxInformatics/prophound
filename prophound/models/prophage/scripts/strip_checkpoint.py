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
"""Strip training checkpoints to weights, for bundling or for `predict` without --trust-checkpoint.

A training checkpoint also holds the optimizer, scaler, RNG states, history and config, so it
can only be read by unpickling arbitrary objects. This writes `{"model": state_dict}` files that
`torch.load(weights_only=True)` reads, a model-only config (`model.yaml`), and `manifest.json`
with every file's SHA-256 and the SHA-256 of the checkpoint it came from. Each stripped file is
read back and every tensor compared with the original before anything is reported done.
docs/tech-debt-audit.md S8; the bundled model is made with:

  python scripts/strip_checkpoint.py --name prophage_v2_pool --out src/../model \\
      runs/prophage_v2_pool_r1/checkpoints/best.pt runs/prophage_v2_pool_r2/checkpoints/best.pt \\
      runs/prophage_v2_pool_r3/checkpoints/best.pt
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import torch
import yaml


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("checkpoints", nargs="+", help="training checkpoints (best.pt)")
    ap.add_argument("--name", required=True, help="model name; files are <name>_r<k>.pt")
    ap.add_argument("--out", required=True, help="directory for the stripped files, model.yaml and manifest.json")
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    entries, model_cfg, data_cfg = [], None, None
    for k, src in enumerate(map(Path, a.checkpoints), 1):
        ck = torch.load(src, map_location="cpu", weights_only=False)
        cfg = ck["cfg"]
        model = {key: v for key, v in cfg["model"].items() if key != "predict_boundary"}
        data = {key: cfg["data"][key] for key in ("window_len", "stride")}
        if model_cfg is None:
            model_cfg, data_cfg = model, data
        elif (model, data) != (model_cfg, data_cfg):
            raise SystemExit(f"{src}: model or window settings differ from the first checkpoint's")
        state = {key: v.detach().cpu().contiguous() for key, v in ck["model"].items()}
        dest = out / f"{a.name}_r{k}.pt"
        torch.save({"model": state}, dest)
        back = torch.load(dest, map_location="cpu", weights_only=True)["model"]
        if list(back) != list(state) or not all(torch.equal(back[key], state[key]) for key in state):
            raise SystemExit(f"{dest}: tensors read back differ from {src}")
        entries.append({"file": dest.name, "sha256": sha256(dest), "source": str(src), "source_sha256": sha256(src),
                        "source_run_name": cfg["run_name"], "tensors": len(state),
                        "parameters": int(sum(v.numel() for v in state.values()))})
        print(f"[strip] {src} -> {dest.name}: {len(state)} tensors identical, {dest.stat().st_size / 1e6:.1f} MB")

    config = {"run_name": a.name, "data": data_cfg, "model": model_cfg}
    (out / "model.yaml").write_text(
        f"# {a.name}: the bundled model's config (window settings and architecture only;\n"
        "# training paths are not needed for inference). Written by scripts/strip_checkpoint.py.\n"
        + yaml.safe_dump(config, sort_keys=False))
    manifest = {"name": a.name, "config": "model.yaml", "config_sha256": sha256(out / "model.yaml"),
                "checkpoints": entries}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"[done] {len(entries)} checkpoint(s), model.yaml and manifest.json in {out}")


if __name__ == "__main__":
    main()
