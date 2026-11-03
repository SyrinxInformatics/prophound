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
"""Run configuration: the YAML schema shared by training and inference.

Defaults live only in the dataclasses; `load_cfg` fills them in and rejects keys the
schema does not have, so a typo in a config fails instead of being ignored.
"""

from __future__ import annotations

from dataclasses import MISSING, dataclass, field, fields, is_dataclass
from typing import List, Optional

import yaml


@dataclass
class EarlyStoppingConfig:
    metric: str = "val_auprc"
    patience: int = 5
    mode: str = "max"  # "max" or "min"


@dataclass
class TrainingConfig:
    seed: int = 13
    batch_size: int = 64
    num_workers: int = 2
    epochs: int = 20
    lr: float = 3e-4
    weight_decay: float = 1e-2
    grad_clip: float = 1.0
    amp: bool = True
    save_every_sec: int = 1800
    save_every_steps: int = 0
    slurm_margin_sec: int = 600
    auto_resume: bool = True
    init_checkpoint: Optional[str] = None
    init_mode: str = "strict"
    early_stopping: EarlyStoppingConfig = field(default_factory=EarlyStoppingConfig)


@dataclass
class DataConfig:
    # training inputs; inference reads only window_len and stride
    manifest_csv: Optional[str] = None
    splits_dir: Optional[str] = None
    window_len: int = 10_000
    stride: int = 1_000
    use_neg_types_in_train: List[str] = field(
        default_factory=lambda: ["far_neg", "positive"]
    )
    use_neg_types_in_val: List[str] = field(
        default_factory=lambda: ["far_neg", "near_neg", "positive"]
    )
    use_neg_types_in_test: List[str] = field(
        default_factory=lambda: ["far_neg", "near_neg", "positive"]
    )


@dataclass
class ModelConfig:
    type: str = "rc_cnn"
    combine: str = "max"
    channels: List[int] = field(default_factory=lambda: [64, 128, 256])
    kernel_sizes: List[int] = field(default_factory=lambda: [15, 11, 7])
    dropout: float = 0.2
    fc_dim: int = 256
    pooling: str = "global_max"
    pool_bins: int = 1


@dataclass
class LossConfig:
    pos_weight: str | float = "auto"


@dataclass
class ProphageTrainConfig:
    run_name: str = "prophage_v0"
    data: DataConfig = field(default_factory=DataConfig)
    training: TrainingConfig = field(default_factory=TrainingConfig)
    model: ModelConfig = field(default_factory=ModelConfig)
    loss: LossConfig = field(default_factory=LossConfig)


# scalar fields are coerced as before (YAML reads "1e-4" as a string); others pass through
_COERCE = {"int": int, "float": float, "bool": bool}
# (section, key) -> why it was removed; a true value is an error, a false one is ignored
_REMOVED = {("model", "predict_boundary"): "the boundary head was removed"}


def _build(cls, raw: dict | None, section: str, where: str):
    raw = raw or {}
    if not isinstance(raw, dict):
        raise ValueError(f"{where}: '{section}' must be a mapping, got {type(raw).__name__}")
    by_name = {f.name: f for f in fields(cls)}
    unknown = []
    for key, value in raw.items():
        if (section, key) in _REMOVED:
            if value:
                raise ValueError(f"{where}: {section}.{key} is no longer supported ({_REMOVED[section, key]})")
        elif key not in by_name:
            unknown.append(key)
    if unknown:
        raise ValueError(f"{where}: unknown key(s) in '{section}': {sorted(unknown)}")
    kwargs = {}
    for name, f in by_name.items():
        if name not in raw:
            continue
        default = f.default_factory() if f.default_factory is not MISSING else f.default
        if is_dataclass(default):
            kwargs[name] = _build(type(default), raw[name], name, where)
        else:
            coerce = _COERCE.get(str(f.type))
            kwargs[name] = coerce(raw[name]) if coerce and raw[name] is not None else raw[name]
    return cls(**kwargs)


def load_cfg(path: str) -> ProphageTrainConfig:
    """Read a run config. `data` must be present; `data.manifest_csv` and
    `data.splits_dir` are needed only for training."""
    with open(path) as f:
        raw = yaml.safe_load(f) or {}
    if "data" not in raw:
        raise ValueError(f"{path}: no 'data' section (it sets window_len and stride)")
    return _build(ProphageTrainConfig, raw, "top level", path)
