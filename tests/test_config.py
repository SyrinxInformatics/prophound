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
"""Run configs: defaults come from the dataclasses, unknown keys are errors, and
inference needs no training paths."""

from pathlib import Path

import pytest
import yaml

from prophound.models.prophage.src.config import load_cfg

CONFIGS = Path(__file__).resolve().parents[1] / "prophound/models/prophage/src/configs"


def write(tmp_path, cfg):
    p = tmp_path / "cfg.yaml"
    p.write_text(yaml.safe_dump(cfg))
    return str(p)


def test_model_only_config_loads_with_defaults(tmp_path):
    cfg = load_cfg(write(tmp_path, {"data": {"window_len": 10_000, "stride": 1_000},
                                    "model": {"pooling": "adaptive_max", "pool_bins": 16}}))
    assert cfg.data.manifest_csv is None and cfg.model.pool_bins == 16
    assert cfg.training.epochs == 20 and cfg.run_name == "prophage_v0"


def test_scalars_are_coerced(tmp_path):
    cfg = load_cfg(write(tmp_path, {"data": {}, "training": {"lr": "1e-4", "early_stopping": {"patience": "3"}}}))
    assert cfg.training.lr == 1e-4 and cfg.training.early_stopping.patience == 3


@pytest.mark.parametrize("cfg, message", [
    ({"data": {"windowlen": 1}}, "unknown key"),
    ({"data": {}, "model": {"predict_boundary": True}}, "no longer supported"),
    ({"model": {}}, "no 'data' section"),
])
def test_bad_configs_fail_clearly(tmp_path, cfg, message):
    with pytest.raises(ValueError, match=message):
        load_cfg(write(tmp_path, cfg))


def test_every_checked_in_config_loads():
    paths = sorted(CONFIGS.glob("*.yaml"))
    assert paths
    for p in paths:
        assert load_cfg(str(p)).data.manifest_csv
