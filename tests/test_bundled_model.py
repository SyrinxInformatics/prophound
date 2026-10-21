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
"""The bundled model: weights-only files that match their manifest, load safely, and are
what `predict` uses when no checkpoint is given."""

import shutil

import numpy as np
import pytest
import torch

from prophound.cli import build_parser
from prophound.models.prophage.src import options
from prophound.models.prophage.src import predict_genome as pg
from prophound.models.prophage.src.config import load_cfg
from prophound.models.prophage.src.model_rc_cnn import build_model
from prophound.models.prophage.src.options import PredictError, PredictOptions, bundled_model


def test_bundled_files_match_the_manifest_and_load_weights_only():
    config, checkpoints, manifest = bundled_model()
    assert len(checkpoints) == 3 == len(manifest["checkpoints"])
    pg.verify_bundled_model()
    states = pg.load_checkpoints(build_model(load_cfg(config), torch.device("cpu")), checkpoints,
                                 torch.device("cpu"), config)
    assert len(states) == 3 and all(len(s) == 10 for s in states)


def test_predict_defaults_to_the_bundled_model():
    args = build_parser().parse_args(["predict", "--input", "g.fna", "--outdir", "o"])
    opts = PredictOptions.from_args(args, 10_000, 1_000)
    config, checkpoints, _ = bundled_model()
    assert opts.bundled_model and opts.checkpoints == checkpoints and opts.config == config


def test_a_damaged_bundled_file_is_reported(tmp_path, monkeypatch):
    copy = tmp_path / "model"
    shutil.copytree(options.BUNDLED_MODEL_DIR, copy)
    with open(copy / "prophage_v2_pool_r2.pt", "r+b") as fh:
        fh.seek(1000)
        fh.write(b"\x00\x01\x02")
    monkeypatch.setattr(options, "BUNDLED_MODEL_DIR", copy)
    monkeypatch.setattr(pg, "BUNDLED_MODEL_DIR", copy)
    with pytest.raises(PredictError, match="missing or damaged"):
        pg.verify_bundled_model()


def test_a_full_training_checkpoint_needs_trust(tmp_path):
    config, checkpoints, _ = bundled_model()
    state = torch.load(checkpoints[0], map_location="cpu", weights_only=True)["model"]
    full = tmp_path / "best.pt"
    torch.save({"model": state, "rng": np.random.get_state()}, full)   # numpy state: not weights-only
    model = build_model(load_cfg(config), torch.device("cpu"))
    with pytest.raises(PredictError, match="--trust-checkpoint"):
        pg.load_checkpoints(model, (str(full),), torch.device("cpu"), config)
    assert len(pg.load_checkpoints(model, (str(full),), torch.device("cpu"), config, trust=True)) == 1
