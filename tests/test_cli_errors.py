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
"""Every expected failure of `predict` ends in one `error:` line, without a traceback."""

import argparse

import pytest

from prophound.models.prophage.src import predict_genome
from prophound.models.prophage.src.options import PredictError


def run_raising(monkeypatch, capsys, exc, verbose=False):
    def boom(args):
        raise exc
    monkeypatch.setattr(predict_genome, "run_predict", boom)
    code = predict_genome.run_cli(argparse.Namespace(verbose=verbose, batch_size=128))
    return code, capsys.readouterr().err


@pytest.mark.parametrize("exc, code, expected", [
    (PredictError("config not found: x.yaml"), 1, "error: config not found: x.yaml\n"),
    (KeyboardInterrupt(), 130, "error: interrupted\n"),
    (RuntimeError("MPS backend out of memory (MPS allocated: 9 GB)"), 1,
     "error: the device ran out of memory; try a smaller --batch-size (it was 128)\n"),
    (ValueError("boom"), 1, "error: unexpected ValueError: boom (rerun with --verbose for the traceback)\n"),
])
def test_failures_print_one_error_line(monkeypatch, capsys, exc, code, expected):
    assert run_raising(monkeypatch, capsys, exc) == (code, expected)


def test_verbose_adds_the_traceback(monkeypatch, capsys):
    code, err = run_raising(monkeypatch, capsys, ValueError("boom"), verbose=True)
    assert code == 1 and "Traceback" in err and err.endswith("error: unexpected ValueError: boom\n")


def test_unavailable_device_is_a_user_error():
    from prophound.models.prophage.src.devices import resolve_device, xpu_available
    if xpu_available():
        pytest.skip("this machine has an XPU")
    with pytest.raises(PredictError, match="--device auto"):
        resolve_device("xpu")


def test_input_errors_are_user_errors(tmp_path):
    from prophound.models.prophage.src.io_checks import validate_fasta
    with pytest.raises(PredictError, match="not found"):
        validate_fasta(str(tmp_path / "missing.fna"))


def test_duration_reads_as_text():
    assert predict_genome._duration(43.46) == "43.5 s"
    assert predict_genome._duration(444) == "7 min 24 s"
    assert predict_genome._duration(3 * 3600 + 120) == "3 h 2 min"
