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
"""The canonical smoke test: the bundled model on the example genome makes the expected calls.

Slow (about 2 minutes on CPU), so it runs only with `-m slow` (`pixi run smoke`). Calls are
compared, not scores: devices agree on calls, not on every last digit (docs/benchmarks.md).
"""

import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"
CALL_COLUMNS = ["call_id", "contig_id", "start", "end", "boundary", "contig_edge", "seed_support"]


@pytest.mark.slow
def test_example_genome_gives_the_expected_calls(tmp_path):
    out = tmp_path / "out"
    proc = subprocess.run([sys.executable, "-c", "from prophound.cli import main; main()", "predict",
                           "--input", str(EXAMPLES / "hin_rd_kw20_region.fna"), "--outdir", str(out), "-q"],
                          capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    got = pd.read_csv(out / "regions.tsv", sep="\t")[CALL_COLUMNS]
    expected = pd.read_csv(EXAMPLES / "expected_calls.tsv", sep="\t")[CALL_COLUMNS]
    pd.testing.assert_frame_equal(got, expected)
