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
"""Ambiguity-code encoding: the 2026-09-24 default flip to `impute` must be a
no-op on clean sequence and must actually remove ambiguity codes elsewhere."""

import random

import torch

from prophound.models.prophage.src.predict_genome import (
    AMBIGUITY_POLICIES,
    DEFAULT_AMBIGUITY_POLICY,
    encode_window,
)


def test_default_is_impute():
    assert DEFAULT_AMBIGUITY_POLICY == "impute"


def test_all_policies_identical_on_clean_sequence():
    random.seed(0)
    seq = "".join(random.choice("ACGT") for _ in range(3000))
    ref = encode_window(seq, "n-channel")
    for pol in AMBIGUITY_POLICIES:
        assert torch.equal(encode_window(seq, pol), ref), pol


def test_impute_leaves_no_ambiguity_and_is_deterministic():
    random.seed(1)
    seq = list("".join(random.choice("ACGT") for _ in range(3000)))
    for i in (10, 11, 500, 2999):
        seq[i] = "N" if i != 500 else "R"
    seq = "".join(seq)
    a, b = encode_window(seq, "impute"), encode_window(seq, "impute")
    assert torch.equal(a, b)
    assert torch.equal(a, encode_window(seq))            # default path
    n = encode_window(seq, "n-channel")
    assert not torch.equal(a, n)
    if a.shape[0] > 4:                                     # N channel, if the encoding has one
        assert a[4].sum() == 0


def test_vectorised_onehot_matches_the_per_base_definition():
    """The lookup-table encoder (Phase 3) gives exactly the tensors of the per-base loop."""
    import random

    import numpy as np
    import torch

    from prophound.models.prophage.src.encoding import BASE_TO_CH, N_CH, onehot_encode

    def per_base(seq):
        x = np.zeros((5, len(seq)), dtype=np.float32)
        for i, b in enumerate(seq):
            x[BASE_TO_CH.get(b, N_CH), i] = 1.0
        return torch.from_numpy(x)

    rng = random.Random(20261020)
    alphabet = "ACGTacgtNnRYKMSWBDHVuU-*?. é€"
    for n in (0, 1, 7, 3000, 10_000):
        for _ in range(20):
            seq = "".join(rng.choice(alphabet) for _ in range(n))
            assert torch.equal(onehot_encode(seq), per_base(seq))
