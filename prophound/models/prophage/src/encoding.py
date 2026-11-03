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
"""DNA windows to model input: one-hot encoding, IUPAC ambiguity handling, batching."""

from __future__ import annotations

import numpy as np
import torch

from .options import AMBIGUITY_POLICIES, DEFAULT_AMBIGUITY_POLICY

BASE_TO_CH = {
    "A": 0,
    "C": 1,
    "G": 2,
    "T": 3,
    "a": 0,
    "c": 1,
    "g": 2,
    "t": 3,
}
N_CH = 4  # N channel index
# byte -> channel; every byte not in BASE_TO_CH is N
_CHANNEL_OF_BYTE = np.full(256, N_CH, dtype=np.intp)
for _base, _ch in BASE_TO_CH.items():
    _CHANNEL_OF_BYTE[ord(_base)] = _ch


def onehot_encode(seq: str) -> torch.Tensor:
    """Encode a sequence as a float32 tensor with channels [A, C, G, T, N]; anything
    other than A/C/G/T (either case) sets the N channel."""
    # latin-1 maps each character to one byte; anything it cannot encode becomes "?", i.e. N
    channels = _CHANNEL_OF_BYTE[np.frombuffer(seq.encode("latin-1", errors="replace"), dtype=np.uint8)]
    x = np.zeros((5, len(seq)), dtype=np.float32)
    x[channels, np.arange(len(seq))] = 1.0
    return torch.from_numpy(x)


def encode_window(sequence: str, policy: str = DEFAULT_AMBIGUITY_POLICY) -> torch.Tensor:
    """One-hot encode a window, with IUPAC ambiguity codes (any non-ACGT base) handled per `policy`.

    - `impute` (default): each ambiguous base becomes one concrete A/C/G/T, drawn with a
      seed taken from the window's own content, so a window always encodes the same way.
    - `n-channel`: the training-time encoding; the base sets the N channel. That channel
      is effectively untrained and dominates under max pooling, so two ambiguous bases can
      take a window's score from 1 to 0 (docs/cli-and-packaging.md "Input handling").
    - `uniform`: 0.25 on each of A/C/G/T. Invents no base, but is itself out of distribution.

    All three are exact no-ops on ambiguity-free sequence.
    """
    x = onehot_encode(sequence)
    if policy == "n-channel":
        return x
    if policy not in AMBIGUITY_POLICIES:
        raise ValueError(f"unknown ambiguity policy {policy!r}; expected one of {AMBIGUITY_POLICIES}")
    ambiguous = x[4] > 0
    if not bool(ambiguous.any()):
        return x
    x[4][ambiguous] = 0.0
    if policy == "uniform":
        x[:4, ambiguous] = 0.25
        return x
    # impute: deterministic per window -- seed from the unambiguous content so
    # the same window always yields the same substitution, across runs and hosts
    idx = torch.nonzero(ambiguous, as_tuple=False).flatten()
    seed = (int(x[:4].sum().item()) * 1000003 + int(idx.sum().item())) & 0x7FFFFFFF
    gen = torch.Generator().manual_seed(seed)
    picks = torch.randint(0, 4, (idx.numel(),), generator=gen)
    x[picks, idx] = 1.0
    return x


def collate_batch(batch):
    xs, ys = zip(*batch)
    x = torch.stack(xs, dim=0)  # (B, 5, L)
    y = torch.stack(ys, dim=0)  # (B,)
    return x, y


def stable_sigmoid(logits: np.ndarray) -> np.ndarray:
    # Use a stable sigmoid to avoid overflow warnings on large logits.
    probs = np.empty_like(logits, dtype=np.float32)
    pos = logits >= 0
    probs[pos] = 1.0 / (1.0 + np.exp(-logits[pos]))
    exp_logits = np.exp(logits[~pos])
    probs[~pos] = exp_logits / (1.0 + exp_logits)
    return probs
