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
"""Adaptive max pooling must give the same values on every device.

PyTorch's MPS kernel (torch 2.13) is wrong for lengths that are not
a multiple of the bin count, including every input shape the shipped model uses.
"""

import pytest
import torch
import torch.nn.functional as F

from prophound.models.prophage.src.model_rc_cnn import AdaptiveMaxPool1dMPSSafe, adaptive_max_pool1d_by_bins

# 750 / 1250 / 2500: the 3 kb, 5 kb and 10 kb windows after the encoder's two halvings
LENGTHS = [1, 2, 15, 16, 17, 31, 750, 1250, 2000, 2499, 2500, 2501]
BINS = [1, 3, 7, 16, 32]


@pytest.mark.parametrize("length", LENGTHS)
@pytest.mark.parametrize("bins", BINS)
def test_by_bins_equals_native_kernel_on_cpu(length, bins):
    x = torch.randn(2, 4, length, generator=torch.Generator().manual_seed(length * 100 + bins))
    assert torch.equal(adaptive_max_pool1d_by_bins(x, bins), F.adaptive_max_pool1d(x, bins))


def test_module_is_native_off_mps():
    x = torch.randn(2, 4, 2500, generator=torch.Generator().manual_seed(0))
    assert torch.equal(AdaptiveMaxPool1dMPSSafe(16)(x), F.adaptive_max_pool1d(x, 16))


@pytest.mark.skipif(not torch.backends.mps.is_available(), reason="needs an MPS device")
@pytest.mark.parametrize("length", [750, 1250, 2500])
def test_module_on_mps_matches_cpu(length):
    x = torch.randn(2, 256, length, generator=torch.Generator().manual_seed(length))
    assert torch.equal(AdaptiveMaxPool1dMPSSafe(16)(x.to("mps")).cpu(), F.adaptive_max_pool1d(x, 16))
