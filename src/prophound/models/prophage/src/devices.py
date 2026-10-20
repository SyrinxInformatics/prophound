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
"""Compute device selection."""

from __future__ import annotations

import torch


def mps_available() -> bool:
    return bool(
        hasattr(torch.backends, "mps")
        and torch.backends.mps.is_available()
        and torch.backends.mps.is_built()
    )


def xpu_available() -> bool:
    return bool(hasattr(torch, "xpu") and torch.xpu.is_available())


def rocm_available() -> bool:
    return bool(torch.cuda.is_available() and getattr(torch.version, "hip", None))


def resolve_device(requested_device: str) -> torch.device:
    requested = requested_device.lower()

    if requested == "auto":
        if torch.cuda.is_available():
            return torch.device("cuda")
        if mps_available():
            return torch.device("mps")
        if xpu_available():
            return torch.device("xpu")
        return torch.device("cpu")

    if requested == "cpu":
        return torch.device("cpu")

    if requested == "cuda":
        if not torch.cuda.is_available():
            raise RuntimeError("Requested device 'cuda' is not available.")
        return torch.device("cuda")

    if requested == "rocm":
        if not rocm_available():
            raise RuntimeError(
                "Requested device 'rocm' is not available through the installed PyTorch backend."
            )
        return torch.device("cuda")

    if requested == "mps":
        if not mps_available():
            raise RuntimeError("Requested device 'mps' is not available.")
        return torch.device("mps")

    if requested in {"xpu", "intel"}:
        if not xpu_available():
            raise RuntimeError("Requested device 'xpu' is not available.")
        return torch.device("xpu")

    raise ValueError(f"Unsupported device request: {requested_device}")


def device_supports_amp(device: torch.device) -> bool:
    return device.type == "cuda"
