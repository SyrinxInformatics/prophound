from __future__ import annotations

from typing import List, Literal, Optional

import torch
import torch.nn as nn


CombineMode = Literal["max", "mean"]
PoolingMode = Literal["global_max", "adaptive_max", "adaptive_mean"]
DEFAULT_CHANNELS = [64, 128, 256]
DEFAULT_KERNEL_SIZES = [15, 11, 7]


class CNNEncoder1D(nn.Module):
    """Encode one-hot DNA windows into a pooled embedding."""

    def __init__(
        self,
        in_channels: int = 5,
        channels: Optional[List[int]] = None,
        kernel_sizes: Optional[List[int]] = None,
        dropout: float = 0.2,
        pooling: PoolingMode = "global_max",
        pool_bins: int = 1,
    ):
        super().__init__()
        channels = channels or DEFAULT_CHANNELS.copy()
        kernel_sizes = kernel_sizes or DEFAULT_KERNEL_SIZES.copy()
        if len(channels) != len(kernel_sizes):
            raise ValueError("channels and kernel_sizes must have same length")
        if pool_bins < 1:
            raise ValueError("pool_bins must be at least 1")

        layers: List[nn.Module] = []
        input_channels = in_channels
        for layer_index, (output_channels, kernel_size) in enumerate(
            zip(channels, kernel_sizes)
        ):
            if kernel_size % 2 == 0:
                raise ValueError("Use odd kernel sizes for symmetric 'same' padding.")
            padding = kernel_size // 2

            layers.append(
                nn.Conv1d(
                    input_channels,
                    output_channels,
                    kernel_size=kernel_size,
                    padding=padding,
                )
            )
            layers.append(nn.ReLU(inplace=True))

            if layer_index < 2:
                layers.append(nn.MaxPool1d(kernel_size=2))

            if dropout and dropout > 0:
                layers.append(nn.Dropout(p=dropout))

            input_channels = output_channels

        self.net = nn.Sequential(*layers)
        self.pooling = pooling
        self.pool_bins = pool_bins
        self.output_channels = channels[-1]
        self.emb_dim = channels[-1] * pool_bins

        if pooling == "global_max":
            if pool_bins != 1:
                raise ValueError("global_max pooling requires pool_bins=1")
            self.pool = None
        elif pooling == "adaptive_max":
            self.pool = nn.AdaptiveMaxPool1d(pool_bins)
        elif pooling == "adaptive_mean":
            self.pool = nn.AdaptiveAvgPool1d(pool_bins)
        else:
            raise ValueError(f"Unknown pooling mode: {pooling}")

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        encoded = self.net(x)
        if self.pool is None:
            return torch.amax(encoded, dim=-1)
        pooled = self.pool(encoded)
        return pooled.reshape(pooled.size(0), -1)


class RCCNN(nn.Module):
    """Score DNA windows with a reverse-complement-aware CNN."""

    def __init__(
        self,
        in_channels: int = 5,
        channels: Optional[List[int]] = None,
        kernel_sizes: Optional[List[int]] = None,
        dropout: float = 0.2,
        fc_dim: int = 256,
        combine: CombineMode = "max",
        pooling: PoolingMode = "global_max",
        pool_bins: int = 1,
    ):
        super().__init__()
        self.combine = combine
        self.encoder = CNNEncoder1D(
            in_channels=in_channels,
            channels=channels,
            kernel_sizes=kernel_sizes,
            dropout=dropout,
            pooling=pooling,
            pool_bins=pool_bins,
        )

        self.head = nn.Sequential(
            nn.Linear(self.encoder.emb_dim, fc_dim),
            nn.ReLU(inplace=True),
            nn.Dropout(p=dropout) if dropout and dropout > 0 else nn.Identity(),
            nn.Linear(fc_dim, 1),
        )

    @staticmethod
    def reverse_complement_onehot(x: torch.Tensor) -> torch.Tensor:
        """Reverse-complement a batch of one-hot encoded DNA windows."""
        reversed_x = torch.flip(x, dims=[-1])
        channel_order = torch.tensor([3, 2, 1, 0, 4], device=x.device)
        return reversed_x.index_select(dim=1, index=channel_order)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        forward_embedding = self.encoder(x)
        reverse_embedding = self.encoder(self.reverse_complement_onehot(x))

        if self.combine == "max":
            combined_embedding = torch.maximum(forward_embedding, reverse_embedding)
        elif self.combine == "mean":
            combined_embedding = 0.5 * (forward_embedding + reverse_embedding)
        else:
            raise ValueError(f"Unknown combine mode: {self.combine}")

        return self.head(combined_embedding).squeeze(-1)
