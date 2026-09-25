"""1D CNN spectrum encoder with absolute-frequency features.

Convolutions downsample the spectrum into tokens. Each token carries the mean
frequency (Hz) of its receptive field, embedded with multi-scale Fourier
features, because absolute frequency is physical information in ZULF spectra.
Optional Transformer encoder layers let distant peaks interact.
"""
from __future__ import annotations

import math

import torch
from torch import nn

from .config import EncoderConfig


class FrequencyEmbedding(nn.Module):
    def __init__(self, d_model: int, features: int, min_period_hz: float, max_period_hz: float):
        super().__init__()
        periods = torch.exp(torch.linspace(math.log(min_period_hz), math.log(max_period_hz), features))
        self.register_buffer("angular", 2 * math.pi / periods, persistent=False)
        self.project = nn.Linear(2 * features + 1, d_model)

    def forward(self, frequency_hz: torch.Tensor) -> torch.Tensor:
        """frequency_hz: (..., T) -> (..., T, d_model)."""
        phase = frequency_hz.unsqueeze(-1) * self.angular
        scaled = (frequency_hz / 1000.0).unsqueeze(-1)
        return self.project(torch.cat([torch.sin(phase), torch.cos(phase), scaled], dim=-1))


class ConvBlock(nn.Module):
    def __init__(self, cin: int, cout: int, kernel: int, stride: int, dropout: float):
        super().__init__()
        groups = max(1, min(8, cout // 8))
        self.conv1 = nn.Conv1d(cin, cout, kernel, stride=stride, padding=kernel // 2)
        self.norm1 = nn.GroupNorm(groups, cout)
        self.conv2 = nn.Conv1d(cout, cout, kernel, padding=kernel // 2)
        self.norm2 = nn.GroupNorm(groups, cout)
        self.skip = nn.Conv1d(cin, cout, 1, stride=stride)
        self.act = nn.GELU()
        self.drop = nn.Dropout(dropout)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        y = self.act(self.norm1(self.conv1(x)))
        y = self.drop(self.norm2(self.conv2(y)))
        return self.act(y + self.skip(x))


class SpectrumEncoder(nn.Module):
    def __init__(self, config: EncoderConfig):
        super().__init__()
        self.config = config
        blocks, cin = [], config.in_channels
        self.strides = []
        for channels, kernel, stride in config.stages:
            blocks.append(ConvBlock(cin, channels, kernel, stride, config.dropout))
            self.strides.append(stride)
            cin = channels
        self.blocks = nn.ModuleList(blocks)
        self.to_model = nn.Linear(cin, config.d_model)
        self.frequency = FrequencyEmbedding(config.d_model, config.frequency_features,
                                            config.frequency_min_period_hz, config.frequency_max_period_hz)
        if config.transformer_layers:
            layer = nn.TransformerEncoderLayer(config.d_model, config.heads, config.feedforward,
                                               config.dropout, batch_first=True, norm_first=True)
            self.transformer = nn.TransformerEncoder(layer, config.transformer_layers, enable_nested_tensor=False)
        else:
            self.transformer = None
        self.norm = nn.LayerNorm(config.d_model)

    def token_frequencies(self, frequency_hz: torch.Tensor, length: int) -> torch.Tensor:
        """Mean frequency of each token's input span (frequency_hz: (F,) or (B, F))."""
        f = frequency_hz if frequency_hz.dim() == 2 else frequency_hz.unsqueeze(0)
        pooled = nn.functional.adaptive_avg_pool1d(f.unsqueeze(1).float(), length).squeeze(1)
        return pooled

    def forward(self, features: torch.Tensor, frequency_hz: torch.Tensor) -> torch.Tensor:
        """features: (B, C, F); frequency_hz: (F,) or (B, F). Returns tokens (B, T, d_model)."""
        x = features
        for block in self.blocks:
            x = block(x)
        tokens = self.to_model(x.transpose(1, 2))
        if tokens.shape[1] > self.config.max_tokens:
            tokens = nn.functional.adaptive_avg_pool1d(tokens.transpose(1, 2), self.config.max_tokens).transpose(1, 2)
        freq = self.token_frequencies(frequency_hz.to(tokens.device), tokens.shape[1])
        tokens = tokens + self.frequency(freq).expand(tokens.shape[0], -1, -1)
        if self.transformer is not None:
            tokens = self.transformer(tokens)
        return self.norm(tokens)
