"""Fourier neural operator (FNO, Li et al., ICLR 2021) for 2D boiling fields.

Each Fourier layer applies a learned linear map to the lowest ``modes``
Fourier coefficients of the hidden state, plus a pointwise linear path, then a
GELU nonlinearity. The domain is zero-padded so the spectral convolution does
not impose periodicity on the wall-bounded boiling domain.
"""

from __future__ import annotations

import torch
from torch import nn

__all__ = ["SpectralConv2d", "FNO2d"]


class SpectralConv2d(nn.Module):
    def __init__(self, in_ch: int, out_ch: int, modes_y: int, modes_x: int) -> None:
        super().__init__()
        self.my, self.mx = modes_y, modes_x
        scale = 1.0 / (in_ch * out_ch)
        self.w1 = nn.Parameter(scale * torch.randn(in_ch, out_ch, modes_y, modes_x, dtype=torch.cfloat))
        self.w2 = nn.Parameter(scale * torch.randn(in_ch, out_ch, modes_y, modes_x, dtype=torch.cfloat))

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        b, _, ny, nx = x.shape
        xf = torch.fft.rfft2(x)
        out = torch.zeros(b, self.w1.shape[1], ny, nx // 2 + 1, dtype=torch.cfloat, device=x.device)
        my, mx = self.my, self.mx
        out[:, :, :my, :mx] = torch.einsum("bixy,ioxy->boxy", xf[:, :, :my, :mx], self.w1)
        out[:, :, -my:, :mx] = torch.einsum("bixy,ioxy->boxy", xf[:, :, -my:, :mx], self.w2)
        return torch.fft.irfft2(out, s=(ny, nx))


class FNO2d(nn.Module):
    """FNO mapping ``in_channels`` fields to ``out_channels`` fields on the same grid.

    Two normalised coordinate channels are appended internally.
    """

    def __init__(self, in_channels: int, out_channels: int, width: int = 32, modes: int = 12,
                 layers: int = 4, padding: int = 8) -> None:
        super().__init__()
        self.padding = padding
        self.lift = nn.Conv2d(in_channels + 2, width, 1)
        self.spectral = nn.ModuleList(SpectralConv2d(width, width, modes, modes) for _ in range(layers))
        self.pointwise = nn.ModuleList(nn.Conv2d(width, width, 1) for _ in range(layers))
        self.project = nn.Sequential(nn.Conv2d(width, 2 * width, 1), nn.GELU(), nn.Conv2d(2 * width, out_channels, 1))

    @staticmethod
    def _grid(x: torch.Tensor) -> torch.Tensor:
        b, _, ny, nx = x.shape
        gy = torch.linspace(0, 1, ny, device=x.device).view(1, 1, ny, 1).expand(b, 1, ny, nx)
        gx = torch.linspace(0, 1, nx, device=x.device).view(1, 1, 1, nx).expand(b, 1, ny, nx)
        return torch.cat([gy, gx], dim=1)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x = self.lift(torch.cat([x, self._grid(x)], dim=1))
        p = self.padding
        if p:
            x = nn.functional.pad(x, (0, p, 0, p))
        n = len(self.spectral)
        for i, (s, w) in enumerate(zip(self.spectral, self.pointwise)):
            x = s(x) + w(x)
            if i < n - 1:
                x = nn.functional.gelu(x)
        if p:
            x = x[..., :-p, :-p]
        return self.project(x)
