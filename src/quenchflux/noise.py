"""Measurement-noise estimation from IR sequences."""

from __future__ import annotations

import numpy as np

__all__ = ["estimate_noise_sigma"]


def estimate_noise_sigma(Y: np.ndarray, mask: np.ndarray | None = None) -> float:
    """Robust estimate of white-noise standard deviation in an IR sequence.

    Uses the second temporal difference ``Y[n+1] - 2 Y[n] + Y[n-1]`` whose
    variance is ``6 sigma^2`` for white noise, and the median absolute
    deviation for robustness against the (smooth) signal and outliers.
    Rapid physical transients inflate the estimate slightly; pass the known
    sigma to the inverse solver when it is available.
    """
    Y = np.asarray(Y, dtype=float)
    d2 = Y[2:] - 2.0 * Y[1:-1] + Y[:-2]
    if mask is not None:
        m = np.asarray(mask, bool)
        valid = m[2:] & m[1:-1] & m[:-2]
        d2 = d2[valid]
    d2 = d2.ravel()
    if d2.size == 0:
        raise ValueError("not enough valid frames to estimate the noise level")
    mad = np.median(np.abs(d2 - np.median(d2)))
    return float(1.4826 * mad / np.sqrt(6.0))
