"""Moving-plate kinematics and camera (lab-frame) to plate-frame mapping.

In a moving-plate quench test the jet and the IR camera are fixed in the
laboratory while the plate travels with velocity ``v`` along its length ``x``.
QuenchFlux solves conduction in the *plate frame* (Lagrangian description), so
no advection term or artificial inflow boundary is needed. The price is that
the jet stagnation point moves across the plate and the camera field of view
becomes a time-dependent observation mask; both are handled here.

Lab position of a plate point:  ``X_lab = X0 + v * t + x_plate``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .grid import PlateGrid

__all__ = ["PlateMotion", "lab_frames_to_plate"]


@dataclass
class PlateMotion:
    """Rigid translation of the plate along x.

    Parameters
    ----------
    velocity:
        Plate speed along +x in the lab frame [m/s]; 0 for a stationary plate.
    x0:
        Lab-frame position of the plate edge ``x_plate = 0`` at ``t = 0`` [m].
    """

    velocity: float = 0.0
    x0: float = 0.0

    def offset(self, t):
        return self.x0 + self.velocity * np.asarray(t, dtype=float)

    def lab_to_plate(self, X_lab, t):
        return np.asarray(X_lab, dtype=float) - self.offset(t)

    def plate_to_lab(self, x_plate, t):
        return np.asarray(x_plate, dtype=float) + self.offset(t)

    def stagnation_point(self, x_jet_lab: float, t):
        """Plate-frame x position of a jet fixed at ``x_jet_lab`` in the lab."""
        return self.lab_to_plate(x_jet_lab, t)

    def view_mask(self, grid: PlateGrid, t, x_view_lab: tuple[float, float]) -> np.ndarray:
        """Boolean mask ``(len(t), ny, nx)`` of plate face cells inside the camera view."""
        t = np.atleast_1d(t)
        X = self.plate_to_lab(grid.x[None, :], t[:, None])  # (nt, nx)
        m = (X >= x_view_lab[0]) & (X <= x_view_lab[1])
        return np.broadcast_to(m[:, None, :], (t.size, grid.ny, grid.nx)).copy()


def lab_frames_to_plate(
    frames: np.ndarray,
    X_px: np.ndarray,
    Y_px: np.ndarray,
    t: np.ndarray,
    grid: PlateGrid,
    motion: PlateMotion,
    y_offset: float = 0.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Resample lab-frame IR frames onto the plate-frame face grid.

    Parameters
    ----------
    frames:
        IR temperatures, shape ``(nt, n_py, n_px)`` [degC] (already converted
        from radiance with the calibrated emissivity).
    X_px, Y_px:
        1D lab-frame pixel-centre coordinates along x (``n_px``) and y
        (``n_py``) [m], both increasing.
    t:
        Frame times [s], shape ``(nt,)``.
    grid, motion:
        Plate grid and kinematics.
    y_offset:
        Lab y coordinate of the plate edge ``y_plate = 0``.

    Returns
    -------
    data, mask:
        ``(nt, ny, nx)`` temperatures and boolean validity mask on the plate grid.
    """
    frames = np.asarray(frames, dtype=float)
    nt = frames.shape[0]
    xs, ys = grid.x, grid.y + y_offset
    data = np.zeros((nt, grid.ny, grid.nx))
    mask = np.zeros((nt, grid.ny, grid.nx), dtype=bool)
    # interpolate along y first (time-independent), then along x per frame
    if grid.ny == 1 and len(Y_px) == 1:
        fy = frames[:, :1, :]
        yvalid = np.ones(1, dtype=bool)
    else:
        yvalid = (ys >= Y_px[0]) & (ys <= Y_px[-1])
        jy = np.clip(np.searchsorted(Y_px, ys) - 1, 0, len(Y_px) - 2)
        wy = np.clip((ys - Y_px[jy]) / (Y_px[jy + 1] - Y_px[jy]), 0.0, 1.0)
        fy = frames[:, jy, :] * (1 - wy)[None, :, None] + frames[:, jy + 1, :] * wy[None, :, None]
    for n in range(nt):
        Xl = motion.plate_to_lab(xs, t[n])
        xvalid = (Xl >= X_px[0]) & (Xl <= X_px[-1])
        jx = np.clip(np.searchsorted(X_px, Xl) - 1, 0, len(X_px) - 2)
        wx = np.clip((Xl - X_px[jx]) / (X_px[jx + 1] - X_px[jx]), 0.0, 1.0)
        data[n] = fy[n][:, jx] * (1 - wx) + fy[n][:, jx + 1] * wx
        mask[n] = yvalid[:, None] & xvalid[None, :]
    data[~mask] = 0.0
    return data, mask
