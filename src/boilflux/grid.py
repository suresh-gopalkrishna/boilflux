"""Structured cell-centred finite-volume grid for a rectangular plate.

Coordinate convention (plate frame)
-----------------------------------
* ``x`` - plate length, also the direction of plate motion   [0, Lx]
* ``y`` - plate width                                         [0, Ly]
* ``z`` - thickness; ``z = 0`` is the quenched (jet) face, ``z = Lz`` the
  infrared-observed face.

Cell arrays have shape ``(nz, ny, nx)`` and are flattened in C order, so the
flat index is ``(iz * ny + iy) * nx + ix``.  Face arrays (unknown heat flux on
``z = 0`` and observations on ``z = Lz``) have shape ``(ny, nx)``.

Setting ``ny = 1`` gives the two-dimensional (length x thickness) model used in
Gopalkrishna et al., Applied Thermal Engineering 304 (2026) 132481.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

__all__ = ["PlateGrid"]


@dataclass
class PlateGrid:
    """Uniform 3D grid of ``nx x ny x nz`` cells over ``Lx x Ly x Lz`` (metres)."""

    Lx: float
    Ly: float
    Lz: float
    nx: int
    ny: int
    nz: int
    _pairs: tuple = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        if min(self.nx, self.ny, self.nz) < 1:
            raise ValueError("nx, ny, nz must be >= 1")
        if self.nz < 2:
            raise ValueError("nz must be >= 2 (through-thickness resolution is required)")
        self._pairs = self._build_pairs()

    # --- geometry -------------------------------------------------------
    @property
    def dx(self) -> float:
        return self.Lx / self.nx

    @property
    def dy(self) -> float:
        return self.Ly / self.ny

    @property
    def dz(self) -> float:
        return self.Lz / self.nz

    @property
    def n_cells(self) -> int:
        return self.nx * self.ny * self.nz

    @property
    def n_face(self) -> int:
        """Number of cells on one z-face (size of the unknown flux vector per step)."""
        return self.nx * self.ny

    @property
    def shape(self) -> tuple[int, int, int]:
        return (self.nz, self.ny, self.nx)

    @property
    def face_shape(self) -> tuple[int, int]:
        return (self.ny, self.nx)

    @property
    def cell_volume(self) -> float:
        return self.dx * self.dy * self.dz

    @property
    def face_area(self) -> float:
        return self.dx * self.dy

    @property
    def x(self) -> np.ndarray:
        """Cell-centre x coordinates [m]."""
        return (np.arange(self.nx) + 0.5) * self.dx

    @property
    def y(self) -> np.ndarray:
        return (np.arange(self.ny) + 0.5) * self.dy

    @property
    def z(self) -> np.ndarray:
        return (np.arange(self.nz) + 0.5) * self.dz

    def face_mesh(self) -> tuple[np.ndarray, np.ndarray]:
        """Return ``(X, Y)`` arrays of shape ``(ny, nx)`` with face-centre coordinates."""
        return np.meshgrid(self.x, self.y)

    # --- index helpers --------------------------------------------------
    def index(self, iz, iy, ix):
        return (np.asarray(iz) * self.ny + np.asarray(iy)) * self.nx + np.asarray(ix)

    @property
    def bottom_cells(self) -> np.ndarray:
        """Flat indices of the cells adjacent to the quenched face z = 0, ordered as the face."""
        return np.arange(self.n_face)

    @property
    def top_cells(self) -> np.ndarray:
        """Flat indices of the cells adjacent to the IR face z = Lz, ordered as the face."""
        return (self.nz - 1) * self.n_face + np.arange(self.n_face)

    def _build_pairs(self):
        """Neighbour pairs (i, j) with geometric factor area/distance for each direction."""
        idx = np.arange(self.n_cells).reshape(self.shape)
        out = []
        # x-direction
        if self.nx > 1:
            out.append((idx[:, :, :-1].ravel(), idx[:, :, 1:].ravel(), self.dy * self.dz / self.dx))
        if self.ny > 1:
            out.append((idx[:, :-1, :].ravel(), idx[:, 1:, :].ravel(), self.dx * self.dz / self.dy))
        out.append((idx[:-1, :, :].ravel(), idx[1:, :, :].ravel(), self.dx * self.dy / self.dz))
        i = np.concatenate([o[0] for o in out])
        j = np.concatenate([o[1] for o in out])
        g = np.concatenate([np.full(o[0].size, o[2]) for o in out])
        return i, j, g

    @property
    def pairs(self):
        return self._pairs

    def refined(self, factor: int, refine_z: bool = True) -> PlateGrid:
        """A grid refined by an integer factor in x and y (and z if ``refine_z``)."""
        fz = factor if refine_z else 1
        ny = self.ny * factor if self.ny > 1 else 1
        return PlateGrid(self.Lx, self.Ly, self.Lz, self.nx * factor, ny, self.nz * fz)

    def coarsen_face(self, field: np.ndarray, factor: int) -> np.ndarray:
        """Block-average a fine face field ``(..., ny*f, nx*f)`` onto this grid's face."""
        fy = factor if self.ny > 1 else 1
        a = np.asarray(field)
        lead = a.shape[:-2]
        a = a.reshape(*lead, self.ny, fy, self.nx, factor)
        return a.mean(axis=(-1, -3))

    def to_dict(self) -> dict:
        return dict(Lx=self.Lx, Ly=self.Ly, Lz=self.Lz, nx=self.nx, ny=self.ny, nz=self.nz)
