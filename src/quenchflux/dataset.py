"""Container for IR quench data (synthetic or experimental) with NPZ I/O."""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .grid import PlateGrid
from .kinematics import PlateMotion

__all__ = ["QuenchDataset"]

FORMAT_VERSION = "1.0"


@dataclass
class QuenchDataset:
    """IR-face temperature data on the plate-frame grid.

    Attributes
    ----------
    grid:
        Plate grid on which ``Y`` and ``mask`` are defined.
    material:
        Name of a built-in material (or a label for a user material).
    dt:
        Frame interval [s]; frames are at ``t_n = n * dt``, n = 1..nt.
    T0:
        Initial (uniform) plate temperature [degC].
    Y:
        Measured IR-face temperatures, shape ``(nt, ny, nx)`` [degC].
    mask:
        Boolean validity mask, same shape (camera field of view, bad pixels).
    sigma:
        Measurement noise standard deviation [K] if known.
    q_true:
        Ground-truth flux (synthetic data only), shape ``(nt, ny, nx)`` [W/m^2].
    motion:
        Plate kinematics.
    meta:
        Free-form metadata (scenario parameters, experiment notes, ...).
    """

    grid: PlateGrid
    material: str
    dt: float
    T0: float
    Y: np.ndarray
    mask: np.ndarray | None = None
    sigma: float | None = None
    q_true: np.ndarray | None = None
    Y_clean: np.ndarray | None = None
    motion: PlateMotion = field(default_factory=PlateMotion)
    meta: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.Y = np.asarray(self.Y, dtype=float)
        if self.Y.shape[1:] != self.grid.face_shape:
            raise ValueError(f"Y has face shape {self.Y.shape[1:]}, grid expects {self.grid.face_shape}")
        if self.mask is None:
            self.mask = np.ones(self.Y.shape, dtype=bool)
        self.mask = np.asarray(self.mask, dtype=bool)

    @property
    def nt(self) -> int:
        return self.Y.shape[0]

    @property
    def t(self) -> np.ndarray:
        return np.arange(1, self.nt + 1) * self.dt

    def save(self, path: str | Path) -> Path:
        path = Path(path)
        header = {
            "format": "quenchflux-dataset",
            "version": FORMAT_VERSION,
            "grid": self.grid.to_dict(),
            "material": self.material,
            "dt": self.dt,
            "T0": self.T0,
            "sigma": self.sigma,
            "motion": {"velocity": self.motion.velocity, "x0": self.motion.x0},
            "meta": self.meta,
        }
        arrays = {"Y": self.Y, "mask": self.mask}
        if self.q_true is not None:
            arrays["q_true"] = self.q_true
        if self.Y_clean is not None:
            arrays["Y_clean"] = self.Y_clean
        np.savez_compressed(path, header=json.dumps(header), **arrays)
        return path if path.suffix == ".npz" else path.with_suffix(path.suffix + ".npz")

    @classmethod
    def load(cls, path: str | Path) -> QuenchDataset:
        with np.load(path, allow_pickle=False) as f:
            header = json.loads(str(f["header"]))
            arrays = {k: f[k] for k in f.files if k != "header"}
        if header.get("format") != "quenchflux-dataset":
            raise ValueError(f"{path} is not a QuenchFlux dataset")
        return cls(
            grid=PlateGrid(**header["grid"]),
            material=header["material"],
            dt=header["dt"],
            T0=header["T0"],
            Y=arrays["Y"],
            mask=arrays.get("mask"),
            sigma=header.get("sigma"),
            q_true=arrays.get("q_true"),
            Y_clean=arrays.get("Y_clean"),
            motion=PlateMotion(**header.get("motion", {})),
            meta=header.get("meta", {}),
        )
