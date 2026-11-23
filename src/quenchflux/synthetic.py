"""Synthetic single-jet quench scenarios with known ground-truth heat flux.

A circular water jet, fixed in the laboratory, impinges on a hot plate that
may translate along its length. The local boundary heat flux follows a
zone-dependent Specht boiling curve evaluated at the local surface temperature,
so film, transition and nucleate boiling and the rewetting front emerge from
the physics instead of being prescribed.

To avoid the "inverse crime", data are generated on a grid (and time step)
refined by ``refine`` and then block-averaged onto the inversion grid.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field

import numpy as np

from .boiling import FFZ_DEFAULT, IZ_DEFAULT, SpechtBoilingCurve
from .dataset import QuenchDataset
from .forward import ForwardModel, SurfaceLoss
from .grid import PlateGrid
from .kinematics import PlateMotion
from .materials import get_material

__all__ = ["JetQuenchScenario"]


@dataclass
class JetQuenchScenario:
    """Single round jet on a (moving) plate.

    Lengths in metres, temperatures in degC, times in seconds.
    """

    material: str = "AA6082"
    Lx: float = 0.120
    Ly: float = 0.060
    Lz: float = 0.005
    nx: int = 48
    ny: int = 24
    nz: int = 5
    dt: float = 0.01
    nt: int = 200
    T0: float = 500.0
    plate_velocity: float = 0.01
    x0: float = 0.0
    x_jet_lab: float = 0.045
    y_jet: float | None = None
    r_iz: float = 0.009
    r_ffz: float = 0.011
    r_wet: float = 0.035
    wet_edge: float = 0.004
    iz: SpechtBoilingCurve = field(default_factory=lambda: SpechtBoilingCurve(**asdict(IZ_DEFAULT)))
    ffz: SpechtBoilingCurve = field(default_factory=lambda: SpechtBoilingCurve(**asdict(FFZ_DEFAULT)))
    noise_sigma: float = 2.0
    refine: int = 2
    camera_view_lab: tuple[float, float] | None = None
    top_loss: bool = False
    seed: int = 0

    # ------------------------------------------------------------------
    @property
    def grid(self) -> PlateGrid:
        return PlateGrid(self.Lx, self.Ly, self.Lz, self.nx, self.ny, self.nz)

    @property
    def motion(self) -> PlateMotion:
        return PlateMotion(self.plate_velocity, self.x0)

    def _flux_law(self, grid: PlateGrid, dt: float):
        X, Y = grid.face_mesh()
        X, Y = X.ravel(), Y.ravel()
        yj = self.Ly / 2 if self.y_jet is None else self.y_jet
        motion = self.motion
        span = max(self.r_ffz - self.r_iz, 1e-12)

        def law(T_cell, n):
            t = n * dt
            xs = motion.stagnation_point(self.x_jet_lab, t)
            if grid.ny == 1:
                r = np.abs(X - xs)
            else:
                r = np.hypot(X - xs, Y - yj)
            w = np.clip((self.r_ffz - r) / span, 0.0, 1.0)  # 1 in IZ, 0 in FFZ
            chi = 0.5 * (1.0 - np.tanh((r - self.r_wet) / self.wet_edge))  # liquid coverage
            q1, d1 = self.iz.flux_and_derivative(T_cell)
            q2, d2 = self.ffz.flux_and_derivative(T_cell)
            q = chi * (w * q1 + (1 - w) * q2)
            dq = chi * (w * d1 + (1 - w) * d2)
            return q, dq

        return law

    def generate(self) -> QuenchDataset:
        """Run the fine-grid forward model and return a noisy coarse-grid dataset."""
        mat = get_material(self.material)
        grid = self.grid
        f = int(self.refine)
        fine = grid.refined(f) if f > 1 else grid
        dt_f = self.dt / f
        nt_f = self.nt * f
        loss = SurfaceLoss() if self.top_loss else None
        fm = ForwardModel(fine, mat, dt_f, nt_f, self.T0, top_loss=loss)
        traj = fm.run(flux_fn=self._flux_law(fine, dt_f))
        top = traj.top_face()[f::f]  # (nt, nyf, nxf) at coarse frame times
        q_f = traj.q.reshape(self.nt, f, *fine.face_shape).mean(axis=1)
        if f > 1:
            Y_clean = grid.coarsen_face(top, f)
            q_true = grid.coarsen_face(q_f, f)
        else:
            Y_clean, q_true = top, q_f
        rng = np.random.default_rng(self.seed)
        Y = Y_clean + self.noise_sigma * rng.standard_normal(Y_clean.shape)
        t = np.arange(1, self.nt + 1) * self.dt
        if self.camera_view_lab is not None:
            mask = self.motion.view_mask(grid, t, self.camera_view_lab)
        else:
            mask = np.ones(Y.shape, dtype=bool)
        meta = {k: v for k, v in asdict(self).items() if k not in ("iz", "ffz")}
        meta["iz"] = self.iz.to_dict()
        meta["ffz"] = self.ffz.to_dict()
        meta["generator"] = "quenchflux.synthetic.JetQuenchScenario"
        return QuenchDataset(
            grid=grid, material=self.material, dt=self.dt, T0=self.T0, Y=Y, mask=mask,
            sigma=self.noise_sigma, q_true=q_true, Y_clean=Y_clean, motion=self.motion, meta=meta,
        )

    @classmethod
    def from_dict(cls, d: dict) -> JetQuenchScenario:
        d = dict(d)
        for key in ("iz", "ffz"):
            if key in d and isinstance(d[key], dict):
                d[key] = SpechtBoilingCurve(**d[key])
        if d.get("camera_view_lab") is not None:
            d["camera_view_lab"] = tuple(float(v) for v in d["camera_view_lab"])
        ints = {"nx", "ny", "nz", "nt", "refine", "seed"}
        for key, val in d.items():
            if isinstance(val, str) and key != "material":
                d[key] = int(val) if key in ints else float(val)
        return cls(**d)
