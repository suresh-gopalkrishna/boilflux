"""Adjoint conjugate-gradient method (CGM) for the 3D nonlinear quench IHCP.

Minimises the Tikhonov functional (cf. Eq. (25) of the ATE 2026 paper)::

    J(p) = 1/2 sum_n sum_i m_ni (H T^n(p) - Y^n)_i^2  +  gamma/2 sum_n sum_j p_nj^2

where ``p = q / q_scale`` is the boundary heat flux in MW/m^2 (``q_scale =
1e6``), ``H`` extracts the IR-face temperatures and ``m`` is the observation
mask (camera field of view of the moving plate, bad pixels).

Each iteration needs one forward, one adjoint and one sensitivity solve:

1. direct problem -> residual ``r^n = m (H T^n - Y^n)``
2. adjoint problem (backward in time) -> gradient ``g = dJ/dp``
3. Polak-Ribiere conjugate direction ``d``
4. sensitivity problem ``theta = dT/dp [d]`` -> optimal step
   ``alpha = -(<r, H theta> + gamma <p, d>) / (||H theta||^2 + gamma ||d||^2)``

Properties are frozen on the current trajectory for steps 2 and 4 (quasi-
linearisation); the nonlinearity is resolved by re-running the direct
problem at every iteration, with a backtracking safeguard.

Stopping uses Morozov's discrepancy principle,
``sum(r^2) <= tau^2 * N_obs * sigma^2``, which acts as iterative
regularisation; ``gamma`` (zeroth-order Tikhonov weight) can be zero or a
fixed positive value.
"""

from __future__ import annotations

import time
from collections.abc import Callable
from dataclasses import dataclass, field

import numpy as np

from ..dataset import BoilDataset
from ..forward import ForwardModel, SurfaceLoss, Trajectory
from ..materials import Material, get_material
from ..noise import estimate_noise_sigma

__all__ = ["AdjointCGM", "CGMOptions", "InverseResult"]

Q_SCALE = 1.0e6  # W/m^2 per unit of the optimisation variable


@dataclass
class CGMOptions:
    max_iter: int = 200
    gamma: float = 0.0
    tau: float = 1.02
    stop_on_discrepancy: bool = True
    rel_tol: float = 1e-7
    beta: str = "polak-ribiere"
    verbose: bool = True
    log_every: int = 5


@dataclass
class InverseResult:
    q: np.ndarray                 # (nt, ny, nx) recovered flux [W/m^2]
    T_top: np.ndarray             # (nt, ny, nx) model IR-face temperatures
    T_surface: np.ndarray         # (nt, ny, nx) estimated quenched-face temperature
    trajectory: Trajectory
    history: dict = field(default_factory=dict)
    sigma: float = float("nan")
    converged: bool = False
    message: str = ""
    runtime_s: float = 0.0
    n_reliable: int = 0           # frames before the terminal (end-time) effect

    @property
    def reliable(self) -> slice:
        """Time slice unaffected by the adjoint terminal condition (see docs/theory)."""
        return slice(0, self.n_reliable)

    def save(self, path) -> None:
        np.savez_compressed(
            path, q=self.q, T_top=self.T_top, T_surface=self.T_surface,
            misfit=np.asarray(self.history.get("misfit", [])),
            sigma=self.sigma, converged=self.converged, n_reliable=self.n_reliable,
        )


class AdjointCGM:
    """Adjoint-CGM inverse solver bound to a dataset.

    Parameters
    ----------
    data:
        :class:`~boilflux.dataset.BoilDataset` with IR-face temperatures.
    material:
        Overrides ``data.material`` (useful for user-defined property tables).
    top_loss:
        Heat loss model on the IR face (default: insulated).
    options:
        :class:`CGMOptions`.
    """

    def __init__(
        self,
        data: BoilDataset,
        material: Material | None = None,
        top_loss: SurfaceLoss | None = None,
        options: CGMOptions | None = None,
        solver: str = "auto",
    ) -> None:
        self.data = data
        self.material = material or get_material(data.material)
        self.opt = options or CGMOptions()
        self.model = ForwardModel(data.grid, self.material, data.dt, data.nt, data.T0,
                                  top_loss=top_loss, solver=solver)
        self.mask = data.mask.astype(float)
        self.n_obs = int(self.mask.sum())
        self.sigma = data.sigma if data.sigma is not None else estimate_noise_sigma(data.Y, data.mask)

    # ------------------------------------------------------------------
    def _forward(self, p):
        traj = self.model.run(q=p * Q_SCALE)
        Ttop = traj.top_face()[1:]
        r = self.mask * (Ttop - self.data.Y)
        return traj, Ttop, r

    def _gradient(self, traj, r, p, gamma):
        return Q_SCALE * self.model.adjoint(traj, r) + gamma * p

    def end_effect_frames(self) -> int:
        """Frames at the end of the record biased by the terminal condition.

        The adjoint vanishes at the final time, so the gradient (and the
        recovered flux) is driven to zero over roughly one through-thickness
        diffusion time ``t_d = 0.5 Lz^2 / alpha(T0)``. Results inside this
        window should be discarded; record a correspondingly longer sequence.
        """
        a = float(self.material.diffusivity(np.asarray(self.data.T0)))
        t_d = 0.5 * self.data.grid.Lz**2 / a
        return int(min(self.data.nt - 1, np.ceil(t_d / self.data.dt)))

    def surface_temperature(self, traj: Trajectory, q: np.ndarray) -> np.ndarray:
        """Quenched-face temperature from the first cell layer and the flux."""
        Tc = traj.bottom_cell()[1:]
        return Tc - q * self.data.grid.dz / (2.0 * self.material.k(Tc))

    def solve(
        self,
        p0: np.ndarray | None = None,
        callback: Callable[[int, dict], None] | None = None,
    ) -> InverseResult:
        """Run the CGM. ``p0`` is an optional initial flux in W/m^2 (warm start)."""
        o = self.opt
        g = self.data.grid
        nt = self.data.nt
        t_start = time.perf_counter()
        p = np.zeros((nt, *g.face_shape)) if p0 is None else np.asarray(p0, float) / Q_SCALE
        gamma = float(o.gamma)
        target = o.tau**2 * self.n_obs * self.sigma**2
        hist = {"misfit": [], "rms": [], "gamma": [], "alpha": []}
        traj, Ttop, r = self._forward(p)
        S = float(np.sum(r**2))
        grad_prev = None
        d = None
        converged, message = False, "max_iter reached"
        for k in range(o.max_iter):
            hist["misfit"].append(S)
            hist["rms"].append(np.sqrt(S / self.n_obs))
            hist["gamma"].append(gamma)
            if o.verbose and (k % o.log_every == 0):
                print(f"[CGM] it {k:4d}  RMS misfit {np.sqrt(S / self.n_obs):8.4f} K"
                      f"  (target {np.sqrt(target / self.n_obs):.4f} K)")
            if callback is not None:
                callback(k, {"p": p, "misfit": S})
            if o.stop_on_discrepancy and S <= target:
                converged, message = True, f"discrepancy principle satisfied at iteration {k}"
                break
            grad = self._gradient(traj, r, p, gamma)
            if d is None:
                d = -grad
            else:
                gg = float(np.sum(grad_prev**2))
                if o.beta == "fletcher-reeves":
                    beta = float(np.sum(grad**2)) / gg
                else:
                    beta = max(0.0, float(np.sum(grad * (grad - grad_prev))) / gg)
                d = -grad + beta * d
            if float(np.sum(grad * d)) >= 0:  # not a descent direction: restart
                d = -grad
            th = self.mask * self.model.sensitivity(traj, d * Q_SCALE)
            num = float(np.sum(r * th)) + gamma * float(np.sum(p * d))
            den = float(np.sum(th**2)) + gamma * float(np.sum(d * d))
            if den <= 0:
                message = "zero curvature along search direction"
                break
            alpha = -num / den
            J_old = 0.5 * S + 0.5 * gamma * float(np.sum(p**2))
            # backtracking safeguard for the nonlinearity
            for _ in range(8):
                p_new = p + alpha * d
                traj_new, Ttop_new, r_new = self._forward(p_new)
                S_new = float(np.sum(r_new**2))
                J_new = 0.5 * S_new + 0.5 * gamma * float(np.sum(p_new**2))
                if J_new <= J_old:
                    break
                alpha *= 0.5
            else:
                message = "line search failed"
                break
            hist["alpha"].append(alpha)
            rel = (J_old - J_new) / max(J_old, 1e-300)
            p, traj, Ttop, r, S = p_new, traj_new, Ttop_new, r_new, S_new
            grad_prev = grad
            if rel < o.rel_tol:
                message = f"stagnation (relative decrease {rel:.2e})"
                hist["misfit"].append(S)
                break
        q = p * Q_SCALE
        res = InverseResult(
            q=q, T_top=Ttop, T_surface=self.surface_temperature(traj, q), trajectory=traj,
            history=hist, sigma=self.sigma, converged=converged, message=message,
            runtime_s=time.perf_counter() - t_start,
            n_reliable=nt - self.end_effect_frames(),
        )
        if o.verbose:
            print(f"[CGM] done: {message}; RMS misfit {np.sqrt(S / self.n_obs):.4f} K;"
                  f" {res.runtime_s:.1f} s")
        return res
