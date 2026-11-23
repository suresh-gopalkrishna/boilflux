"""Nonlinear transient heat-conduction solver (direct problem).

Discretisation
--------------
Cell-centred finite volumes in space, backward (implicit) Euler in time.
Temperature-dependent properties are evaluated at the previous time level
(semi-implicit linearisation), so each step requires one sparse SPD solve::

    A_n T^n = C_n T^{n-1} - B q^n + c_n ,     A_n = C_n + K(T^{n-1}) + L_n

with ``C_n = diag(rho c V / dt)`` at ``T^{n-1}``, ``K`` the conduction matrix
(harmonic-mean face conductances), ``B`` the face-area operator of the quenched
face, and ``L_n, c_n`` the linearised optional heat loss on the IR face.

Sign convention: a positive heat flux ``q`` extracts heat from the plate
through the quenched face ``z = 0`` (as in the ATE 2026 paper).

Because every step is affine in ``T^n`` with coefficients frozen at the stored
trajectory, the exact discrete adjoint and sensitivity problems of the
linearised map are available (see :mod:`quenchflux.inverse.cgm`).
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from .grid import PlateGrid
from .materials import Material

__all__ = ["ForwardModel", "SurfaceLoss", "Trajectory"]

SIGMA_SB = 5.670374419e-8


@dataclass
class SurfaceLoss:
    """Convective + radiative heat loss on the IR face z = Lz.

    ``q_loss = h (T - T_amb) + eps * sigma * ((T+273.15)^4 - (T_amb+273.15)^4)``
    """

    h: float = 10.0
    emissivity: float = 0.9
    T_amb: float = 25.0

    def flux(self, T: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        Tk = T + 273.15
        Ta = self.T_amb + 273.15
        q = self.h * (T - self.T_amb) + self.emissivity * SIGMA_SB * (Tk**4 - Ta**4)
        dq = self.h + 4.0 * self.emissivity * SIGMA_SB * Tk**3
        return q, dq


@dataclass
class Trajectory:
    """Result of a forward run.

    Attributes
    ----------
    T:  (nt+1, n_cells) temperatures at t_0..t_nt [degC]
    q:  (nt, ny, nx) heat flux applied on the quenched face at t_1..t_nt [W/m^2]
    t:  (nt+1,) time levels [s]
    """

    T: np.ndarray
    q: np.ndarray
    t: np.ndarray
    grid: PlateGrid

    def field(self, n: int) -> np.ndarray:
        return self.T[n].reshape(self.grid.shape)

    def top_face(self) -> np.ndarray:
        """IR-face temperatures (top-cell values), shape (nt+1, ny, nx)."""
        g = self.grid
        return self.T[:, g.top_cells].reshape(-1, *g.face_shape)

    def bottom_cell(self) -> np.ndarray:
        g = self.grid
        return self.T[:, g.bottom_cells].reshape(-1, *g.face_shape)


class ForwardModel:
    """Implicit FV solver for transient nonlinear conduction in a plate.

    Parameters
    ----------
    grid, material:
        Discretisation and properties.
    dt, nt:
        Time-step size [s] and number of steps.
    T0:
        Initial temperature, scalar or array of shape ``grid.shape`` [degC].
    top_loss:
        Optional :class:`SurfaceLoss` on the IR face (default: insulated).
    solver:
        ``"cg"`` (Jacobi-preconditioned conjugate gradients, default),
        ``"direct"`` (sparse LU), or ``"auto"``.
    """

    def __init__(
        self,
        grid: PlateGrid,
        material: Material,
        dt: float,
        nt: int,
        T0: float | np.ndarray,
        top_loss: SurfaceLoss | None = None,
        solver: str = "auto",
        rtol: float = 1e-11,
    ) -> None:
        self.grid = grid
        self.material = material
        self.dt = float(dt)
        self.nt = int(nt)
        T0 = np.asarray(T0, dtype=float)
        if T0.ndim:
            self.T0 = np.broadcast_to(T0, grid.shape).ravel().copy()
        else:
            self.T0 = np.full(grid.n_cells, float(T0))
        self.top_loss = top_loss
        if solver == "auto":
            solver = "direct" if grid.n_cells <= 4000 else "cg"
        if solver not in ("cg", "direct"):
            raise ValueError("solver must be 'cg', 'direct' or 'auto'")
        self.solver = solver
        self.rtol = rtol
        i, j, g = grid.pairs
        self._i, self._j, self._g = i, j, g
        n = grid.n_cells
        # Sparsity pattern: off-diagonals (i,j),(j,i) then the diagonal
        self._rows = np.concatenate([i, j, np.arange(n)])
        self._cols = np.concatenate([j, i, np.arange(n)])
        self._area = grid.face_area

    # --- assembly -------------------------------------------------------
    def coefficients(self, T_prev: np.ndarray):
        """Frozen coefficients for the step starting from ``T_prev``.

        Returns ``(A, cap, c)`` with ``A`` a CSR SPD matrix, ``cap = rho c V / dt``
        and ``c`` the affine loss term (zeros if no loss).
        """
        g = self.grid
        mat = self.material
        k = mat.k(T_prev)
        cap = mat.rho_cp(T_prev) * g.cell_volume / self.dt
        ki, kj = k[self._i], k[self._j]
        G = self._g * 2.0 * ki * kj / (ki + kj)  # harmonic mean conductance
        n = g.n_cells
        diag = cap.copy()
        np.add.at(diag, self._i, G)
        np.add.at(diag, self._j, G)
        c = np.zeros(n)
        if self.top_loss is not None:
            top = g.top_cells
            ql, dql = self.top_loss.flux(T_prev[top])
            diag[top] += dql * self._area
            c[top] = -(ql - dql * T_prev[top]) * self._area
        data = np.concatenate([-G, -G, diag])
        A = sp.csr_matrix((data, (self._rows, self._cols)), shape=(n, n))
        return A, cap, c

    def solve(self, A, b, x0=None, extra_diag=None):
        if extra_diag is not None:
            A = A + sp.diags(extra_diag)
        if self.solver == "direct":
            return spla.spsolve(A.tocsc(), b)
        M = sp.diags(1.0 / A.diagonal())
        x, info = spla.cg(A, b, x0=x0, rtol=self.rtol, atol=0.0, M=M, maxiter=2000)
        if info != 0:  # fall back to a direct solve
            x = spla.spsolve(A.tocsc(), b)
        return x

    # --- time integration ----------------------------------------------
    def run(
        self,
        q: np.ndarray | None = None,
        flux_fn: Callable | None = None,
    ) -> Trajectory:
        """Integrate the direct problem.

        Either give the boundary flux history ``q`` with shape ``(nt, ny, nx)``
        (W/m^2, at t_1..t_nt), or a nonlinear boundary law
        ``flux_fn(T_cell, n) -> (q, dq_dT)`` evaluated on the quenched-face
        cells; the law is linearised implicitly about ``T^{n-1}``.
        """
        g = self.grid
        if (q is None) == (flux_fn is None):
            raise ValueError("give exactly one of q or flux_fn")
        bot = g.bottom_cells
        A_face = self._area
        T = np.empty((self.nt + 1, g.n_cells))
        T[0] = self.T0
        q_out = np.empty((self.nt, g.n_face))
        if q is not None:
            q = np.asarray(q, dtype=float).reshape(self.nt, g.n_face)
        for n in range(1, self.nt + 1):
            Tp = T[n - 1]
            A, cap, c = self.coefficients(Tp)
            b = cap * Tp + c
            if q is not None:
                b[bot] -= q[n - 1] * A_face
                T[n] = self.solve(A, b, x0=Tp)
                q_out[n - 1] = q[n - 1]
            else:
                qs, dq = flux_fn(Tp[bot], n)
                extra = np.zeros(g.n_cells)
                extra[bot] = dq * A_face
                b[bot] -= (qs - dq * Tp[bot]) * A_face
                T[n] = self.solve(A, b, x0=Tp, extra_diag=extra)
                q_out[n - 1] = qs + dq * (T[n][bot] - Tp[bot])
        t = np.arange(self.nt + 1) * self.dt
        return Trajectory(T=T, q=q_out.reshape(self.nt, *g.face_shape), t=t, grid=g)

    # --- linearised problems (used by the inverse solver) ---------------
    def sensitivity(self, traj: Trajectory, dq: np.ndarray) -> np.ndarray:
        """Directional derivative of the IR-face temperatures w.r.t. the flux.

        Solves ``A_n th^n = C_n th^{n-1} - B dq^n`` with coefficients frozen on
        ``traj``. Returns the top-face response, shape (nt, ny, nx).
        """
        g = self.grid
        dq = np.asarray(dq).reshape(self.nt, g.n_face)
        th = np.zeros(g.n_cells)
        out = np.empty((self.nt, g.n_face))
        for n in range(1, self.nt + 1):
            A, cap, _ = self.coefficients(traj.T[n - 1])
            b = cap * th
            b[g.bottom_cells] -= dq[n - 1] * self._area
            th = self.solve(A, b, x0=th)
            out[n - 1] = th[g.top_cells]
        return out.reshape(self.nt, *g.face_shape)

    def adjoint(self, traj: Trajectory, weighted_residual: np.ndarray) -> np.ndarray:
        """Gradient of ``0.5 * sum(w * r^2)`` with respect to the flux history.

        ``weighted_residual`` = ``w * (H T^n - Y^n)`` with shape (nt, ny, nx).
        Solves the discrete adjoint backward in time and returns ``dJ/dq``
        with shape (nt, ny, nx) [K^2 per W/m^2 in the misfit units].
        """
        g = self.grid
        wr = np.asarray(weighted_residual).reshape(self.nt, g.n_face)
        mu = np.zeros(g.n_cells)
        cap_next = None
        grad = np.empty((self.nt, g.n_face))
        for n in range(self.nt, 0, -1):
            A, cap, _ = self.coefficients(traj.T[n - 1])
            b = np.zeros(g.n_cells) if cap_next is None else cap_next * mu
            b[g.top_cells] += wr[n - 1]
            mu = self.solve(A, b, x0=mu)
            grad[n - 1] = -mu[g.bottom_cells] * self._area
            cap_next = cap
        return grad.reshape(self.nt, *g.face_shape)
