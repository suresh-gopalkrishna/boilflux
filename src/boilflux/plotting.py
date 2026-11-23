"""Matplotlib figures for datasets and inverse results."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import numpy as np

__all__ = [
    "plot_boiling_curves",
    "plot_centerline",
    "plot_convergence",
    "plot_flux_maps",
    "plot_histories",
    "save_all",
]


def _plt():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"figure.dpi": 130, "font.size": 9, "axes.grid": False})
    return plt


def _centre_row(grid):
    return grid.ny // 2


def plot_flux_maps(data, result, times: Sequence[float], path=None):
    """Recovered vs true flux on the quenched face at selected times."""
    plt = _plt()
    g = data.grid
    rows = 2 if data.q_true is not None else 1
    fig, axes = plt.subplots(rows, len(times), figsize=(2.6 * len(times), 1.9 * rows + 0.6),
                             squeeze=False, constrained_layout=True)
    ext = [0, g.Lx * 1e3, 0, g.Ly * 1e3]
    vmax = (data.q_true.max() if data.q_true is not None else result.q.max()) / 1e6
    im = None
    for j, tt in enumerate(times):
        n = int(np.clip(round(tt / data.dt) - 1, 0, data.nt - 1))
        fields = [("recovered", result.q[n])]
        if data.q_true is not None:
            fields.insert(0, ("true", data.q_true[n]))
        for i, (lab, f) in enumerate(fields):
            ax = axes[i, j]
            im = ax.imshow(f / 1e6, origin="lower", extent=ext, vmin=0, vmax=vmax,
                           cmap="inferno", aspect="equal")
            ax.set_title(f"{lab}, t = {(n + 1) * data.dt:.2f} s")
            ax.set_xlabel("x [mm]")
            if j == 0:
                ax.set_ylabel("y [mm]")
    fig.colorbar(im, ax=axes, label="q [MW/m²]", shrink=0.9)
    return _finish(fig, path)


def plot_centerline(data, result, path=None):
    """Space-time map of q along the plate centre line (x vs t)."""
    plt = _plt()
    g = data.grid
    iy = _centre_row(g)
    panels = [("recovered", result.q[:, iy, :])]
    if data.q_true is not None:
        panels.insert(0, ("true", data.q_true[:, iy, :]))
    fig, axes = plt.subplots(1, len(panels), figsize=(4.2 * len(panels), 3.2),
                             constrained_layout=True, squeeze=False)
    vmax = max(p[1].max() for p in panels) / 1e6
    ext = [0, g.Lx * 1e3, data.dt, data.nt * data.dt]
    im = None
    for ax, (lab, f) in zip(axes[0], panels):
        im = ax.imshow(f / 1e6, origin="lower", extent=ext, aspect="auto", cmap="inferno",
                       vmin=0, vmax=vmax)
        if data.motion.velocity != 0 and "x_jet_lab" in data.meta:
            t = data.t
            ax.plot(data.motion.stagnation_point(data.meta["x_jet_lab"], t) * 1e3, t, "c--", lw=0.8,
                    label="stagnation point")
            ax.legend(loc="upper right", fontsize=7)
        n_rel = getattr(result, "n_reliable", 0)
        if lab == "recovered" and 0 < n_rel < data.nt:
            ax.axhspan(n_rel * data.dt, data.nt * data.dt, color="w", alpha=0.35, lw=0)
            ax.text(2, n_rel * data.dt, " end-time effect", color="w", fontsize=7, va="bottom")
        ax.set_xlim(ext[0], ext[1])
        ax.set_title(f"{lab}: centre line")
        ax.set_xlabel("x [mm] (plate frame)")
        ax.set_ylabel("t [s]")
    fig.colorbar(im, ax=axes, label="q [MW/m²]")
    return _finish(fig, path)


def plot_histories(data, result, x_points_mm: Sequence[float], path=None):
    """Heat-flux time histories at points on the centre line."""
    plt = _plt()
    g = data.grid
    iy = _centre_row(g)
    fig, ax = plt.subplots(figsize=(6.2, 3.4), constrained_layout=True)
    colors = plt.cm.viridis(np.linspace(0, 0.9, len(x_points_mm)))
    for c, xm in zip(colors, x_points_mm):
        ix = int(np.clip(xm * 1e-3 / g.dx, 0, g.nx - 1))
        if data.q_true is not None:
            ax.plot(data.t, data.q_true[:, iy, ix] / 1e6, color=c, lw=1.6, alpha=0.55)
        ax.plot(data.t, result.q[:, iy, ix] / 1e6, "--", color=c, lw=1.1, label=f"x = {g.x[ix]*1e3:.1f} mm")
    n_rel = getattr(result, "n_reliable", 0)
    if 0 < n_rel < data.nt:
        ax.axvspan(n_rel * data.dt, data.nt * data.dt, color="0.85", lw=0, label="end-time effect")
    ax.set_xlabel("t [s]")
    ax.set_ylabel("q [MW/m²]")
    ax.set_title("solid: true   dashed: recovered" if data.q_true is not None else "recovered")
    ax.legend(fontsize=7, ncol=2)
    return _finish(fig, path)


def plot_boiling_curves(data, result, x_points_mm: Sequence[float], path=None):
    """Local boiling curves q(T_s) extracted from the reconstruction."""
    plt = _plt()
    g = data.grid
    iy = _centre_row(g)
    fig, ax = plt.subplots(figsize=(5.4, 3.4), constrained_layout=True)
    colors = plt.cm.plasma(np.linspace(0, 0.85, len(x_points_mm)))
    sl = slice(0, getattr(result, "n_reliable", 0) or data.nt)
    for c, xm in zip(colors, x_points_mm):
        ix = int(np.clip(xm * 1e-3 / g.dx, 0, g.nx - 1))
        ax.plot(result.T_surface[sl, iy, ix], result.q[sl, iy, ix] / 1e6, color=c, lw=1.1,
                label=f"x = {g.x[ix]*1e3:.1f} mm")
    ax.set_xlabel("quenched-face temperature T_s [°C]")
    ax.set_ylabel("q [MW/m²]")
    ax.set_title("Reconstructed local boiling curves")
    ax.legend(fontsize=7)
    return _finish(fig, path)


def plot_convergence(result, path=None):
    plt = _plt()
    fig, ax = plt.subplots(figsize=(4.6, 3.0), constrained_layout=True)
    ax.semilogy(result.history["rms"], "k.-", ms=3)
    ax.axhline(result.sigma, color="tab:red", ls="--", lw=0.9, label=f"noise σ = {result.sigma:.2f} K")
    ax.set_xlabel("CGM iteration")
    ax.set_ylabel("RMS temperature misfit [K]")
    ax.legend()
    return _finish(fig, path)


def save_all(data, result, outdir, times: Sequence[float] | None = None,
             x_points_mm: Sequence[float] | None = None) -> list[Path]:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    n_rel = getattr(result, "n_reliable", 0) or data.nt
    tmax = n_rel * data.dt
    times = times or [0.15 * tmax, 0.4 * tmax, 0.7 * tmax, tmax]
    if x_points_mm is None:
        xs0 = data.meta.get("x_jet_lab", data.grid.Lx / 2) - data.motion.x0
        Lmm = data.grid.Lx * 1e3
        x_points_mm = [xs0 * 1e3 + d for d in (-15, -5, 0, 5, 15) if 0 <= xs0 * 1e3 + d <= Lmm]
        if len(x_points_mm) < 4:
            x_points_mm = [xs0 * 1e3 + d for d in (0, 5, 10, 15, 25) if 0 <= xs0 * 1e3 + d <= Lmm]
    paths = [
        plot_flux_maps(data, result, times, outdir / "flux_maps.png"),
        plot_centerline(data, result, outdir / "centerline_xt.png"),
        plot_histories(data, result, x_points_mm, outdir / "histories.png"),
        plot_boiling_curves(data, result, x_points_mm, outdir / "boiling_curves.png"),
        plot_convergence(result, outdir / "convergence.png"),
    ]
    return paths


def _finish(fig, path):
    if path is None:
        return fig
    fig.savefig(path, dpi=150)
    import matplotlib.pyplot as plt

    plt.close(fig)
    return Path(path)
