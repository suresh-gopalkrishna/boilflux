"""Boiling-regime analysis of reconstructed quench data.

From the local heat flux ``q(t)`` and quenched-face temperature ``T_s(t)``
recovered by the inverse solver, this module identifies, at every surface
point, the characteristic points of the transient boiling curve and the
regime sequence a quenched surface passes through::

    film boiling --(rewetting / Leidenfrost point, minimum heat flux)-->
    transition boiling --(maximum heat flux, departure from nucleate boiling)-->
    nucleate boiling --(T_s <= T_sat)--> single-phase forced convection

Definitions used (see docs/physics):

* **Maximum heat flux (MHF / CHF)** ``q_max`` and its surface temperature
  ``T_dnb`` (departure from nucleate boiling): the global maximum of the
  smoothed ``q(t)``.
* **Rewetting (Leidenfrost) point** ``T_leid``, ``t_wet``: the last instant
  before the MHF at which ``q`` is within ``frac`` of the pre-peak minimum,
  ``q <= q_min + frac (q_max - q_min)``; ``q_min`` is the minimum heat flux
  (MHF point of the film-boiling branch).
* Regimes are assigned by **time** relative to ``t_wet`` and ``t_max``,
  which is robust to non-monotonic surface temperatures near a moving jet.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np

__all__ = ["Regime", "BoilingFeatures", "extract_features", "feature_maps", "wetting_front",
           "plot_boiling_curve_regimes", "plot_feature_maps"]


class Regime:
    CONVECTION = 0
    NUCLEATE = 1
    TRANSITION = 2
    FILM = 3
    INACTIVE = 4
    NAMES = {0: "single-phase convection", 1: "nucleate boiling", 2: "transition boiling", 3: "film boiling",
             4: "no significant boiling (dry / unresolved)"}
    COLORS = {0: "#4cc3ff", 1: "#2a9d8f", 2: "#f4a261", 3: "#b3261e", 4: "#d9d9d9"}


@dataclass
class BoilingFeatures:
    q_max: float          # maximum heat flux [W/m^2]
    T_dnb: float          # surface temperature at q_max [degC]
    t_max: float          # time of q_max [s]
    q_min: float          # minimum heat flux on the film branch before rewetting [W/m^2]
    T_leid: float         # rewetting (Leidenfrost) temperature [degC]
    t_wet: float          # rewetting time [s]
    wetted: bool          # False if the peak lies at the end of the record (not rewetted)
    regime: np.ndarray    # (nt,) regime code per sample

    def to_dict(self) -> dict:
        d = asdict(self)
        d.pop("regime")
        return d


def _smooth(q: np.ndarray, w: int) -> np.ndarray:
    if w <= 1:
        return q
    k = np.ones(w) / w
    pad = np.pad(q, (w // 2, w - 1 - w // 2), mode="edge")
    return np.convolve(pad, k, mode="valid")


def extract_features(q: np.ndarray, Ts: np.ndarray, t: np.ndarray, T_sat: float = 100.0,
                     frac: float = 0.05, smooth: int = 5) -> BoilingFeatures:
    """Characteristic boiling-curve points of one surface location."""
    q = np.asarray(q, float)
    Ts = np.asarray(Ts, float)
    t = np.asarray(t, float)
    qs = _smooth(q, smooth)
    i_max = int(np.argmax(qs))
    q_max = float(qs[i_max])
    pre = qs[: i_max + 1]
    i_min = int(np.argmin(pre))
    q_min = float(pre[i_min])
    thr = q_min + frac * (q_max - q_min)
    below = np.nonzero(pre <= thr)[0]
    i_wet = int(below[-1]) if below.size else 0
    wetted = i_max < len(q) - max(2, smooth)
    regime = np.full(len(q), Regime.FILM, dtype=np.int8)
    regime[i_wet:i_max] = Regime.TRANSITION
    regime[i_max:] = Regime.NUCLEATE
    after = np.arange(len(q)) >= i_max
    regime[after & (Ts <= T_sat)] = Regime.CONVECTION
    if not wetted:
        regime[:] = Regime.FILM
        regime[i_wet:] = Regime.TRANSITION
    return BoilingFeatures(q_max=q_max, T_dnb=float(Ts[i_max]), t_max=float(t[i_max]), q_min=q_min,
                           T_leid=float(Ts[i_wet]), t_wet=float(t[i_wet]), wetted=bool(wetted),
                           regime=regime)


def feature_maps(q: np.ndarray, Ts: np.ndarray, t: np.ndarray, min_activity: float = 0.15, **kw) -> dict:
    """Apply :func:`extract_features` to every face cell of ``(nt, ny, nx)`` arrays.

    Cells whose flux never rises by more than ``min_activity`` times the
    global maximum heat flux (``q_max - q_min``) are marked
    :attr:`Regime.INACTIVE`: no liquid contact, or activity below the
    resolution of the reconstruction.
    """
    nt, ny, nx = q.shape
    q_glob = float(np.max(q))
    keys = ("q_max", "T_dnb", "t_max", "q_min", "T_leid", "t_wet")
    out = {k: np.empty((ny, nx)) for k in keys}
    out["wetted"] = np.zeros((ny, nx), bool)
    out["regime"] = np.empty((nt, ny, nx), np.int8)
    for j in range(ny):
        for i in range(nx):
            f = extract_features(q[:, j, i], Ts[:, j, i], t, **kw)
            for k in keys:
                out[k][j, i] = getattr(f, k)
            active = (f.q_max - f.q_min) >= min_activity * q_glob
            out["wetted"][j, i] = f.wetted and active
            out["regime"][:, j, i] = f.regime if active else Regime.INACTIVE
    return out


def wetting_front(t_wet_line: np.ndarray, wetted_line: np.ndarray, x: np.ndarray, t: np.ndarray) -> dict:
    """Wetting-front positions and velocities along a line.

    Returns, for every time, the upstream (min x) and downstream (max x)
    extent of the rewetted region, and their mean velocities (linear fit over
    the interval in which each front moves).
    """
    xl = np.full(t.size, np.nan)
    xr = np.full(t.size, np.nan)
    for n, tt in enumerate(t):
        wet = wetted_line & (t_wet_line <= tt)
        if wet.any():
            xl[n], xr[n] = x[wet].min(), x[wet].max()
    def _mean_speed(xf):
        ok = np.isfinite(xf)
        if ok.sum() < 3 or np.ptp(xf[ok]) == 0:
            return 0.0
        # fit over the interval in which the front actually moves
        moving = ok & (xf != xf[ok][0]) & (xf != xf[ok][-1])
        sel = moving if moving.sum() >= 3 else ok
        return float(np.polyfit(t[sel], xf[sel], 1)[0])

    return {"t": t, "x_upstream": xl, "x_downstream": xr,
            "v_upstream": _mean_speed(xl), "v_downstream": _mean_speed(xr)}


# ---------------------------------------------------------------- plots
def _plt():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"figure.dpi": 130, "font.size": 9})
    return plt


def plot_boiling_curve_regimes(q, Ts, t, path=None, title="Local transient boiling curve", **kw):
    """Boiling curve coloured by regime with MHF and rewetting points marked."""
    plt = _plt()
    f = extract_features(q, Ts, t, **kw)
    fig, ax = plt.subplots(figsize=(5.6, 3.6), constrained_layout=True)
    for code, name in Regime.NAMES.items():
        m = f.regime == code
        if m.any():
            ax.scatter(Ts[m], np.asarray(q)[m] / 1e6, s=7, color=Regime.COLORS[code], label=name)
    ax.plot(Ts, np.asarray(q) / 1e6, color="0.6", lw=0.5, zorder=0)
    ax.plot(f.T_dnb, f.q_max / 1e6, "k^", ms=7, label=f"MHF {f.q_max/1e6:.2f} MW/m² @ {f.T_dnb:.0f} °C")
    ax.plot(f.T_leid, np.interp(f.t_wet, t, q) / 1e6, "kv", ms=7,
            label=f"rewetting (Leidenfrost) @ {f.T_leid:.0f} °C")
    ax.set_xlabel("surface temperature T_s [°C]")
    ax.set_ylabel("q [MW/m²]")
    ax.set_title(title)
    ax.legend(fontsize=7, loc="upper right")
    return _finish(fig, path)


def plot_feature_maps(maps: dict, grid, t_snapshot_index: int, dt: float, path=None):
    """Maps of q_max, T_dnb, rewetting time and the regime field at one instant."""
    plt = _plt()
    from matplotlib.colors import ListedColormap

    ext = [0, grid.Lx * 1e3, 0, grid.Ly * 1e3]
    fig, ax = plt.subplots(2, 2, figsize=(9.0, 4.6), constrained_layout=True)
    panels = [(np.where(maps["wetted"], maps["q_max"] / 1e6, np.nan), "maximum heat flux q_max [MW/m²]",
               "inferno"),
              (np.where(maps["wetted"], maps["T_dnb"], np.nan), "T at q_max (DNB) [°C]", "viridis"),
              (np.where(maps["wetted"], maps["t_wet"], np.nan), "rewetting time t_wet [s]", "cividis")]
    for a, (f, lab, cm) in zip(ax.flat[:3], panels):
        im = a.imshow(f, origin="lower", extent=ext, cmap=cm, aspect="equal")
        fig.colorbar(im, ax=a, shrink=0.85).set_label(lab, fontsize=8)
        a.set_xlabel("x [mm]")
        a.set_ylabel("y [mm]")
    cmap = ListedColormap([Regime.COLORS[k] for k in range(5)])
    a = ax.flat[3]
    a.imshow(maps["regime"][t_snapshot_index], origin="lower", extent=ext, cmap=cmap, vmin=-0.5, vmax=4.5,
             aspect="equal")
    a.set_title(f"regimes at t = {(t_snapshot_index + 1) * dt:.2f} s", fontsize=9)
    a.set_xlabel("x [mm]")
    from matplotlib.patches import Patch

    a.legend(handles=[Patch(color=Regime.COLORS[k], label=v) for k, v in Regime.NAMES.items()],
             fontsize=6, loc="upper right")
    return _finish(fig, path)


def _finish(fig, path):
    if path is None:
        return fig
    fig.savefig(path, dpi=150)
    import matplotlib.pyplot as plt

    plt.close(fig)
    return path
