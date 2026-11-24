"""Regenerate the figures of the BoilFlux book (docs/book).

    python docs/scripts/make_book_figures.py

Every figure is computed, either from closed-form results derived in the book
or with BoilFlux itself, so the numbers in the text can be reproduced.
"""
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import boilflux as bf  # noqa: E402

OUT = Path("docs/assets/figures/book")
OUT.mkdir(parents=True, exist_ok=True)

INK, BLUE, BRICK, GRAY, TEAL = "#1b2733", "#2b6cb0", "#b5462b", "#8a94a3", "#2c7f78"
plt.rcParams.update({
    "font.size": 9, "axes.edgecolor": INK, "axes.labelcolor": INK, "xtick.color": INK,
    "ytick.color": INK, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": "#e3e6ea", "grid.linewidth": 0.6, "legend.frameon": False,
    "savefig.dpi": 160, "font.family": "DejaVu Sans",
})

# Representative constant properties of AA6082 (book, Chapter 2)
K, RHO, CP = 180.0, 2700.0, 926.0
ALPHA = K / (RHO * CP)


def save(fig, name):
    fig.savefig(OUT / name, bbox_inches="tight")
    plt.close(fig)
    print("wrote", OUT / name)


# --- Fig. 2.3: attenuation through the plate -------------------------------
def fig_attenuation():
    fig, ax = plt.subplots(1, 2, figsize=(7.4, 2.9), constrained_layout=True)
    f = np.logspace(-2, 2, 400)
    for L, c in ((2e-3, TEAL), (5e-3, BLUE), (10e-3, BRICK)):
        lam = np.sqrt(1j * 2 * np.pi * f / ALPHA)
        ax[0].loglog(f, 1 / np.abs(np.cosh(lam * L)), color=c, lw=1.8, label=f"L = {L*1e3:.0f} mm")
    ax[0].set(xlabel="frequency of the surface flux  f  [Hz]", ylabel=r"$|\hat T(L)/\hat T(0)|$",
              ylim=(1e-6, 1.5), title="(a) temporal content")
    ax[0].legend()
    wl = np.logspace(-3.3, -1, 400)
    for L, c in ((2e-3, TEAL), (5e-3, BLUE), (10e-3, BRICK)):
        ax[1].loglog(wl * 1e3, 1 / np.cosh(2 * np.pi / wl * L), color=c, lw=1.8, label=f"L = {L*1e3:.0f} mm")
    ax[1].set(xlabel=r"spatial wavelength  $\Lambda$  [mm]", ylabel=r"$|\hat T(L)/\hat T(0)|$ (steady)",
              ylim=(1e-6, 1.5), title="(b) spatial content")
    save(fig, "attenuation.png")


# --- Fig. 2.2 / 7.1: slab with constant extracted flux ----------------------
def slab_series(zeta, Fo, n_terms=400):
    n = np.arange(1, n_terms + 1)[:, None]
    s = np.sum(2 / (n * np.pi) ** 2 * np.cos(n * np.pi * zeta) * np.exp(-(n * np.pi) ** 2 * Fo), axis=0)
    return -Fo + zeta - zeta**2 / 2 - 1 / 3 + s


def fig_slab():
    L, q = 5e-3, 1.0e6
    mat = bf.Material.constant(K, RHO, CP)
    g = bf.PlateGrid(L, L, L, nx=1, ny=1, nz=40)
    t_ref = L**2 / ALPHA
    dt = t_ref / 2000
    nt = 2000
    fm = bf.ForwardModel(g, mat, dt, nt, 0.0, solver="direct")
    tr = fm.run(q=np.full((nt, 1, 1), q))
    zeta = g.z / L
    fig, ax = plt.subplots(figsize=(4.6, 3.1), constrained_layout=True)
    scale = q * L / K
    err = 0.0
    for Fo, c in ((0.02, TEAL), (0.1, BLUE), (0.3, INK), (1.0, BRICK)):
        n = int(round(Fo * t_ref / dt))
        th_num = tr.T[n].reshape(g.shape)[:, 0, 0]
        zz = np.linspace(0, 1, 200)
        ax.plot(zz, scale * slab_series(zz, Fo), color=c, lw=1.6, label=f"Fo = {Fo:g}")
        ax.plot(zeta[::3], th_num[::3], "o", ms=3.2, mfc="white", mec=c)
        err = max(err, np.max(np.abs(th_num - scale * slab_series(zeta, Fo))))
    ax.set(xlabel=r"$\zeta = z/L$  (0: quenched face, 1: IR face)", ylabel=r"$T - T_0$  [K]")
    ax.legend(title="lines: Eq. (2.14)   dots: BoilFlux", fontsize=8, title_fontsize=8)
    save(fig, "slab_verification.png")
    print(f"slab: max |FV - exact| = {err:.3f} K  (scale qL/k = {scale:.1f} K)")


# --- 1D inverse problem used in Chapters 8, 9 and 11 ------------------------
def build_1d(nt=80, nz=10, L=5e-3, dt=0.01):
    mat = bf.Material.constant(K, RHO, CP)
    g = bf.PlateGrid(L, L, L, nx=1, ny=1, nz=nz)
    fm = bf.ForwardModel(g, mat, dt, nt, 0.0, solver="direct")
    X = np.empty((nt, nt))
    for j in range(nt):
        qj = np.zeros((nt, 1, 1))
        qj[j] = 1.0
        X[:, j] = fm.run(q=qj).top_face()[1:, 0, 0]
    t = np.arange(1, nt + 1) * dt
    q_true = 4e6 * np.exp(-0.5 * ((t - 0.3) / 0.07) ** 2) + 1e6 * (t > 0.45) * np.exp(-(t - 0.45) / 0.1)
    return X, t, q_true


def cgls(X, y, kmax):
    q = np.zeros(X.shape[1])
    r = y.copy()
    s = X.T @ r
    p = s.copy()
    gam = s @ s
    out = []
    for _ in range(kmax):
        w = X @ p
        a = gam / (w @ w)
        q = q + a * p
        r = r - a * w
        s = X.T @ r
        gn = s @ s
        p = s + gn / gam * p
        gam = gn
        out.append((q.copy(), np.linalg.norm(r)))
    return out


def fig_inverse():
    X, t, q_true = build_1d()
    rng = np.random.default_rng(1)
    sigma = 0.3
    y_clean = X @ q_true
    y = y_clean + sigma * rng.standard_normal(y_clean.size)
    U, s, Vt = np.linalg.svd(X)
    print(f"1D IHCP: s_1 = {s[0]:.3e}, s_n = {s[-1]:.3e}, cond = {s[0]/s[-1]:.2e}")

    # Fig. 8.2: singular values and Picard plot
    fig, ax = plt.subplots(figsize=(4.8, 3.1), constrained_layout=True)
    i = np.arange(1, s.size + 1)
    beta = np.abs(U.T @ y)
    ax.semilogy(i, s, "-", color=INK, lw=1.6, label=r"singular values $s_i$")
    ax.semilogy(i, beta, "o", ms=3, color=BLUE, label=r"$|u_i^{T} y|$")
    ax.semilogy(i, beta / s, ".", ms=4, color=BRICK, label=r"$|u_i^{T} y|/s_i$")
    ax.axhline(sigma, color=GRAY, ls="--", lw=1)
    ax.text(1, sigma * 0.02, r"noise level $\sigma$", color=GRAY, fontsize=8)
    ax.set(xlabel="index i", ylim=(1e-10, 1e10))
    ax.legend(fontsize=8, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.16))
    save(fig, "picard.png")

    # Fig. 9.1: Tikhonov filter factors
    fig, ax = plt.subplots(figsize=(4.6, 2.9), constrained_layout=True)
    for gam, c in ((1e-14, TEAL), (1e-12, BLUE), (1e-10, BRICK)):
        ax.semilogy(i, s**2 / (s**2 + gam), "-", color=c, lw=1.6, label=rf"$\gamma = 10^{{{int(np.log10(gam))}}}$")
    ax.set(xlabel="index i", ylabel=r"filter factor $f_i = s_i^2/(s_i^2+\gamma)$", ylim=(1e-6, 2))
    ax.legend(fontsize=8)
    save(fig, "filter_factors.png")

    # Fig. 9.2: reconstructions
    N = y.size
    delta = 1.02 * np.sqrt(N) * sigma
    q_ls = Vt.T @ (beta * np.sign(U.T @ y) / s)
    gammas = np.logspace(-16, -6, 200)
    res = []
    for gam in gammas:
        f = s**2 / (s**2 + gam)
        qg = Vt.T @ (f * (U.T @ y) / s)
        res.append(np.linalg.norm(X @ qg - y))
    res = np.array(res)
    gam_dp = gammas[np.argmin(np.abs(res - delta))]
    f = s**2 / (s**2 + gam_dp)
    q_tik = Vt.T @ (f * (U.T @ y) / s)
    it = cgls(X, y, 60)
    rn = np.array([r for _, r in it])
    k_dp = int(np.argmax(rn <= delta))
    q_cg = it[k_dp][0]
    errs = np.array([np.linalg.norm(q - q_true) / np.linalg.norm(q_true) for q, _ in it])
    print(f"Tikhonov (discrepancy): gamma = {gam_dp:.2e}, rel err = "
          f"{np.linalg.norm(q_tik-q_true)/np.linalg.norm(q_true):.3f}")
    print(f"CGLS stopped at k = {k_dp+1}, rel err = {errs[k_dp]:.3f}; best k = {np.argmin(errs)+1}"
          f" ({errs.min():.3f}); unregularised rel err = "
          f"{np.linalg.norm(q_ls-q_true)/np.linalg.norm(q_true):.2e}")

    fig, ax = plt.subplots(1, 2, figsize=(7.6, 2.9), constrained_layout=True)
    ax[0].plot(t, q_true / 1e6, color=INK, lw=2.2, label="true")
    ax[0].plot(t, q_tik / 1e6, color=BLUE, lw=1.5, label="Tikhonov, discrepancy")
    ax[0].plot(t, q_cg / 1e6, color=BRICK, lw=1.5, ls="--", label=f"CG stopped (k = {k_dp+1})")
    ax[0].set(xlabel="t [s]", ylabel="q [MW/m²]", ylim=(-1.5, 5.5), title="(a) reconstructions")
    ax[0].legend(fontsize=8)
    k = np.arange(1, len(it) + 1)
    ax[1].semilogy(k, errs, color=BRICK, lw=1.6, label="relative error in q")
    ax[1].semilogy(k, rn / rn[0], color=BLUE, lw=1.6, label="relative residual")
    ax[1].axhline(delta / rn[0], color=GRAY, ls="--", lw=1)
    ax[1].axvline(k_dp + 1, color=GRAY, ls=":", lw=1)
    ax[1].text(k_dp + 2, 0.5, "discrepancy\nstop", color=GRAY, fontsize=8)
    ax[1].set(xlabel="CG iteration k", title="(b) semi-convergence")
    ax[1].legend(fontsize=8, loc="center right")
    save(fig, "regularisation_1d.png")


# --- Fig. 3.2: van der Waals spinodal ---------------------------------------
def fig_vdw():
    v = np.linspace(0.42, 5, 2000)
    fig, ax = plt.subplots(figsize=(4.6, 3.1), constrained_layout=True)
    for Tr, c in ((0.80, TEAL), (27 / 32, BRICK), (0.90, BLUE), (1.0, INK)):
        p = 8 * Tr / (3 * v - 1) - 3 / v**2
        lab = r"$T_r = 27/32$" if abs(Tr - 27 / 32) < 1e-9 else rf"$T_r = {Tr:.2f}$"
        ax.plot(v, p, color=c, lw=1.6 if c != BRICK else 2.0, label=lab)
    vs = np.linspace(0.5, 5, 600)
    ps = (8 * ((3 * vs - 1) ** 2 / (4 * vs**3)) / (3 * vs - 1) - 3 / vs**2)
    ax.plot(vs, ps, color=GRAY, ls="--", lw=1.2, label="spinodal")
    ax.plot([2 / 3], [0], "o", color=BRICK, ms=5)
    ax.annotate(r"liquid spinodal at $p_r = 0$" "\n" r"$v_r = 2/3,\ T_r = 27/32$", (2 / 3, 0), (1.5, -1.6),
                arrowprops=dict(arrowstyle="->", color=INK), fontsize=8)
    ax.axhline(0, color=INK, lw=0.6)
    ax.set(xlim=(0.4, 4), ylim=(-3, 2.2), xlabel=r"reduced volume $v_r$", ylabel=r"reduced pressure $p_r$")
    ax.legend(fontsize=8, loc="upper right")
    save(fig, "vdw_spinodal.png")


# --- Fig. 5.2: conduction-controlled rewetting -------------------------------
def fig_rewetting():
    th = np.logspace(-2, 2, 400)
    fig, ax = plt.subplots(figsize=(4.6, 3.0), constrained_layout=True)
    ax.loglog(th, 1 / np.sqrt(th * (th + 1)), color=INK, lw=2, label=r"Eq. (5.9): $1/\sqrt{\Theta(\Theta+1)}$")
    ax.loglog(th, 1 / th, color=BRICK, ls="--", lw=1.2, label=r"$\Theta \gg 1$: $1/\Theta$")
    ax.loglog(th, 1 / np.sqrt(th), color=BLUE, ls="--", lw=1.2, label=r"$\Theta \ll 1$: $\Theta^{-1/2}$")
    ax.set(xlabel=r"dry-wall superheat ratio $\Theta = (T_w - T_0)/(T_0 - T_\mathrm{sat})$",
           ylabel=r"$u\,/\,(\alpha\sqrt{h/(k\delta)})$")
    ax.legend(fontsize=8)
    save(fig, "rewetting_speed.png")


# --- Fig. 14.2: heat-equation step as a Fourier multiplier ------------------
def fig_multiplier():
    k = np.arange(0, 33)
    fig, ax = plt.subplots(figsize=(4.6, 3.0), constrained_layout=True)
    for nu, c in ((1e-4, TEAL), (1e-3, BLUE), (1e-2, BRICK)):
        ax.semilogy(k, np.exp(-((2 * np.pi * k) ** 2) * nu), "o-", ms=2.5, lw=1.3, color=c,
                    label=rf"$\alpha\Delta t/\ell^2 = 10^{{{int(np.log10(nu))}}}$")
    ax.axvspan(0, 12, color=BLUE, alpha=0.07, lw=0)
    ax.text(1, 2e-6, "modes kept by the FNO\n(k ≤ 12)", fontsize=8, color=BLUE)
    ax.set(xlabel="Fourier mode k", ylabel=r"multiplier $e^{-\alpha (2\pi k/\ell)^2 \Delta t}$", ylim=(1e-8, 2))
    ax.legend(fontsize=8, loc="lower left")
    save(fig, "heat_multiplier.png")


if __name__ == "__main__":
    fig_attenuation()
    fig_slab()
    fig_inverse()
    fig_vdw()
    fig_rewetting()
    fig_multiplier()
