"""Regenerate the figures and reference numbers of the physics guide.

python docs/scripts/make_physics_figures.py
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import boilflux as qf

OUT = "docs/assets/figures/physics/"
g = 9.81
# saturated water at 1 atm
rl, rv, hfg, sig, Tsat_K = 958.4, 0.597, 2.257e6, 0.0589, 373.15
Tc_K = 647.1
L_RT = np.sqrt(sig / (g * (rl - rv)))
q_chf_zuber = 0.131 * rv * hfg * (sig * g * (rl - rv) / rv**2) ** 0.25
q_chf_ld = 0.149 * rv * hfg * (sig * g * (rl - rv) / rv**2) ** 0.25
q_min_ber = 0.09 * rv * hfg * (sig * g * (rl - rv) / (rl + rv) ** 2) ** 0.25
T_spiegler = 27 / 32 * Tc_K - 273.15
lam_T = 2 * np.pi * np.sqrt(3) * L_RT
print(f"capillary length            {L_RT*1e3:.2f} mm")
print(f"most dangerous RT wavelength {lam_T*1e3:.1f} mm")
print(f"CHF Zuber (C=0.131)          {q_chf_zuber/1e6:.2f} MW/m2")
print(f"CHF Lienhard-Dhir (0.149)    {q_chf_ld/1e6:.2f} MW/m2")
print(f"MHF Berenson                 {q_min_ber/1e3:.1f} kW/m2")
print(f"Spiegler Leidenfrost estimate {T_spiegler:.0f} C")


def effusivity(k, rho, cp):
    return np.sqrt(k * rho * cp)


e_w = effusivity(0.68, 958, 4216)  # water near 100 C
for name, mat, T in (("AA6082", "aa6082", 500.0), ("Nickel", "nickel", 800.0), ("C45 steel", "c45", 800.0)):
    m = qf.get_material(mat)
    e_s = effusivity(m.k(T), m.rho(T), m.cp(T))
    Ti = (e_s * T + e_w * 20.0) / (e_s + e_w)
    print(f"contact temperature {name:9s} T0={T:.0f} C, water 20 C: e_s={e_s:.0f}, T_i={Ti:.0f} C")

# --- annotated boiling curve -------------------------------------------
bc = qf.IZ_DEFAULT
T = np.linspace(25, 650, 2000)
q = bc(T) / 1e6
fig, ax = plt.subplots(figsize=(7.2, 4.0), constrained_layout=True)
spans = [(20, 100, "#4cc3ff", "single-phase\nconvection"),
         (100, bc.T_dnb, "#2a9d8f", "nucleate\nboiling"),
         (bc.T_dnb, bc.T_leid, "#f4a261", "transition\nboiling"),
         (bc.T_leid, 650, "#b3261e", "film\nboiling")]
for a, b, c, lab in spans:
    ax.axvspan(a, b, color=c, alpha=0.13, lw=0)
    ax.text((a + b) / 2, 7.05, lab, ha="center", va="top", fontsize=8, color=c)
ax.plot(T, q, "k-", lw=2)
ax.plot(bc.T_dnb, bc.q_max / 1e6, "k^", ms=8)
ax.annotate("maximum heat flux (CHF / DNB)\nq_max at T_DNB", (bc.T_dnb, bc.q_max / 1e6), (225, 4.4),
            arrowprops=dict(arrowstyle="->"), fontsize=8)
ax.plot(bc.T_leid, bc.q_leid / 1e6, "kv", ms=8)
ax.annotate("Leidenfrost / minimum-heat-flux\npoint (rewetting temperature)", (bc.T_leid, bc.q_leid / 1e6),
            (470, 3.6), arrowprops=dict(arrowstyle="->"), fontsize=8)
ax.annotate("", xy=(140, 0.6), xytext=(600, 0.6), arrowprops=dict(arrowstyle="->", color="0.4", lw=1.2))
ax.text(370, 0.75, "path of a quenched surface (cooling)", ha="center", fontsize=8, color="0.3")
ax.set_xlim(20, 650)
ax.set_ylim(0, 7.2)
ax.set_xlabel("surface temperature T_s [°C]  (water, T_sat = 100 °C)")
ax.set_ylabel("heat flux q [MW/m²]")
ax.set_title("Boiling curve of a jet-impingement quench (Specht model, impingement zone)", fontsize=10)
fig.savefig(OUT + "boiling_curve_regimes.png", dpi=160)

# --- contact temperature vs Leidenfrost ------------------------------------
fig, ax = plt.subplots(figsize=(6.0, 3.4), constrained_layout=True)
T0 = np.linspace(200, 1000, 200)
for name, mat, c in (("AA6082", "aa6082", "#4cc3ff"), ("Nickel", "nickel", "#f26b1d"), ("C45", "c45", "#555")):
    m = qf.get_material(mat)
    e_s = effusivity(m.k(T0), m.rho(T0), m.cp(T0))
    ax.plot(T0, (e_s * T0 + e_w * 20) / (e_s + e_w), color=c, label=name)
ax.axhline(T_spiegler, ls="--", color="k", lw=0.9)
ax.text(210, T_spiegler + 12, "homogeneous-nucleation (Spiegler) limit of water ≈ %.0f °C" % T_spiegler, fontsize=8)
ax.plot(T0, T0, ":", color="0.6", lw=0.8)
ax.set_xlabel("wall temperature before contact T_w [°C]")
ax.set_ylabel("contact temperature T_i [°C]")
ax.set_title("Instantaneous solid–water contact temperature (water at 20 °C)")
ax.legend(fontsize=8)
fig.savefig(OUT + "contact_temperature.png", dpi=160)
