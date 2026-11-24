# Appendix A. Notation and units

## Symbols

| Symbol | Meaning | Unit |
|---|---|---|
| \(T\), \(T_s\), \(T_w\) | temperature; surface (quenched-face) temperature; wall temperature | K or °C |
| \(T_0\) | initial temperature (Part I–IV); rewetting temperature (Chapter 5 only) | K or °C |
| \(T_\text{sat}\), \(T_\text{DNB}\), \(T_\text{Leid}\) | saturation, maximum-heat-flux and Leidenfrost (rewetting) temperatures | K or °C |
| \(q\) | heat flux leaving the plate through the quenched face (positive = cooling) | W/m² |
| \(\mathbf q\) | heat-flux vector (Chapter 1); discrete flux unknowns (Parts III–IV) | W/m² |
| \(\mathbf p = \mathbf q/Q\) | scaled optimisation variable, \(Q = 10^6\) W/m² | – |
| \(k\), \(\rho\), \(c_p\) | conductivity, density, specific heat capacity | W/(m K), kg/m³, J/(kg K) |
| \(\alpha = k/(\rho c_p)\) | thermal diffusivity | m²/s |
| \(e = \sqrt{k\rho c_p}\) | thermal effusivity | J/(m² K s\(^{1/2}\)) |
| \(H(T)\) | volumetric enthalpy | J/m³ |
| \(h\) | heat-transfer coefficient | W/(m² K) |
| \(h_{fg}\), \(h'_{fg}\) | latent heat of vaporisation; effective latent heat | J/kg |
| \(\sigma\) | surface tension (Part II); noise standard deviation (Parts IV–VI) | N/m; K |
| \(\sigma_\text{SB}\) | Stefan–Boltzmann constant, \(5.670\times10^{-8}\) | W/(m² K⁴) |
| \(\varepsilon\) | emissivity | – |
| \(L, L_x, L_y, L_z\) | plate thickness; plate dimensions | m |
| \(L_c\) | capillary length \(\sqrt{\sigma/(g(\rho_l-\rho_v))}\) | m |
| \(\Omega\), \(\Gamma_j\), \(\Gamma_c\), \(\Gamma_e\) | plate; quenched face \(z=0\); IR face \(z=L_z\); edges | – |
| \(\hat{\mathbf n}\), \(\partial_n\) | outward unit normal; normal derivative | – |
| \(v\) | plate velocity | m/s |
| \(\mathrm{Fo}, \mathrm{Bi}, \mathrm{Pe}\) | Fourier, Biot and Péclet numbers, Eq. (2.1) | – |
| \(t_d = L^2/\alpha\) | through-thickness diffusion time | s |
| \(\mathbf A_n, \mathbf C_n, \mathbf K, \mathbf B, \mathbf H, \mathbf M_n\) | step matrix, capacity/\(\Delta t\), conduction, flux, observation and mask operators | – |
| \(\boldsymbol\theta\), \(\boldsymbol\mu\), \(\lambda\) | sensitivity; discrete adjoint; continuous adjoint | – |
| \(\mathbf X\), \(s_i, \mathbf u_i, \mathbf v_i\) | sensitivity matrix; its singular values and vectors | – |
| \(\gamma\), \(\tau\), \(\delta\) | Tikhonov weight; discrepancy safety factor; noise norm \(\sqrt m\,\sigma\) | – |
| \(f_i\) | filter factor | – |
| \(\theta\) (Part V) | network parameters | – |

## Water at 1 atm (saturation, 100 °C)

Values from standard steam tables, as used in Chapters 3–4 and in `docs/scripts/make_physics_figures.py`.

| Property | Liquid | Vapour |
|---|---:|---:|
| density \(\rho\) [kg/m³] | 957.9 | 0.596 |
| specific heat \(c_p\) [J/(kg K)] | 4217 | 2029 |
| conductivity \(k\) [W/(m K)] | 0.680 | 0.0248 |
| viscosity \(\mu\) [Pa s] | \(2.79\times10^{-4}\) | \(1.20\times10^{-5}\) |
| Prandtl number | 1.76 | 1.00 |

Latent heat \(h_{fg} = 2.257\times10^6\) J/kg; surface tension \(\sigma = 0.0589\) N/m; critical point \(T_c = 647.1\) K,
\(p_c = 22.06\) MPa. Derived: capillary length 2.5 mm; most dangerous Taylor wavelength 27 mm.

## Representative solid properties

| Material | \(T\) [°C] | \(k\) [W/(m K)] | \(\rho\) [kg/m³] | \(c_p\) [J/(kg K)] | \(\alpha\) [m²/s] |
|---|---:|---:|---:|---:|---:|
| AA6082 | 20 | 172 | 2700 | 894 | \(7.1\times10^{-5}\) |
| AA6082 | 500 | 191 | 2610 | 1120 | \(6.5\times10^{-5}\) |
| Nickel | 20 | 90.7 | 8900 | 445 | \(2.3\times10^{-5}\) |
| Nickel | 800 | 69.0 | 8670 | 555 | \(1.4\times10^{-5}\) |

These are the BoilFlux tables in `boilflux.materials` (representative literature values). Replace them with measured
properties for quantitative work. Chapter 2 uses the constant set \(k = 180\) W/(m K), \(\rho = 2700\) kg/m³,
\(c_p = 926\) J/(kg K) (\(\alpha = 7.2\times10^{-5}\) m²/s) as a representative AA6082 value.
