# 1. Conservation of energy and Fourier's law

<p class="bf-lead">Everything in BoilFlux rests on one balance law and one constitutive law. This chapter
derives the nonlinear heat equation of a quenched plate from them, states the boundary conditions of the
quenching problem, transforms it to the frame of a moving plate, and shows that the direct problem is
well posed.</p>

## 1.1 The continuum hypothesis and temperature

A metal plate contains about \(10^{29}\) atoms per cubic metre. We do not follow them. Instead we assign to
every point \(\mathbf x\) of the body \(\Omega\subset\mathbb R^3\) and every time \(t\) a **temperature**
\(T(\mathbf x,t)\): the local mean thermal energy of a volume that is large compared with the atomic spacing
(and with the electron and phonon mean free paths, below 100 nm in metals at these temperatures) but small
compared with the length scales we resolve (0.1 mm and above). This is the continuum hypothesis. It holds
everywhere in a quenched plate, so \(T\) can be treated as a smooth field.

## 1.2 The integral energy balance

Consider a fixed control volume \(V\subset\Omega\) with boundary \(\partial V\) and outward unit normal
\(\hat{\mathbf n}\). The solid is at rest in the frame we use (Section 1.6 treats a moving plate). The first law
of thermodynamics states that the internal energy in \(V\) changes only by heat crossing \(\partial V\) and by
volumetric sources:

\[
\frac{\mathrm d}{\mathrm dt}\int_V \rho\,e \,\mathrm dV
= -\oint_{\partial V} \mathbf q\cdot\hat{\mathbf n}\,\mathrm dS + \int_V \dot s\,\mathrm dV .
\tag{1.1}
\]

Here \(e\) is the specific internal energy [J/kg], \(\mathbf q\) the **heat-flux vector** [W/m²] and
\(\dot s\) a volumetric source [W/m³] (zero in quenching unless latent heat is released, Section 1.5). The minus
sign makes outflow reduce the energy. Mechanical work is negligible: thermal strains in a quenched plate are
below 1 %, and the work they do is many orders of magnitude smaller than the heat flows.

For a solid, \(\mathrm de = c_p\,\mathrm dT\) to excellent accuracy (the difference between \(c_p\) and \(c_v\) is
below 1 % for metals). We also neglect the change of the plate's geometry by thermal expansion and use the
volumetric heat capacity \(\rho(T)c_p(T)\) at the fixed grid volume. With \(V\) fixed in time, the derivative
can be taken inside the integral:

\[
\int_V \rho c_p\,\frac{\partial T}{\partial t}\,\mathrm dV
= -\oint_{\partial V} \mathbf q\cdot\hat{\mathbf n}\,\mathrm dS + \int_V \dot s\,\mathrm dV .
\tag{1.2}
\]

Equation (1.2) is exactly the form that the finite-volume method of [Chapter 6](06-finite-volumes.md)
discretises: each grid cell is a control volume.

## 1.3 Fourier's law

The balance (1.2) has two unknown fields, \(T\) and \(\mathbf q\). A constitutive law closes it.

!!! definition "Definition 1.1 (Fourier's law)"
    In an isotropic solid the heat flux is proportional to the temperature gradient,

    \[
    \mathbf q = -k(T)\,\nabla T ,
    \tag{1.3}
    \]

    with the **thermal conductivity** \(k>0\) [W/(m K)]. In an anisotropic solid \(\mathbf q=-\mathbf K\nabla T\)
    with a symmetric conductivity tensor \(\mathbf K\).

Fourier's law is empirical, but its sign is fixed by the second law. The local entropy production rate per
unit volume in pure conduction is

\[
\dot\sigma = \mathbf q\cdot\nabla\!\left(\frac1T\right) = -\frac{\mathbf q\cdot\nabla T}{T^2}
= \frac{\nabla T\cdot\mathbf K\nabla T}{T^2}\ \ge 0 ,
\tag{1.4}
\]

which must hold for every gradient. Hence \(\mathbf K\) is positive semidefinite (and positive definite for a
real conductor), and in the isotropic case \(k>0\): **heat flows from hot to cold**. The polycrystalline metals
treated here are isotropic at the scale of the grid.

## 1.4 The heat equation

Apply the divergence theorem ([Appendix B](appendix-b-vector-calculus.md), Eq. B.1) to the surface integral in
(1.2) and insert (1.3):

\[
\int_V \Big[\rho c_p\,\frac{\partial T}{\partial t} - \nabla\cdot\big(k\nabla T\big) - \dot s\Big]\mathrm dV = 0 .
\]

This holds for **every** control volume \(V\subset\Omega\). If the integrand is continuous and non-zero at some
point, it keeps its sign in a small ball around that point, and the integral over that ball cannot vanish.
The integrand is therefore zero everywhere:

!!! theorem "Result 1.1 (nonlinear heat equation)"
    \[
    \rho(T)\,c_p(T)\,\frac{\partial T}{\partial t} = \nabla\cdot\big(k(T)\,\nabla T\big) + \dot s
    \qquad\text{in } \Omega\times(0,t_f].
    \tag{1.5}
    \]

Expanding the divergence with the chain rule, \(\nabla k(T) = k'(T)\nabla T\), shows where the nonlinearity
sits:

\[
\rho c_p\,\frac{\partial T}{\partial t} = k\,\Delta T + k'(T)\,\lvert\nabla T\rvert^2 + \dot s .
\tag{1.6}
\]

For constant properties (1.5) reduces to the linear heat equation
\(\partial_t T = \alpha\,\Delta T\) with the **thermal diffusivity**

\[
\alpha = \frac{k}{\rho c_p}\quad[\mathrm{m^2/s}] .
\tag{1.7}
\]

For AA6082, \(\alpha\approx 7.2\times10^{-5}\) m²/s; for nickel near 800 °C, \(\alpha\approx 1.4\times10^{-5}\) m²/s.
The diffusivity sets every time scale in this book.

## 1.5 Two useful transformations

### Kirchhoff transformation

Define the **Kirchhoff temperature**

\[
U(T) = \frac{1}{k_0}\int_{T_\text{ref}}^{T} k(\tau)\,\mathrm d\tau ,
\tag{1.8}
\]

with a reference conductivity \(k_0\). Since \(\nabla U = (k/k_0)\nabla T\),

\[
\nabla\cdot(k\nabla T) = k_0\,\Delta U .
\]

Similarly \(\partial_t U = (k/k_0)\,\partial_t T\), so \(\rho c_p\,\partial_t T = (k_0/\alpha)\,\partial_t U\), and (1.5) becomes

\[
\frac{k_0}{\alpha(T)}\,\frac{\partial U}{\partial t} = k_0\,\Delta U + \dot s ,
\qquad \alpha(T) = \frac{k(T)}{\rho(T)c_p(T)} .
\tag{1.9}
\]

The diffusion operator is now linear; the nonlinearity has moved into the coefficient \(1/\alpha(T)\), and the
flux boundary condition \(-k\,\partial_n T = q\) becomes the linear condition \(-k_0\,\partial_n U = q\). We use
this in [Chapter 10](10-adjoint.md) to linearise the problem exactly.

### Enthalpy form and latent heat

Steels release latent heat during solid-state transformations (austenite to ferrite, pearlite, bainite or
martensite), and nickel shows a peak in \(c_p\) at its Curie point (about 358 °C). With the volumetric
enthalpy

\[
H(T) = \int_{T_\text{ref}}^{T}\rho c_p(\tau)\,\mathrm d\tau ,
\qquad \frac{\partial H}{\partial t} = \nabla\cdot(k\nabla T),
\tag{1.10}
\]

a transformation that releases latent heat \(\mathcal L\) over a temperature interval appears as a peak in the
**effective** \(\rho c_p = \mathrm dH/\mathrm dT\). This is how BoilFlux handles it: the property table carries
the peak. It is exact when the transformation is a function of temperature only. Transformations that depend
on cooling rate (bainite, martensite start) need a kinetic model and are outside the current scope.

## 1.6 Boundary and initial conditions

The heat equation needs one condition at every point of the boundary and an initial temperature. With
outward normal \(\hat{\mathbf n}\) and \(\partial_n T = \nabla T\cdot\hat{\mathbf n}\), the outward heat flux through the
boundary is \(\mathbf q\cdot\hat{\mathbf n} = -k\,\partial_n T\).

| Type | Condition | Physical meaning |
|---|---|---|
| Dirichlet | \(T = T_b\) | temperature imposed |
| Neumann | \(-k\,\partial_n T = q\) | outward heat flux imposed (\(q=0\): insulated) |
| Robin | \(-k\,\partial_n T = h\,(T-T_\infty)\) | convection to a fluid at \(T_\infty\) |
| Radiation | \(-k\,\partial_n T = \varepsilon\sigma_\text{SB}(T^4-T_\infty^4)\) | grey-body exchange with surroundings (kelvin) |
| Boiling | \(-k\,\partial_n T = q_b(T_s;\ \mathbf x, t)\) | nonlinear boiling curve (Chapter 4) |

Radiation is linearised exactly by factorising the difference of fourth powers,

\[
T^4 - T_\infty^4 = (T^2+T_\infty^2)(T+T_\infty)(T-T_\infty)
\quad\Longrightarrow\quad
h_\text{rad} = \varepsilon\sigma_\text{SB}\,(T^2+T_\infty^2)(T+T_\infty) .
\tag{1.11}
\]

BoilFlux instead linearises the total loss about the previous time level (Section 7.3), which is also exact in
the limit \(\Delta t\to0\).

!!! definition "Definition 1.2 (the direct quenching problem)"
    The plate \(\Omega=(0,L_x)\times(0,L_y)\times(0,L_z)\) has the quenched face \(\Gamma_j=\{z=0\}\), the
    infrared-observed face \(\Gamma_c=\{z=L_z\}\) and four edges \(\Gamma_e\). Given the initial temperature
    \(T_0(\mathbf x)\) and the surface heat flux \(q(x,y,t)\) on \(\Gamma_j\), find \(T\) such that

    \[
    \begin{aligned}
    \rho c_p(T)\,\partial_t T &= \nabla\cdot\big(k(T)\nabla T\big) && \text{in }\Omega\times(0,t_f],\\
    -k(T)\,\partial_n T &= q(x,y,t) && \text{on }\Gamma_j,\\
    -k(T)\,\partial_n T &= q_\text{loss}(T) && \text{on }\Gamma_c,\\
    -k(T)\,\partial_n T &= 0 && \text{on }\Gamma_e,\\
    T(\mathbf x,0) &= T_0(\mathbf x) && \text{in }\Omega .
    \end{aligned}
    \tag{1.12}
    \]

    On \(\Gamma_j\) the outward normal is \(-\hat{\mathbf z}\), so the condition reads \(k\,\partial_z T = q\) at
    \(z=0\): a positive \(q\) extracts heat and makes the plate colder at the quenched face than inside.

In the **inverse** problem (Part IV), \(q\) is unknown and the temperature on \(\Gamma_c\) is measured instead.
Problem (1.12) is then the *model* that links them.

## 1.7 Moving plates and the choice of frame

In a moving-plate experiment the jet and the camera are fixed in the laboratory and the plate moves along
\(x\) with velocity \(v\). Let \(X\) be the laboratory coordinate and \(x\) the coordinate attached to the
plate, with

\[
X = X_0 + v\,t + x .
\tag{1.13}
\]

A material point has fixed \(x\). Let \(\tilde T(X,y,z,t)=T(x,y,z,t)\) be the same temperature field described in
laboratory coordinates. By the chain rule, at fixed \(x\) (following the material),

\[
\left.\frac{\partial T}{\partial t}\right|_x
= \left.\frac{\partial \tilde T}{\partial t}\right|_X + v\,\frac{\partial \tilde T}{\partial X} .
\tag{1.14}
\]

The energy balance (1.2) applies to a material volume, which is at rest in the plate frame. In the plate frame
the equation is therefore (1.5) unchanged, while in the laboratory frame it acquires an advection term:

\[
\text{plate frame: }\ \rho c_p\,\partial_t T = \nabla\cdot(k\nabla T),
\qquad
\text{lab frame: }\ \rho c_p\big(\partial_t \tilde T + v\,\partial_X \tilde T\big) = \nabla\cdot(k\nabla\tilde T).
\tag{1.15}
\]

The physics is identical; only the description differs. The price of the plate frame is that the jet, fixed at
\(X_\text{jet}\), moves across the plate:

\[
x_s(t) = X_\text{jet} - X_0 - v\,t ,
\tag{1.16}
\]

and the camera window becomes a time-dependent set of observed plate points. BoilFlux uses the plate frame
because the discrete operator then stays symmetric (no upwinding of the advection term), the plate's real
edges are boundaries of the domain, and the adjoint of Chapter 10 is exact.

!!! note "In BoilFlux"
    `boilflux.kinematics.PlateMotion` implements (1.13) and (1.16); `lab_frames_to_plate` resamples camera
    frames to plate coordinates; `PlateMotion.view_mask` gives the observed set.

## 1.8 Well-posedness of the direct problem

A problem is **well posed** in the sense of Hadamard if a solution exists, is unique and depends
continuously on the data. The direct problem (1.12) is well posed; the inverse problem of Part IV is not.
Here is the uniqueness and stability argument for constant properties. The nonlinear case follows with the
Kirchhoff variable of Section 1.5 when \(\alpha\) is bounded above and below.

!!! theorem "Theorem 1.2 (energy estimate)"
    Let \(T_1, T_2\) solve (1.12) with constant \(k,\rho c_p\), no loss on \(\Gamma_c\), initial data
    \(T_{0,1}, T_{0,2}\) and fluxes \(q_1, q_2\). Then \(w = T_1 - T_2\) satisfies, for all \(t\le t_f\),

    \[
    \rho c_p\,\lVert w(t)\rVert^2 \le \rho c_p\,\lVert w(0)\rVert^2 + \frac{1}{\epsilon}\int_0^t\!\!\int_{\Gamma_j}(q_1-q_2)^2\,\mathrm dS\,\mathrm dt'
    + \epsilon\int_0^t\!\!\int_{\Gamma_j} w^2\,\mathrm dS\,\mathrm dt'
    \tag{1.17}
    \]

    for any \(\epsilon>0\), with \(\lVert\cdot\rVert\) the \(L^2(\Omega)\) norm. In particular, equal data give
    \(w\equiv0\): the solution is unique.

??? derivation "Proof"
    The difference \(w\) satisfies the linear heat equation with flux \(\delta q = q_1-q_2\) on \(\Gamma_j\) and
    zero flux elsewhere. Multiply by \(w\) and integrate over \(\Omega\):

    \[
    \frac{\rho c_p}{2}\frac{\mathrm d}{\mathrm dt}\int_\Omega w^2\,\mathrm dV
    = \int_\Omega w\,\nabla\cdot(k\nabla w)\,\mathrm dV
    = \oint_{\partial\Omega} w\,k\,\partial_n w\,\mathrm dS - \int_\Omega k\lvert\nabla w\rvert^2\,\mathrm dV ,
    \]

    using Green's first identity (B.3). On \(\Gamma_j\), \(k\,\partial_n w = -\delta q\); elsewhere it vanishes. Drop
    the non-positive volume term and apply Young's inequality \(\lvert ab\rvert\le a^2/(2\epsilon)+\epsilon b^2/2\):

    \[
    \frac{\rho c_p}{2}\frac{\mathrm d}{\mathrm dt}\lVert w\rVert^2 \le \int_{\Gamma_j}\lvert w\,\delta q\rvert\,\mathrm dS
    \le \frac{1}{2\epsilon}\int_{\Gamma_j}\delta q^2\,\mathrm dS + \frac{\epsilon}{2}\int_{\Gamma_j} w^2\,\mathrm dS .
    \]

    Integrate from 0 to \(t\) and multiply by 2. With equal data, \(\delta q = 0\), so
    \(\lVert w(t)\rVert\le\lVert w(0)\rVert = 0\). \(\square\)

The estimate says that small changes of the flux produce small changes of the temperature. The converse,
which is what the inverse problem needs, is false: Chapter 2 shows that very large, rapidly varying changes
of \(q\) produce almost no change at the opposite face.

A second property, used to check the discretisation in Chapter 6, is the **maximum principle**: for the
linear heat equation without sources and with insulated boundaries, \(\max_\Omega T(\cdot,t)\) does not increase
and \(\min_\Omega T(\cdot,t)\) does not decrease in time. A quench (\(q\ge0\) on \(\Gamma_j\)) can therefore never
heat any point of the plate above its initial maximum.

## 1.9 The global energy balance

Integrating (1.5) over the whole plate and using the boundary conditions of (1.12) gives

\[
\frac{\mathrm d}{\mathrm dt}\int_\Omega H(T)\,\mathrm dV
= -\int_{\Gamma_j} q\,\mathrm dS - \int_{\Gamma_c} q_\text{loss}\,\mathrm dS .
\tag{1.18}
\]

The total heat extracted by the jet up to time \(t\) equals the drop in the plate's enthalpy minus the loss on
the IR face. This is a useful sanity check of any reconstruction: the integral of the recovered \(q\) must match
the enthalpy drop computed from the measured temperatures, within the measurement error.

## Exercises

!!! exercise "1.1"
    Show that the steady one-dimensional solution of (1.5) with \(k(T)=k_0(1+\beta T)\) between \(T(0)=T_a\) and
    \(T(L)=T_b\) is linear in the Kirchhoff variable \(U\), and find \(T(z)\). For which sign of \(\beta\) does the
    profile bulge above the straight line?

!!! exercise "1.2"
    A 5 mm AA6082 plate (\(\rho c_p = 2.5\) MJ/(m³ K)) cools from 500 °C to 100 °C in 4 s under a uniform flux.
    Use (1.18) to estimate the mean flux. *Answer:* \(q\approx 1.25\) MW/m².

!!! exercise "1.3"
    Repeat the derivation of (1.14) for a plate whose position \(X_0 + s(t)\) follows a measured encoder signal
    \(s(t)\). What replaces (1.16)?
