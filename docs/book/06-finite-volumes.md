# 6. Finite-volume discretisation

<p class="bf-lead">This chapter turns the energy balance (1.2) into the sparse linear system that BoilFlux
solves at every time step. Each step of the derivation corresponds to a line of
<code>boilflux.forward</code>. The chapter ends with three properties the inverse solver relies on: the system
matrix is symmetric positive definite, it obeys a discrete maximum principle for any time step, and it
conserves energy exactly.</p>

## 6.1 Grid and control volumes

The plate is divided into \(n_x\times n_y\times n_z\) equal boxes of size \(\Delta x\times\Delta y\times\Delta z\), with
\(\Delta x = L_x/n_x\) and so on. Each box is a **control volume** (cell) with its unknown temperature \(T_P\) stored at
the centre. Cells are numbered with \(x\) fastest, then \(y\), then \(z\):

\[
P = (i_z\,n_y + i_y)\,n_x + i_x ,\qquad V = \Delta x\,\Delta y\,\Delta z .
\tag{6.1}
\]

The first layer \(i_z = 0\) touches the quenched face \(\Gamma_j\), the last \(i_z = n_z-1\) the IR face \(\Gamma_c\). Each face
of a cell has area \(A_x = \Delta y\Delta z\), \(A_y = \Delta x\Delta z\) or \(A_z = \Delta x\Delta y\). Setting \(n_y = 1\) gives the 2D
length × thickness model of Gopalkrishna et al. (2026).

!!! note "In BoilFlux"
    `boilflux.grid.PlateGrid` defines the grid; `PlateGrid.pairs` lists every interior face as a pair of cells
    \((P,N)\) with its geometric factor \(A/\Delta\).

## 6.2 Integrating the balance over a cell

Apply (1.2) to cell \(P\), with no source. The surface integral splits into a sum over the six faces \(f\):

\[
\int_{V_P}\rho c_p\,\frac{\partial T}{\partial t}\,\mathrm dV = -\sum_{f\in\partial P} F_f ,
\qquad F_f = \int_f \mathbf q\cdot\hat{\mathbf n}_f\,\mathrm dS ,
\tag{6.2}
\]

where \(F_f\) is the heat leaving \(P\) through face \(f\) [W]. This is still exact. Two approximations follow: the
volume integral is approximated with the cell-centre value, \(\rho c_p(T_P)\,V\,\mathrm dT_P/\mathrm dt\), and each face flux
is expressed through the neighbouring cell temperatures.

## 6.3 Interior faces: the harmonic mean

Take two neighbouring cells \(P\) and \(N\) a distance \(\Delta\) apart, with conductivities \(k_P\) and \(k_N\) evaluated at
their own temperatures. Assume the temperature is linear within each half-cell, with face temperature \(T_f\).
The flux leaving \(P\) through the face must equal the flux entering \(N\) (no heat accumulates in a face):

\[
F = k_P A\,\frac{T_P - T_f}{\Delta/2} = k_N A\,\frac{T_f - T_N}{\Delta/2} .
\]

Eliminating the unknown \(T_f\) is the familiar addition of thermal resistances in series:

\[
F = \frac{T_P - T_N}{\dfrac{\Delta}{2k_P A} + \dfrac{\Delta}{2k_N A}}
= G_{PN}\,(T_P - T_N),
\qquad
G_{PN} = \frac{A}{\Delta}\cdot\frac{2k_Pk_N}{k_P + k_N} .
\tag{6.3}
\]

The face conductance \(G_{PN}\) [W/K] uses the **harmonic mean** of the cell conductivities. It is exact for a
layered material with a jump of \(k\) at the face, while the arithmetic mean would overestimate the flux when the
contrast is large. For smoothly varying \(k(T)\) both means agree to second order.

## 6.4 Boundary faces

**Quenched face (prescribed flux).** For a cell in the first layer, the face on \(\Gamma_j\) carries the outward heat
flow \(F = q\,A_z\) by the boundary condition of (1.12). It enters the equation as a known term.

**IR face (loss).** On \(\Gamma_c\) the outward flow is \(q_\text{loss}(T)\,A_z\), linearised in Section 7.3. Without loss,
\(F=0\).

**Edges.** Insulated: \(F = 0\).

**Surface temperature.** The temperature of the quenched face itself is not a cell unknown. With the linear
half-cell profile and \(k\,\partial_z T = q\) at \(z=0\),

\[
T_s \approx T_P - \frac{\Delta z}{2}\,\frac{q}{k(T_P)} .
\tag{6.4}
\]

**Observation at the IR face.** The camera measures the face temperature \(T(L_z)\), while the model provides the
cell-centre value at \(z = L_z - \Delta z/2\). A Taylor expansion about the face gives

\[
T_P = T(L_z) - \frac{\Delta z}{2}\,\partial_z T(L_z) + \frac{\Delta z^2}{8}\,\partial_{zz}T(L_z) + O(\Delta z^3).
\tag{6.5}
\]

On an insulated face \(\partial_z T = 0\), so the model error of using \(T_P\) is only \(O(\Delta z^2)\). With a loss \(q_\text{loss}\),
the first-order term is \(\Delta z\,q_\text{loss}/(2k)\): for 1 mm cells, \(k = 180\) W/(m K) and a loss of 20 kW/m² it is
0.06 K, well below camera noise.

!!! note "In BoilFlux"
    `AdjointCGM.surface_temperature` implements (6.4); `Trajectory.top_face` returns the cell-centre values
    \(T_P\) of the last layer as the model of the IR data.

## 6.5 The semi-discrete system

Collect all cell temperatures in the vector \(\mathbf T\in\mathbb R^{N}\), \(N = n_xn_yn_z\), and the face fluxes of the quenched
face in \(\mathbf q\in\mathbb R^{n_xn_y}\). Summing (6.3) over all interior faces, (6.2) becomes the system of ordinary
differential equations

\[
\mathbf M(\mathbf T)\,\frac{\mathrm d\mathbf T}{\mathrm dt} = -\mathbf K(\mathbf T)\,\mathbf T - \mathbf B\,\mathbf q - \mathbf f_\text{loss}(\mathbf T) ,
\tag{6.6}
\]

with the diagonal **capacity matrix** \(\mathbf M = \operatorname{diag}(\rho c_p(T_P)V)\), the **conduction matrix**

\[
\mathbf K = \sum_{\text{faces }(P,N)} G_{PN}\,(\mathbf e_P - \mathbf e_N)(\mathbf e_P - \mathbf e_N)^\top ,
\tag{6.7}
\]

the **flux operator** \(\mathbf B\in\mathbb R^{N\times n_xn_y}\), which places \(A_z\) in the rows of the first-layer cells, and the loss
vector \(\mathbf f_\text{loss}\), non-zero in the last layer. Here \(\mathbf e_P\) is the \(P\)-th unit vector. In components, \(K_{PP} = \sum_N G_{PN}\)
and \(K_{PN} = -G_{PN}\) for neighbours.

!!! theorem "Lemma 6.1 (properties of the conduction matrix)"
    For any conductances \(G_{PN}>0\):

    1. \(\mathbf K\) is symmetric;
    2. \(\mathbf x^\top\mathbf K\mathbf x = \sum_{(P,N)} G_{PN}(x_P - x_N)^2 \ge 0\), so \(\mathbf K\) is positive semidefinite;
    3. \(\mathbf K\mathbf 1 = \mathbf 0\) and \(\mathbf 1^\top\mathbf K = \mathbf 0^\top\): interior conduction moves heat but neither creates nor destroys it.

??? derivation "Proof"
    Each term of (6.7) is a symmetric rank-one matrix with a positive weight, which proves 1. For 2, note
    \(\mathbf x^\top(\mathbf e_P - \mathbf e_N)(\mathbf e_P - \mathbf e_N)^\top\mathbf x = (x_P - x_N)^2\). For 3, \((\mathbf e_P - \mathbf e_N)^\top\mathbf 1 = 0\) for every
    face. \(\square\)

\(\mathbf K\) is the weighted **graph Laplacian** of the grid. This structure is why the adjoint of Chapter 10 needs no
transposition work: \(\mathbf K^\top = \mathbf K\).

## 6.6 Accuracy

For smooth solutions on a uniform grid, the two-point flux (6.3) is a centred difference, and the cell-centre
quadrature of the volume integral is second-order accurate. The local truncation error of the interior equations
is therefore \(O(\Delta x^2+\Delta y^2+\Delta z^2)\). At a flux boundary the half-cell approximation is locally first order,
but because the boundary layer of cells has measure \(O(\Delta z)\), the global error remains \(O(\Delta^2)\) in space. Time
integration (Chapter 7) adds \(O(\Delta t)\).

Figure 2.1 is the corresponding verification: with 40 cells through 5 mm and \(\Delta t = t_d/2000\), BoilFlux reproduces
the exact slab solution (2.14) to 0.014 K on a scale of 28 K.

## Exercises

!!! exercise "6.1"
    Show that the arithmetic mean \((k_P+k_N)/2\) is never smaller than the harmonic mean, and that for
    \(k_P/k_N = 10\) it overestimates the flux across a material interface by a factor of 3.

!!! exercise "6.2"
    Write \(\mathbf K\) explicitly for a 1D grid of 4 cells with equal conductances and check Lemma 6.1.

!!! exercise "6.3"
    Derive a second-order surface-temperature formula that uses the first two cell layers instead of (6.4).
