# 7. Time integration, nonlinearity and verification

<p class="bf-lead">The semi-discrete system (6.6) is stiff: its fastest modes decay thousands of times faster
than the ones we observe. This chapter explains why BoilFlux integrates it with the backward Euler method,
how the temperature-dependent properties and boundary laws are linearised so that each step is one symmetric
linear solve, and how the result is verified.</p>

## 7.1 The backward Euler step

Integrate (6.6) from \(t_{n-1}\) to \(t_n = t_{n-1}+\Delta t\) and evaluate the right-hand side at the new time level
(implicit, or backward, Euler):

\[
\mathbf M\,\frac{\mathbf T^n - \mathbf T^{n-1}}{\Delta t} = -\mathbf K\,\mathbf T^n - \mathbf B\,\mathbf q^n - \mathbf f_\text{loss}(\mathbf T^n) .
\tag{7.1}
\]

The flux \(\mathbf q^n\) is the flux applied during the step \((t_{n-1}, t_n]\). The local truncation error is \(O(\Delta t^2)\) per step,
so the global error is \(O(\Delta t)\).

## 7.2 Why an implicit method: stability and stiffness

Consider the 1D heat equation with constant properties on a uniform grid, and a Fourier mode
\(T_j^n = g^n e^{\mathrm i j\varphi}\) with \(\varphi\in[0,\pi]\). With the mesh ratio \(r = \alpha\Delta t/\Delta x^2\) and
\(s = \sin^2(\varphi/2)\), the three classical one-step methods have the amplification factors

\[
g_\text{FE} = 1 - 4rs,\qquad
g_\text{BE} = \frac{1}{1 + 4rs},\qquad
g_\text{CN} = \frac{1 - 2rs}{1 + 2rs} .
\tag{7.2}
\]

!!! theorem "Result 7.1 (stability of the three methods)"
    - Forward Euler is stable (\(\lvert g\rvert\le1\) for all \(\varphi\)) only if \(r\le\tfrac12\) in 1D (in 3D,
      \(\alpha\Delta t\sum_d \Delta_d^{-2}\le\tfrac12\)).
    - Backward Euler satisfies \(0 < g \le 1\) for every \(r\): it is unconditionally stable and never oscillates.
    - Crank–Nicolson satisfies \(\lvert g\rvert\le1\) for every \(r\), but \(g\to-1\) as \(rs\to\infty\): stiff modes are
      barely damped and flip sign every step.

For the benchmark grid (\(\Delta z = 1\) mm, AA6082), forward Euler would need \(\Delta t\le 7\) ms in the thickness
direction alone, and less in 3D. A quench has sharp changes of \(q\) at rewetting; Crank–Nicolson turns them
into persistent oscillations of the first cell layer, which the inverse solver would interpret as flux. Backward
Euler damps them monotonically (it is **L-stable**: \(g\to0\) for stiff modes), and its discrete adjoint is again a
backward Euler step (Chapter 10). Its first-order accuracy in time is adequate because the time step is set by
the camera frame rate and the flux is piecewise constant between frames.

## 7.3 Linearising the nonlinearity

Equation (7.1) is nonlinear in \(\mathbf T^n\): \(\mathbf M\) and \(\mathbf K\) depend on temperature through \(\rho c_p(T)\) and \(k(T)\), and the
loss through \(T^4\). BoilFlux uses a **semi-implicit** linearisation that keeps one linear solve per step.

**Lagged properties.** Evaluate \(\mathbf M\) and \(\mathbf K\) at the known \(\mathbf T^{n-1}\). The error is
\(O(\Delta t\,\lVert\partial_t T\rVert\,\lvert\mathrm d\ln k/\mathrm dT\rvert)\) per step, of the same order as the time discretisation.

**Linearised loss.** Expand the loss to first order about \(\mathbf T^{n-1}\):

\[
q_\text{loss}(T^n) \approx q_\text{loss}(T^{n-1}) + q'_\text{loss}(T^{n-1})\,(T^n - T^{n-1}) .
\tag{7.3}
\]

The part proportional to \(T^n\) moves to the left-hand side as a diagonal matrix \(\mathbf L_n\) with entries
\(q'_\text{loss}A_z\) in the last layer; the rest is a known vector \(\mathbf c_n\) with entries
\(-(q_\text{loss} - q'_\text{loss}T^{n-1})A_z\).

!!! theorem "Result 7.2 (one BoilFlux time step)"
    \[
    \mathbf A_n\mathbf T^n = \mathbf C_n\mathbf T^{n-1} - \mathbf B\,\mathbf q^n + \mathbf c_n,
    \qquad
    \mathbf A_n = \mathbf C_n + \mathbf K(\mathbf T^{n-1}) + \mathbf L_n,
    \qquad
    \mathbf C_n = \frac{\mathbf M(\mathbf T^{n-1})}{\Delta t} .
    \tag{7.4}
    \]

    \(\mathbf A_n\) is symmetric and positive definite: \(\mathbf C_n\) is diagonal and positive, \(\mathbf K\) is positive semidefinite
    (Lemma 6.1), and \(\mathbf L_n\) is diagonal and non-negative because \(q'_\text{loss} = h + 4\varepsilon\sigma_\text{SB}T^3>0\).

The step is **affine** in \(\mathbf T^n\), \(\mathbf T^{n-1}\) and \(\mathbf q^n\) once the coefficients are frozen. This is what makes exact
sensitivity and adjoint problems available in Chapter 10.

**Boiling boundary law (synthetic data only).** When the flux is a law \(q_b(T)\) instead of a prescribed field,
the same linearisation is applied on the quenched face,
\(q^n \approx q_b(T^{n-1}_P) + q_b'(T^{n-1}_P)\,(T^n_P - T^{n-1}_P)\), with the derivative (4.13). The applied flux is
recorded and becomes the ground truth.

!!! warning "Which temperature the synthetic boiling law sees"
    `JetQuenchScenario` evaluates the boiling curve at the first-layer **cell** temperature \(T_P\), not at the
    surface temperature (6.4). The two differ by \(q\,\Delta z/(2k)\): at 6.5 MW/m² on the refined grid
    (\(\Delta z = 0.5\) mm) about 9 K for AA6082 and about 20 K for nickel. The ground truth is self-consistent, since
    the recorded \(q\) is the flux actually applied. But when the surface temperature at which a reconstruction
    peaks is compared with the curve's nominal \(T_\text{DNB}\), part of the difference comes from this offset rather
    than from the inversion.

## 7.4 Discrete conservation and the maximum principle

!!! theorem "Theorem 7.3 (discrete energy balance)"
    The step (7.4) satisfies exactly

    \[
    \sum_P \frac{\rho c_p(T_P^{n-1})\,V}{\Delta t}\,\big(T_P^n - T_P^{n-1}\big)
    = -A_z\sum_{\Gamma_j} q^n - \sum_{\Gamma_c}\big(\text{linearised loss}\big)\,A_z .
    \tag{7.5}
    \]

??? derivation "Proof"
    Multiply (7.4) from the left by \(\mathbf 1^\top\). By Lemma 6.1, \(\mathbf 1^\top\mathbf K = \mathbf 0^\top\), so the conduction terms drop out.
    \(\mathbf 1^\top\mathbf B\mathbf q^n = A_z\sum q^n\), and the loss terms combine to the linearised loss. \(\square\)

Equation (7.5) is the discrete form of (1.18). The heat removed through the quenched face is accounted for to
round-off; the only approximation is that the enthalpy change is evaluated with the lagged heat capacity, which
differs from \(\int\rho c_p\,\mathrm dT\) by \(O(\lvert\Delta T\rvert^2)\) per step.

!!! theorem "Theorem 7.4 (discrete maximum principle, any \(\Delta t\))"
    Without loss and with \(\mathbf q^n\ge0\) (extraction), every step (7.4) satisfies
    \(\max_P T_P^n \le \max_P T_P^{n-1}\). With \(\mathbf q^n = 0\) it also satisfies \(\min_P T^n_P\ge\min_P T^{n-1}_P\).

??? derivation "Proof"
    \(\mathbf A_n\) has positive diagonal entries, non-positive off-diagonal entries, and is strictly diagonally dominant
    because \(A_{PP} - \sum_{N\ne P}\lvert A_{PN}\rvert = (\mathbf C_n)_{PP} > 0\). Such a matrix is a non-singular **M-matrix**, and its
    inverse has non-negative entries ([Appendix D](appendix-d-linear-algebra.md)). Write the step as
    \(\mathbf T^n = \mathbf A_n^{-1}\mathbf C_n\mathbf T^{n-1} - \mathbf A_n^{-1}\mathbf B\mathbf q^n\). The matrix \(\mathbf P = \mathbf A_n^{-1}\mathbf C_n\) is non-negative, and since
    \(\mathbf A_n\mathbf 1 = \mathbf C_n\mathbf 1\) (rows of \(\mathbf K\) sum to zero) we have \(\mathbf P\mathbf 1 = \mathbf 1\): each \(T^n_P\) of the first term is a
    convex combination of the old temperatures. The second term is non-positive because \(\mathbf A_n^{-1}\ge0\),
    \(\mathbf B\ge0\) and \(\mathbf q^n\ge0\). \(\square\)

No numerical overshoot can therefore appear in a reconstruction's forward model, whatever the time step: a
property that explicit and Crank–Nicolson schemes do not have.

## 7.5 Solving the linear systems

Each step requires the solution of an SPD system. For up to about 4000 cells BoilFlux uses a sparse direct
factorisation. Above that it uses the **conjugate-gradient method** with Jacobi (diagonal) preconditioning,
warm-started from the previous time level. CG converges at a rate set by the condition number
\(\kappa\) of the preconditioned matrix:
\(\lVert\mathbf e_k\rVert_{\mathbf A}\le 2\big((\sqrt\kappa-1)/(\sqrt\kappa+1)\big)^k\lVert\mathbf e_0\rVert_{\mathbf A}\) ([Appendix D](appendix-d-linear-algebra.md)).

A bound on \(\kappa\) follows from Gershgorin's theorem. Scaling \(\mathbf A_n\) by \(\mathbf C_n^{-1}\) (equivalent to Jacobi
preconditioning up to the loss term), the eigenvalues lie in \([1,\,1 + 2\max_P K_{PP}/C_{PP}]\), and for constant
properties \(K_{PP}/C_{PP} \le 2\alpha\Delta t\sum_d\Delta_d^{-2}\). Hence

\[
\kappa \le 1 + 4\,\alpha\Delta t\left(\frac{1}{\Delta x^2}+\frac{1}{\Delta y^2}+\frac{1}{\Delta z^2}\right).
\tag{7.6}
\]

For the BF-Bench grid (2.5 mm × 2.5 mm × 1 mm cells, \(\Delta t = 10\) ms, AA6082) this gives \(\kappa\le 4.8\): a handful
of CG iterations per step reach the tolerance \(10^{-11}\). The time step regularises the system, which is another
benefit of the implicit method.

## 7.6 Verification

A solver is **verified** when it is shown to solve its equations correctly; it is **validated** when the equations are
shown to describe reality. BoilFlux verifies the direct solver in three ways.

1. **Exact solutions.** Figure 2.1 compares the solver with the slab solution (2.14).
2. **Manufactured solutions.** Choose any smooth field, for example
   \(T^\star = T_0 + A\cos(\pi x/L_x)\cos(\pi z/(2L_z))\,e^{-t/\tau}\), insert it into (1.12) to compute the source and flux
   that make it exact, and check that the error falls as \(\Delta^2 + \Delta t\) when the grid and step are refined.
3. **Discrete identities.** The energy balance (7.5) and the symmetry of \(\mathbf A_n\) are tested to round-off, and the
   adjoint–sensitivity transpose identity of Chapter 10 to \(10^{-8}\).

**The inverse crime.** Verification of the *inverse* solver needs synthetic data. If those were generated by the
same discretisation used for inversion, model errors would vanish and the inversion would look better than it can
on real data. BoilFlux generates data on a grid refined by 2 in every direction and in time (Section 5.5).

## Exercises

!!! exercise "7.1"
    Derive the three amplification factors (7.2) and verify the stability limits of Result 7.1.

!!! exercise "7.2"
    Show that for a single cell with linear loss \(h\) the lagged step (7.4) is exactly the backward Euler step,
    and that for radiation it differs from it by \(O(\Delta t^2)\) per step.

!!! exercise "7.3"
    Evaluate (7.6) for nickel at 800 °C on the same grid. How many CG iterations does the bound suggest for a
    tolerance of \(10^{-11}\)?
