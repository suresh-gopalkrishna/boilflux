# 10. Sensitivity and adjoint problems

<p class="bf-lead">To minimise the misfit over a flux field with hundreds of thousands of values, we need its
gradient. Finite differences would need one forward solve per value; the adjoint method needs one backward
solve in total. This chapter derives the adjoint method in general, then the continuous and the discrete
adjoint of the quench problem, proves the identity that BoilFlux tests, and explains the end-time effect.</p>

## 10.1 The cost of a gradient

For the BF-Bench grid, the flux has \(n = n_x n_y N_t = 48\times24\times200 = 230\,400\) unknowns. A one-sided finite-difference
gradient needs \(n+1\) forward solves, each of which is 200 sparse linear solves: about \(4.6\times10^7\) solves per
gradient. The **adjoint method** obtains the exact gradient of the discrete objective with one forward and one
backward solve, independently of \(n\).

## 10.2 The adjoint method in general

Let the state \(\mathbf T\) depend on parameters \(\mathbf q\) through a constraint \(\mathbf R(\mathbf T,\mathbf q) = \mathbf 0\) that has a unique solution
\(\mathbf T(\mathbf q)\), and let the objective be \(J(\mathbf q) = j(\mathbf T(\mathbf q), \mathbf q)\). Differentiating the constraint,

\[
\frac{\partial\mathbf R}{\partial\mathbf T}\,\frac{\mathrm d\mathbf T}{\mathrm d\mathbf q} + \frac{\partial\mathbf R}{\partial\mathbf q} = \mathbf 0
\quad\Longrightarrow\quad
\frac{\mathrm d\mathbf T}{\mathrm d\mathbf q} = -\Big(\frac{\partial\mathbf R}{\partial\mathbf T}\Big)^{-1}\frac{\partial\mathbf R}{\partial\mathbf q} .
\tag{10.1}
\]

By the chain rule,

\[
\nabla J^\top = \frac{\partial j}{\partial\mathbf q} + \frac{\partial j}{\partial\mathbf T}\,\frac{\mathrm d\mathbf T}{\mathrm d\mathbf q}
= \frac{\partial j}{\partial\mathbf q} - \underbrace{\frac{\partial j}{\partial\mathbf T}\Big(\frac{\partial\mathbf R}{\partial\mathbf T}\Big)^{-1}}_{\boldsymbol\lambda^\top}\frac{\partial\mathbf R}{\partial\mathbf q} .
\tag{10.2}
\]

There are two ways to evaluate the product. **Forward (sensitivity) mode** computes \(\mathrm d\mathbf T/\mathrm d\mathbf q\) column by column:
one linear solve with \(\partial\mathbf R/\partial\mathbf T\) per parameter. **Reverse (adjoint) mode** computes the row vector \(\boldsymbol\lambda^\top\) first:

!!! theorem "Result 10.1 (adjoint method)"
    Solve the **adjoint equation**

    \[
    \Big(\frac{\partial\mathbf R}{\partial\mathbf T}\Big)^{\!\top}\boldsymbol\lambda = \Big(\frac{\partial j}{\partial\mathbf T}\Big)^{\!\top}
    \tag{10.3}
    \]

    once; then \(\nabla J = (\partial j/\partial\mathbf q)^\top - (\partial\mathbf R/\partial\mathbf q)^\top\boldsymbol\lambda\). Equivalently, \(\nabla J\) is the derivative with respect to \(\mathbf q\)
    of the **Lagrangian** \(\mathcal L(\mathbf T,\mathbf q,\boldsymbol\lambda) = j(\mathbf T,\mathbf q) - \boldsymbol\lambda^\top\mathbf R(\mathbf T,\mathbf q)\) at a point where \(\partial\mathcal L/\partial\mathbf T = \mathbf 0\).

The cost is one solve with the transposed Jacobian, whatever the number of parameters. The same idea, applied
to a neural network, is back-propagation (Section 12.3).

## 10.3 The continuous sensitivity problem

Let \(T\) solve (1.12) for the flux \(q\), and perturb the flux to \(q + \varepsilon\,\delta q\). The **sensitivity**
\(\theta = \partial T/\partial\varepsilon\,|_{\varepsilon=0}\) is the directional derivative of the temperature field.

**Full linearisation.** With the enthalpy form (1.10), \(\partial_t H(T) = \nabla\cdot(k(T)\nabla T)\). Differentiating, and using
\(H' = \rho c_p\) and \(\nabla(k\theta) = k\nabla\theta + \theta k'(T)\nabla T\):

\[
\partial_t\big(\rho c_p(T)\,\theta\big) = \Delta\big(k(T)\,\theta\big),\qquad
-\partial_n\big(k\theta\big) = \delta q\ \text{on }\Gamma_j,\qquad
\partial_n(k\theta) = 0\ \text{elsewhere},\qquad \theta(\cdot,0) = 0 .
\tag{10.4}
\]

(With a loss on \(\Gamma_c\), the condition there becomes \(-\partial_n(k\theta) = q'_\text{loss}\theta\).) This is the exact derivative.

**Frozen coefficients.** If \(\rho c_p\) and \(k\) are evaluated on the current trajectory and held fixed as functions of
\((\mathbf x,t)\), the derivative of the property dependence is dropped:

\[
\rho c_p\,\partial_t\theta = \nabla\cdot(k\nabla\theta),\qquad -k\,\partial_n\theta = \delta q\ \text{on }\Gamma_j,\qquad\dots
\tag{10.5}
\]

The two agree for constant properties. BoilFlux uses the discrete counterpart of (10.5) (Section 10.5), and
Section 10.6 quantifies the difference.

## 10.4 The continuous adjoint problem

Take the objective

\[
J(q) = \frac12\int_0^{t_f}\!\!\int_{\Gamma_c} m\,(T - Y)^2\,\mathrm dS\,\mathrm dt + \frac{\gamma}{2}\int_0^{t_f}\!\!\int_{\Gamma_j}q^2\,\mathrm dS\,\mathrm dt ,
\tag{10.6}
\]

with the observation mask \(m\in\{0,1\}\). Its directional derivative is

\[
\delta J = \int_0^{t_f}\!\!\int_{\Gamma_c} m\,(T-Y)\,\theta\,\mathrm dS\,\mathrm dt + \gamma\int_0^{t_f}\!\!\int_{\Gamma_j} q\,\delta q\,\mathrm dS\,\mathrm dt .
\tag{10.7}
\]

The first term depends on \(\delta q\) only implicitly, through \(\theta\). The adjoint turns it into an explicit integral over
\(\Gamma_j\). Multiply the frozen sensitivity equation (10.5) by a function \(\lambda(\mathbf x,t)\), integrate over space and time,
integrate by parts once in time and twice in space (Green's second identity, B.4):

\[
0 = \Big[\int_\Omega\lambda\rho c_p\theta\,\mathrm dV\Big]_0^{t_f}
+ \int_0^{t_f}\!\!\int_\Omega\theta\Big[-\partial_t(\rho c_p\lambda) - \nabla\cdot(k\nabla\lambda)\Big]\mathrm dV\,\mathrm dt
+ \int_0^{t_f}\!\!\oint_{\partial\Omega}\big(\theta\,k\,\partial_n\lambda - \lambda\,k\,\partial_n\theta\big)\,\mathrm dS\,\mathrm dt .
\]

Now choose \(\lambda\) to make every term that contains \(\theta\) in an unknown way vanish or become the misfit term:

!!! theorem "Result 10.2 (continuous adjoint problem, frozen coefficients)"
    Let \(\lambda\) solve, **backwards in time**,

    \[
    \begin{aligned}
    -\partial_t(\rho c_p\lambda) &= \nabla\cdot(k\nabla\lambda) && \text{in }\Omega\times[0,t_f),\\
    k\,\partial_n\lambda &= m\,(T - Y) && \text{on }\Gamma_c,\\
    k\,\partial_n\lambda &= 0 && \text{on }\Gamma_j\cup\Gamma_e,\\
    \lambda(\cdot,t_f) &= 0 && \text{in }\Omega .
    \end{aligned}
    \tag{10.8}
    \]

    Then

    \[
    \nabla J(q) = -\lambda\big|_{\Gamma_j} + \gamma\,q .
    \tag{10.9}
    \]

??? derivation "Verification"
    With (10.8), the time boundary term vanishes (\(\theta=0\) at \(t=0\), \(\lambda=0\) at \(t_f\)) and so does the volume term.
    In the boundary term, \(k\,\partial_n\theta = -\delta q\) on \(\Gamma_j\) and 0 elsewhere, while \(k\,\partial_n\lambda = m(T-Y)\) on \(\Gamma_c\) and 0
    elsewhere. What remains is

    \[
    0 = \int_0^{t_f}\!\!\int_{\Gamma_c} m(T-Y)\,\theta\,\mathrm dS\,\mathrm dt + \int_0^{t_f}\!\!\int_{\Gamma_j}\lambda\,\delta q\,\mathrm dS\,\mathrm dt ,
    \]

    so the misfit term of (10.7) equals \(-\int\!\!\int_{\Gamma_j}\lambda\,\delta q\). Then
    \(\delta J = \int\!\!\int_{\Gamma_j}(-\lambda + \gamma q)\,\delta q\), which identifies the gradient. \(\square\)

**Interpretation.** The adjoint is a heat equation run backwards from \(t_f\). The misfit on the camera face acts as a
heat flux injected into the plate; the adjoint "heat" diffuses back to the quenched face, and its value there is the
gradient. If the model is too hot where the camera looks (\(T>Y\)), \(\lambda>0\) and the gradient \(-\lambda<0\): increase the
extracted flux, which is the physically correct correction.

**Full linearisation.** Repeating the argument with (10.4) instead of (10.5) gives the non-divergence form
\(-\rho c_p\,\partial_t\lambda = k\,\Delta\lambda\) with \(\partial_n\lambda = m(T-Y)/k\) on \(\Gamma_c\) and the same gradient (10.9). The two adjoints differ
by the terms \(\lambda\,\partial_t(\rho c_p)\) and \(\nabla k\cdot\nabla\lambda\), which are proportional to
\((\mathrm d\rho c_p/\mathrm dT)\,\partial_t T\) and \((\mathrm dk/\mathrm dT)\,\nabla T\): small when the properties vary slowly with temperature.

## 10.5 The discrete adjoint of BoilFlux

The continuous adjoint (10.8) would have to be discretised separately, and its discretisation would not give the
exact gradient of the discrete objective (**optimise-then-discretise**). BoilFlux instead differentiates the discrete
equations themselves (**discretise-then-optimise**), so the gradient is exact for the discrete problem and
consistent with the line search of Chapter 11.

The forward model is the sequence of steps (7.4), with coefficients frozen on the current trajectory:

\[
\mathbf A_n\mathbf T^n - \mathbf C_n\mathbf T^{n-1} + \mathbf B\mathbf q^n - \mathbf c_n = \mathbf 0,\qquad n = 1,\dots,N_t .
\tag{10.10}
\]

The observation operator \(\mathbf H\) extracts the last cell layer, \(\mathbf M_n\) is the 0/1 mask of observed pixels at step \(n\), and the
masked residual is \(\mathbf r^n = \mathbf M_n(\mathbf H\mathbf T^n - \mathbf Y^n)\). The objective is
\(J = \tfrac12\sum_n\lVert\mathbf r^n\rVert^2 + \tfrac\gamma2\sum_n\lVert\mathbf q^n\rVert^2\). Form the Lagrangian with one multiplier vector \(\boldsymbol\mu^n\) per step:

\[
\mathcal L = J - \sum_{n=1}^{N_t}\boldsymbol\mu^{n\top}\big(\mathbf A_n\mathbf T^n - \mathbf C_n\mathbf T^{n-1} + \mathbf B\mathbf q^n - \mathbf c_n\big) .
\tag{10.11}
\]

The state \(\mathbf T^n\) appears in step \(n\) (through \(\mathbf A_n\)) and in step \(n+1\) (through \(\mathbf C_{n+1}\)). Setting
\(\partial\mathcal L/\partial\mathbf T^n = \mathbf 0\), with \(\mathbf A_n\) and \(\mathbf C_n\) symmetric and \(\mathbf M_n^2 = \mathbf M_n\):

\[
\mathbf H^\top\mathbf r^n - \mathbf A_n\boldsymbol\mu^n + \mathbf C_{n+1}\boldsymbol\mu^{n+1} = \mathbf 0 .
\]

!!! theorem "Result 10.3 (discrete adjoint and gradient)"
    Set \(\boldsymbol\mu^{N_t+1} = \mathbf 0\) and solve backwards for \(n = N_t,\dots,1\)

    \[
    \mathbf A_n\,\boldsymbol\mu^n = \mathbf H^\top\mathbf r^n + \mathbf C_{n+1}\,\boldsymbol\mu^{n+1} .
    \tag{10.12}
    \]

    Then the gradient of \(J\) with respect to the flux of step \(n\) is

    \[
    \frac{\partial J}{\partial\mathbf q^n} = -\mathbf B^\top\boldsymbol\mu^n + \gamma\,\mathbf q^n .
    \tag{10.13}
    \]

    Each adjoint step costs one solve with the same SPD matrix \(\mathbf A_n\) as the forward step.

In BoilFlux the optimisation variable is \(\mathbf p = \mathbf q/Q\) with \(Q = 10^6\) W/m², so by the chain rule
\(\partial J/\partial\mathbf p^n = -Q\,\mathbf B^\top\boldsymbol\mu^n + \gamma\,\mathbf p^n\), with \(\gamma\) acting on \(\mathbf p\) as in (9.11).

!!! note "In BoilFlux"
    `ForwardModel.adjoint` implements (10.12) and returns \(-\mathbf B^\top\boldsymbol\mu^n\);
    `AdjointCGM._gradient` adds the scaling and the Tikhonov term.

The discrete sensitivity problem follows by differentiating (10.10) with respect to \(\mathbf q\) in a direction \(\mathbf d\):

\[
\mathbf A_n\boldsymbol\theta^n = \mathbf C_n\boldsymbol\theta^{n-1} - \mathbf B\,\mathbf d^n,\qquad \boldsymbol\theta^0 = \mathbf 0 .
\tag{10.14}
\]

!!! theorem "Theorem 10.4 (adjoint–sensitivity duality)"
    For any weights \(\mathbf w^n\) on the observed face and any flux direction \(\mathbf d^n\), let \(\boldsymbol\theta\) solve (10.14) and let \(\boldsymbol\mu\)
    solve (10.12) with \(\mathbf H^\top\mathbf w^n\) in place of \(\mathbf H^\top\mathbf r^n\). Then

    \[
    \sum_{n=1}^{N_t}\mathbf w^{n\top}\mathbf H\boldsymbol\theta^n = \sum_{n=1}^{N_t}\big(-\mathbf B^\top\boldsymbol\mu^n\big)^{\!\top}\mathbf d^n .
    \tag{10.15}
    \]

    In words: the adjoint map is exactly the transpose of the sensitivity map.

??? derivation "Proof"
    Multiply (10.14) by \(\boldsymbol\mu^{n\top}\) and sum over \(n\):
    \(\sum_n\boldsymbol\mu^{n\top}\mathbf A_n\boldsymbol\theta^n - \sum_n\boldsymbol\mu^{n\top}\mathbf C_n\boldsymbol\theta^{n-1} = -\sum_n\boldsymbol\mu^{n\top}\mathbf B\mathbf d^n\).
    Shift the index of the second sum (\(n\to n+1\)), use \(\boldsymbol\theta^0 = \mathbf 0\), \(\boldsymbol\mu^{N_t+1} = \mathbf 0\) and the symmetry of \(\mathbf A_n, \mathbf C_n\):
    the left side becomes \(\sum_n\boldsymbol\theta^{n\top}(\mathbf A_n\boldsymbol\mu^n - \mathbf C_{n+1}\boldsymbol\mu^{n+1}) = \sum_n\boldsymbol\theta^{n\top}\mathbf H^\top\mathbf w^n\) by (10.12). \(\square\)

## 10.6 What exactly is exact

Because the coefficients in (10.10) are frozen, Result 10.3 is the exact gradient of the **quasi-linearised** objective:
the one in which \(\mathbf A_n\) and \(\mathbf C_n\) do not respond to a change in \(\mathbf q\). The full discrete gradient has extra terms from
\(\partial\mathbf A_n/\partial\mathbf T^{n-1}\) and \(\partial\mathbf C_n/\partial\mathbf T^{n-1}\), the discrete counterparts of the terms identified at the end of Section 10.4.
The BoilFlux test suite measures the difference:

| Test | Material | Result |
|---|---|---|
| gradient vs central finite differences | constant properties | agreement to a relative \(10^{-6}\): exact |
| transpose identity (10.15) | nickel (strongly temperature-dependent) | agreement to \(10^{-8}\): exact by construction |
| gradient vs central finite differences | AA6082 (weakly temperature-dependent) | agreement within 5 % |

The inexactness affects only the search direction, not convergence: the CGM of Chapter 11 re-evaluates the true
nonlinear objective after every step and halves the step if it increased. Using the full linearisation (10.4) would
make the gradient exact for temperature-dependent properties too, at the cost of assembling the property
derivatives; it is on the [roadmap](../roadmap.md).

## 10.7 The end-time effect

The terminal condition \(\boldsymbol\mu^{N_t+1} = \mathbf 0\) (continuous: \(\lambda(t_f) = 0\)) has a physical consequence. Adjoint heat injected at
the camera face near \(t_f\) has no time to diffuse back to the quenched face, so the gradient (10.13) is nearly zero for
the last steps. The flux there is hardly updated from its initial guess. Physically: a flux change in the last
moments of the record has not yet reached the camera, so the data contain no information about it.

The affected window is about one through-thickness diffusion time. BoilFlux reports as reliable all frames except
the last \(\lceil 0.5\,L_z^2/(\alpha(T_0)\,\Delta t)\rceil\) and computes benchmark scores on that window. Experiments should
record correspondingly beyond the period of interest.

## Exercises

!!! exercise "10.1"
    Derive the full-linearisation adjoint stated at the end of Section 10.4 and show that its difference from (10.8)
    vanishes for constant properties.

!!! exercise "10.2"
    Add the linearised loss \(-k\,\partial_n\theta = q'_\text{loss}\theta\) on \(\Gamma_c\) to (10.5) and show that the adjoint condition there
    becomes \(k\,\partial_n\lambda + q'_\text{loss}\lambda = m(T-Y)\).

!!! exercise "10.3"
    Implement the Taylor test: for a random direction \(\mathbf d\), check that
    \(\lvert J(\mathbf q+\varepsilon\mathbf d) - J(\mathbf q) - \varepsilon\nabla J^\top\mathbf d\rvert = O(\varepsilon^2)\) with constant properties, and \(O(\varepsilon)\) with nickel.
