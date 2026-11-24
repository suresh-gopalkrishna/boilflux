# Appendix B. Vector calculus and integral theorems

Let \(\Omega\subset\mathbb R^3\) be a bounded domain with piecewise smooth boundary \(\partial\Omega\) and outward unit normal
\(\hat{\mathbf n}\). All functions are assumed smooth enough for the integrals to exist.

## B.1 Divergence theorem

For a vector field \(\mathbf F\),

\[
\int_\Omega\nabla\cdot\mathbf F\,\mathrm dV = \oint_{\partial\Omega}\mathbf F\cdot\hat{\mathbf n}\,\mathrm dS .
\tag{B.1}
\]

It converts the boundary flux of the energy balance (1.1) into a volume integral (Section 1.4) and, read backwards,
turns each finite volume's balance into a sum of face fluxes (Section 6.2).

## B.2 Product rule and Green's identities

For a scalar \(u\) and a vector field \(\mathbf F\),

\[
\nabla\cdot(u\,\mathbf F) = \nabla u\cdot\mathbf F + u\,\nabla\cdot\mathbf F .
\tag{B.2}
\]

Applying (B.1) with \(\mathbf F = k\nabla w\) gives **Green's first identity** (with a variable coefficient \(k\)):

\[
\int_\Omega u\,\nabla\cdot(k\nabla w)\,\mathrm dV = \oint_{\partial\Omega} u\,k\,\partial_n w\,\mathrm dS - \int_\Omega k\,\nabla u\cdot\nabla w\,\mathrm dV .
\tag{B.3}
\]

Exchanging \(u\) and \(w\) and subtracting gives **Green's second identity**:

\[
\int_\Omega\big(u\,\nabla\cdot(k\nabla w) - w\,\nabla\cdot(k\nabla u)\big)\,\mathrm dV
= \oint_{\partial\Omega} k\,\big(u\,\partial_n w - w\,\partial_n u\big)\,\mathrm dS .
\tag{B.4}
\]

Identity (B.4) says that the operator \(\nabla\cdot(k\nabla\cdot)\) is **self-adjoint** up to boundary terms. It is the reason the
adjoint heat equation of Chapter 10 has the same spatial operator as the forward one, and it is the continuous form of
the symmetry \(\mathbf K^\top = \mathbf K\) of Lemma 6.1.

## B.3 Integration by parts in time

For functions \(u(t), w(t)\) on \([0,t_f]\),

\[
\int_0^{t_f} u\,\partial_t w\,\mathrm dt = \big[u\,w\big]_0^{t_f} - \int_0^{t_f} w\,\partial_t u\,\mathrm dt .
\tag{B.5}
\]

In the adjoint derivation, the boundary terms vanish because the sensitivity starts at zero (\(\theta(0) = 0\)) and the adjoint
ends at zero (\(\lambda(t_f) = 0\)). The minus sign in the remaining term turns the forward heat equation into the
backward one. Its discrete counterpart is the index shift (summation by parts) in the proof of Theorem 10.4.

## B.4 Moving frames and the material derivative

If a material point moves with velocity \(\mathbf v\), the rate of change of a field \(\phi(\mathbf X, t)\) following the point is the
**material derivative**

\[
\frac{\mathrm D\phi}{\mathrm Dt} = \frac{\partial\phi}{\partial t} + \mathbf v\cdot\nabla\phi .
\tag{B.6}
\]

For a rigid translation with constant \(\mathbf v = v\hat{\mathbf x}\), (B.6) is the chain rule (1.14). Reynolds' transport theorem,
\(\tfrac{\mathrm d}{\mathrm dt}\int_{V(t)}\phi\,\mathrm dV = \int_{V(t)}\big(\mathrm D\phi/\mathrm Dt + \phi\,\nabla\cdot\mathbf v\big)\,\mathrm dV\) for a material volume \(V(t)\),
reduces to it for rigid motion (\(\nabla\cdot\mathbf v = 0\)).

## B.5 Fourier transform

With \(\hat f(\kappa) = \int_{-\infty}^\infty f(x)\,e^{-\mathrm i\kappa x}\,\mathrm dx\):

| Property | Statement |
|---|---|
| derivative | \(\widehat{f'} = \mathrm i\kappa\hat f\), \(\widehat{f''} = -\kappa^2\hat f\) |
| convolution | \(\widehat{f * g} = \hat f\,\hat g\), with \((f*g)(x) = \int f(x-y)\,g(y)\,\mathrm dy\) |
| Gaussian | \(f = (4\pi a)^{-1/2}e^{-x^2/(4a)}\ \Longleftrightarrow\ \hat f = e^{-a\kappa^2}\) |
| Parseval | \(\int\lvert f\rvert^2\,\mathrm dx = \tfrac{1}{2\pi}\int\lvert\hat f\rvert^2\,\mathrm d\kappa\) |

These give the heat-kernel result (14.3) and underlie the attenuation analysis of Section 2.6.
