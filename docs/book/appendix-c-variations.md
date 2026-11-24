# Appendix C. Calculus of variations and adjoint operators

## C.1 Directional derivatives of functionals

Let \(J\) map functions (or vectors) \(q\) in a Hilbert space with inner product \(\langle\cdot,\cdot\rangle\) to real numbers. The **directional
(Gâteaux) derivative** of \(J\) at \(q\) in the direction \(\delta q\) is

\[
\delta J(q;\delta q) = \lim_{\varepsilon\to0}\frac{J(q + \varepsilon\,\delta q) - J(q)}{\varepsilon} .
\tag{C.1}
\]

If \(\delta J(q;\cdot)\) is linear and bounded, the Riesz representation theorem gives a unique element \(\nabla J(q)\), the
**gradient**, with \(\delta J(q;\delta q) = \langle\nabla J(q),\,\delta q\rangle\) for all \(\delta q\). The gradient depends on the inner product.
In Chapter 10 the inner product on fluxes is \(\int_0^{t_f}\!\int_{\Gamma_j} q_1q_2\,\mathrm dS\,\mathrm dt\); in the discrete setting it is the
Euclidean product of the flux values. This is why the discrete gradient (10.13) carries the face area in \(\mathbf B\), while the
continuous gradient (10.9) does not.

**Example.** For \(J(q) = \tfrac12\lVert\mathcal Xq - y\rVert^2\) with a bounded linear operator \(\mathcal X\),
\(\delta J = \langle\mathcal Xq - y,\,\mathcal X\delta q\rangle = \langle\mathcal X^*(\mathcal Xq - y),\,\delta q\rangle\), so \(\nabla J = \mathcal X^*(\mathcal Xq - y)\).

## C.2 Adjoint operators

The **adjoint** of a bounded linear operator \(\mathcal X\) between Hilbert spaces is the operator \(\mathcal X^*\) defined by

\[
\langle\mathcal X u,\,w\rangle = \langle u,\,\mathcal X^* w\rangle\quad\text{for all } u, w .
\tag{C.2}
\]

For matrices with the Euclidean inner product, \(\mathcal X^* = \mathbf X^\top\). For the sensitivity operator of the quench problem,
\(u = \delta q\) on \(\Gamma_j\times(0,t_f)\) and \(\mathcal X u = \theta|_{\Gamma_c}\). Section 10.4 computes \(\mathcal X^*\) explicitly: it maps a function \(w\) on
\(\Gamma_c\times(0,t_f)\) to \(-\lambda|_{\Gamma_j}\), where \(\lambda\) solves the backward heat equation (10.8) with \(m(T-Y)\) replaced by \(w\).
Theorem 10.4 is (C.2) for the discrete operators.

## C.3 Constrained optimisation and Lagrange multipliers

To minimise \(j(T, q)\) subject to \(R(T,q) = 0\), introduce the Lagrangian

\[
\mathcal L(T, q, \lambda) = j(T,q) - \langle\lambda,\,R(T,q)\rangle .
\tag{C.3}
\]

At a constrained minimum, all partial derivatives of \(\mathcal L\) vanish (the KKT conditions for equality constraints):

\[
\partial_\lambda\mathcal L = 0\ \ (\text{state equation}),\qquad
\partial_T\mathcal L = 0\ \ (\text{adjoint equation}),\qquad
\partial_q\mathcal L = 0\ \ (\text{optimality}).
\tag{C.4}
\]

The adjoint method of Section 10.2 uses the first two conditions to evaluate the gradient \(\partial_q\mathcal L\) at any \(q\), not only at
the optimum. A gradient method then drives it to zero. All-at-once methods, including PINNs as penalty methods
(Section 13.3), work on all three conditions or their penalised forms simultaneously.

## C.4 The Euler–Lagrange equation

For \(J(u) = \int_a^b F(x, u, u')\,\mathrm dx\) with fixed end values, \(\delta J = \int_a^b(F_u\,\delta u + F_{u'}\,\delta u')\,\mathrm dx\). Integrating the
second term by parts and using \(\delta u(a) = \delta u(b) = 0\),

\[
\delta J = \int_a^b\Big(F_u - \frac{\mathrm d}{\mathrm dx}F_{u'}\Big)\delta u\,\mathrm dx
\quad\Longrightarrow\quad
F_u - \frac{\mathrm d}{\mathrm dx}F_{u'} = 0 .
\tag{C.5}
\]

The same manipulation (vary, integrate by parts, collect the boundary terms) produces the adjoint equations and their
boundary and terminal conditions in Chapter 10.
