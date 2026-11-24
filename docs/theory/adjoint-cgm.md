# Adjoint conjugate-gradient method

## Discretisation

Cell-centred finite volumes on a uniform \(n_x\times n_y\times n_z\) grid, with harmonic-mean face
conductances, and backward Euler in time. Properties are evaluated at the previous time level,
so each step is a single sparse symmetric positive-definite solve:

\[
A_n\,\mathbf T^n = C_n\,\mathbf T^{n-1} - B\,\mathbf q^n + \mathbf c_n,
\qquad A_n = C_n + K(\mathbf T^{n-1}) + L_n ,
\]

with \(C_n = \operatorname{diag}(\rho c_p V/\Delta t)\), \(K\) the conduction matrix, \(B\) the face-area
operator on \(\Gamma_j\), and \(L_n, \mathbf c_n\) the linearised IR-face loss. The observation
operator \(H\) extracts the IR-face cell temperatures.

## Discrete adjoint

With the coefficients frozen on the current trajectory, the Lagrangian of
\(J_d = \tfrac12\sum_n \lVert M_n(H\mathbf T^n - \mathbf Y^n)\rVert^2\) gives the adjoint recursion,
solved **backwards** in time:

\[
\boldsymbol\mu^{N_t+1} = 0,\qquad
A_n\,\boldsymbol\mu^n = H^\top M_n\,\mathbf r^n + C_{n+1}\,\boldsymbol\mu^{n+1},
\qquad n = N_t,\dots,1,
\]

and the gradient with respect to the optimisation variable \(\mathbf p = \mathbf q/Q\), \(Q = 10^6\) W/m²,

\[
\frac{\partial J}{\partial \mathbf p^n} = -Q\,B^\top\boldsymbol\mu^n + \gamma\,\mathbf p^n .
\]

This is the discrete counterpart of Eqs. (27)–(28) of the reference paper. Because the coefficients are frozen on
the current trajectory, it is the exact gradient of the **quasi-linearised** problem. The test suite checks three
things: agreement with finite differences to a relative \(10^{-6}\) for constant properties; that the adjoint and
sensitivity operators are exact transposes (to \(10^{-8}\), also for strongly temperature-dependent nickel); and
agreement with finite differences within 5 % for AA6082, where the frozen properties make the gradient slightly
inexact. The line search re-evaluates the true objective, so this affects only the search direction. Derivation:
[book, Chapter 10](../book/10-adjoint.md).

## Sensitivity problem and step size

For a search direction \(\mathbf d\), the sensitivity \(\boldsymbol\theta = \partial\mathbf T/\partial\mathbf q\,[\mathbf d]\)
solves

\[
A_n\,\boldsymbol\theta^n = C_n\,\boldsymbol\theta^{n-1} - B\,\mathbf d^n ,\qquad \boldsymbol\theta^0 = 0,
\]

and the step that minimises \(J\) along \(\mathbf d\) is

\[
\alpha = -\frac{\sum_n \langle \mathbf r^n, M_n H\boldsymbol\theta^n\rangle + \gamma\langle \mathbf q, \mathbf d\rangle}
{\sum_n \lVert M_n H\boldsymbol\theta^n\rVert^2 + \gamma\lVert\mathbf d\rVert^2}.
\]

## Algorithm

1. Solve the direct problem for \(\mathbf q_k\) and form the residuals \(\mathbf r^n\).
2. Stop if \(\sum \mathbf r^2 \le \tau^2 N_\text{obs}\,\sigma^2\) (discrepancy principle).
3. Solve the adjoint problem and form the gradient \(\mathbf g_k\).
4. Direction \(\mathbf d_k = -\mathbf g_k + \beta_k \mathbf d_{k-1}\), with the Polak–Ribière
   \(\beta_k = \max\!\big(0, \langle \mathbf g_k, \mathbf g_k-\mathbf g_{k-1}\rangle/\lVert\mathbf g_{k-1}\rVert^2\big)\).
5. Solve the sensitivity problem, compute \(\alpha_k\), and set \(\mathbf q_{k+1} = \mathbf q_k + \alpha_k\mathbf d_k\).
   Because properties depend on \(T\), the true objective is re-evaluated, and \(\alpha_k\) is halved if it increased.

The noise level \(\sigma\) is taken from the dataset, or estimated robustly from the second
temporal differences of the IR sequence (`estimate_noise_sigma`). The default \(\tau = 1.02\).

## End-time effect

Since \(\boldsymbol\mu^{N_t+1}=0\), the gradient vanishes at the final time and the recovered flux is
pulled towards the initial guess over the last diffusion time. BoilFlux reports
`result.reliable`, which excludes the last \(\lceil 0.5\,L_z^2/(\alpha(T_0)\,\Delta t)\rceil\) frames. Benchmark scores are
computed on this window.

## Cost

Each iteration needs three linear solves per time step: adjoint, sensitivity, and the forward solve at the
trial step, which is reused as the direct solve of the next iteration. Each halving of the step adds one forward
solve. Systems up to a few thousand cells use a sparse direct solver. Larger ones use
Jacobi-preconditioned conjugate gradients warm-started from the previous step. For the
48 × 24 × 5 benchmark grid with 200 frames, one CGM iteration takes about one second on a laptop CPU.
