# 11. The conjugate-gradient method in BoilFlux

<p class="bf-lead">This chapter assembles the pieces of Part IV into the algorithm that
<code>boilflux.AdjointCGM</code> runs: conjugate-gradient directions from the adjoint gradient, an exact step
length from the sensitivity problem, a safeguard for the nonlinearity, a robust estimate of the noise level, and
the discrepancy principle as the stopping rule.</p>

## 11.1 Conjugate gradients for a quadratic

Minimise \(\phi(\mathbf p) = \tfrac12\mathbf p^\top\mathbf A\mathbf p - \mathbf b^\top\mathbf p\) with \(\mathbf A\) symmetric positive definite. Its gradient is
\(\mathbf g = \mathbf A\mathbf p - \mathbf b\). Steepest descent zig-zags when \(\mathbf A\) is ill-conditioned. Conjugate gradients instead choose
search directions that are **\(\mathbf A\)-conjugate**, \(\mathbf d_i^\top\mathbf A\mathbf d_j = 0\) for \(i\neq j\).

!!! theorem "Result 11.1 (exact line search and the CG recurrence)"
    1. Along a direction \(\mathbf d\), \(\phi(\mathbf p + \alpha\mathbf d)\) is a parabola in \(\alpha\) with minimum at

        \[
        \alpha = -\frac{\mathbf g^\top\mathbf d}{\mathbf d^\top\mathbf A\mathbf d} .
        \tag{11.1}
        \]

    2. With \(\mathbf d_0 = -\mathbf g_0\) and \(\mathbf d_k = -\mathbf g_k + \beta_k\mathbf d_{k-1}\), the choice
       \(\beta_k = \lVert\mathbf g_k\rVert^2/\lVert\mathbf g_{k-1}\rVert^2\) makes all directions mutually conjugate and all gradients mutually
       orthogonal.
    3. Consequently \(\mathbf p_k\) minimises \(\phi\) over \(\mathbf p_0 + \operatorname{span}\{\mathbf d_0,\dots,\mathbf d_{k-1}\}\), and the exact minimiser is
       reached in at most \(n\) steps.

??? derivation "Sketch of proof"
    Item 1: \(\mathrm d\phi(\mathbf p+\alpha\mathbf d)/\mathrm d\alpha = \mathbf g^\top\mathbf d + \alpha\,\mathbf d^\top\mathbf A\mathbf d = 0\). For items 2–3, show by induction that \(\mathbf g_k\)
    is orthogonal to \(\mathbf d_0,\dots,\mathbf d_{k-1}\) (exact line searches along conjugate directions) and that the
    directions span the Krylov space \(\mathcal K_k = \operatorname{span}\{\mathbf g_0,\mathbf A\mathbf g_0,\dots,\mathbf A^{k-1}\mathbf g_0\}\). Requiring
    \(\mathbf d_k^\top\mathbf A\mathbf d_{k-1} = 0\) gives \(\beta_k = \mathbf g_k^\top\mathbf A\mathbf d_{k-1}/\mathbf d_{k-1}^\top\mathbf A\mathbf d_{k-1}\). Using
    \(\mathbf A\mathbf d_{k-1} = (\mathbf g_k - \mathbf g_{k-1})/\alpha_{k-1}\) and the orthogonality of gradients reduces it to the stated ratio.
    Conjugacy with the older directions follows because \(\mathbf A\mathbf d_j\in\mathcal K_{j+2}\perp\mathbf g_k\). See
    [Appendix D](appendix-d-linear-algebra.md). \(\square\)

For the linearised misfit, \(\mathbf A = \mathbf X^\top\mathbf X + \gamma\mathbf I\), and CG is exactly the regularising iteration of Section 9.4.

## 11.2 Nonlinear conjugate gradients

For a non-quadratic objective, the same recurrence is used with the current gradient, but the equivalent forms of
\(\beta_k\) are no longer equivalent. BoilFlux offers

\[
\beta_k^\text{FR} = \frac{\lVert\mathbf g_k\rVert^2}{\lVert\mathbf g_{k-1}\rVert^2},\qquad
\beta_k^\text{PR+} = \max\!\left(0,\ \frac{\mathbf g_k^\top(\mathbf g_k - \mathbf g_{k-1})}{\lVert\mathbf g_{k-1}\rVert^2}\right) .
\tag{11.2}
\]

The default is Polak–Ribière with the non-negativity cut (PR+). When successive gradients are nearly equal, which
happens when the iteration stalls, \(\beta^\text{PR}\approx0\) and the method restarts with steepest descent automatically.
If a direction is not a descent direction (\(\mathbf g_k^\top\mathbf d_k\ge0\)), BoilFlux also restarts with \(\mathbf d_k = -\mathbf g_k\).

## 11.3 The step length

The objective (9.11) depends on \(\mathbf p\) through the nonlinear forward model. Linearise the temperatures along the
direction \(\mathbf d\): \(\mathbf T^n(\mathbf p+\alpha\mathbf d)\approx\mathbf T^n(\mathbf p) + \alpha\boldsymbol\theta^n\), where \(\boldsymbol\theta\) solves the sensitivity problem (10.14) for the flux
direction \(Q\mathbf d\). The masked residuals become \(\mathbf r^n + \alpha\,\mathbf M_n\mathbf H\boldsymbol\theta^n\), and

\[
J(\alpha) \approx \tfrac12\sum_n\lVert\mathbf r^n + \alpha\mathbf M_n\mathbf H\boldsymbol\theta^n\rVert^2 + \tfrac\gamma2\lVert\mathbf p + \alpha\mathbf d\rVert^2
\]

is a parabola in \(\alpha\). Setting its derivative to zero:

!!! theorem "Result 11.2 (optimal step)"
    \[
    \alpha = -\,\frac{\sum_n\langle\mathbf r^n,\ \mathbf M_n\mathbf H\boldsymbol\theta^n\rangle + \gamma\,\langle\mathbf p, \mathbf d\rangle}
    {\sum_n\lVert\mathbf M_n\mathbf H\boldsymbol\theta^n\rVert^2 + \gamma\,\lVert\mathbf d\rVert^2} .
    \tag{11.3}
    \]

    For constant properties this is the exact minimiser of \(J\) along \(\mathbf d\), the same as (11.1).

**Nonlinear safeguard.** With temperature-dependent properties, (11.3) minimises a model of \(J\). BoilFlux therefore
runs the full nonlinear forward model at \(\mathbf p + \alpha\mathbf d\), and halves \(\alpha\) (up to eight times) until the true objective has
not increased. The accepted forward solution is reused as the direct solve of the next iteration.

## 11.4 Estimating the noise level

The discrepancy principle needs \(\sigma\). When the dataset does not provide it, BoilFlux estimates it from the IR
sequence itself. Model each pixel's record as \(Y^n = T(t_n) + e^n\) with white noise \(e^n\). The second difference

\[
d^n = Y^{n+1} - 2Y^n + Y^{n-1} = \underbrace{T(t_{n+1}) - 2T(t_n) + T(t_{n-1})}_{\approx\,\Delta t^2\,\ddot T} + \big(e^{n+1} - 2e^n + e^{n-1}\big)
\tag{11.4}
\]

removes the smooth signal to \(O(\Delta t^2)\), and the noise part has variance \((1 + 4 + 1)\sigma^2 = 6\sigma^2\). To be robust against
outliers and against the remaining signal in fast transients, BoilFlux uses the median absolute deviation (MAD)
instead of the standard deviation. For a normal distribution, \(\operatorname{MAD} = \Phi^{-1}(3/4)\,\mathrm{sd}\approx0.6745\,\mathrm{sd}\), so

\[
\hat\sigma = \frac{1.4826\ \operatorname{MAD}(d)}{\sqrt6} .
\tag{11.5}
\]

Rapid physical transients (the rewetting surge) leave some signal in \(d^n\) and inflate the estimate slightly, which
makes the discrepancy stop slightly earlier and the result slightly smoother. When the camera noise is known from
calibration, pass it explicitly.

## 11.5 The algorithm

!!! definition "Algorithm 11.1 (`AdjointCGM.solve`)"
    **Input:** IR data \(\mathbf Y^n\), mask \(\mathbf M_n\), material, initial temperature, noise level \(\sigma\) (or estimate it by
    (11.5)), \(\tau = 1.02\), \(\gamma\ge0\).

    1. \(\mathbf p\leftarrow\mathbf p_0\) (zero or a warm start). Solve the forward problem; form \(\mathbf r^n\) and \(S = \sum_n\lVert\mathbf r^n\rVert^2\).
    2. **For** \(k = 0, 1, \dots\):
        1. **Stop** if \(S\le\tau^2N_\text{obs}\sigma^2\) (discrepancy principle), with \(N_\text{obs}\) the number of observed values.
        2. Solve the adjoint problem (10.12) backwards; gradient \(\mathbf g_k = -Q\mathbf B^\top\boldsymbol\mu + \gamma\mathbf p\).
        3. Direction \(\mathbf d_k = -\mathbf g_k + \beta_k\mathbf d_{k-1}\) with \(\beta^\text{PR+}\) (11.2); restart if not a descent direction.
        4. Solve the sensitivity problem (10.14) for \(Q\mathbf d_k\); step \(\alpha_k\) from (11.3).
        5. Forward solve at \(\mathbf p + \alpha_k\mathbf d_k\); halve \(\alpha_k\) until \(J\) does not increase.
        6. Accept: \(\mathbf p\leftarrow\mathbf p+\alpha_k\mathbf d_k\). Stop if the relative decrease of \(J\) is below \(10^{-7}\) (stagnation).
    3. **Output:** \(\mathbf q = Q\mathbf p\), the surface temperature (6.4), the model IR temperatures, the convergence history and the
       reliable window of Section 10.7.

## 11.6 Cost and memory

One iteration requires, per time step, one adjoint solve, one sensitivity solve and one forward solve (plus one per
halving of the step), all with SPD matrices of the same pattern. The forward trajectory must be stored for the
adjoint pass: \((N_t+1)\times N\) values, 9 MB for the BF-Bench grid (5760 cells, 200 steps). One iteration takes about one
second on a two-core laptop CPU for that grid; the BF-Bench cases stop by the discrepancy principle after 34–61
iterations (24–131 s, see [Benchmark](../benchmark.md)).

## 11.7 Practical guidance

- **Grid.** Resolve the plate thickness with at least 4–5 cells, and in-plane features no finer than about the
  thickness (Chapter 8): finer in-plane grids add unknowns that the data cannot determine.
- **Time step.** Use the camera frame interval. Shorter steps add exponentially weak singular components.
- **Record length.** Record at least \(0.5\,L^2/\alpha\) beyond the period of interest (end-time effect).
- **Noise level.** Prefer a calibrated \(\sigma\). An overestimated \(\sigma\) stops too early (over-smoothed flux); an
  underestimated one fits noise (oscillations near the end of the iteration).
- **Warm starts.** A good initial flux (for example from a previous experiment or a neural-operator prior,
  Chapter 16) reduces the iteration count but also biases the result where the data are uninformative. Report
  which was used.

## Exercises

!!! exercise "11.1"
    Show that for the linear problem with \(\gamma = 0\) and exact line search, Fletcher–Reeves CG applied to \(J\) produces the
    same iterates as CGLS in Section 9.4.

!!! exercise "11.2"
    Derive the variance \(6\sigma^2\) in (11.4) and show that the signal term contributes less than 10 % of \(\hat\sigma\) when
    \(\Delta t^2\lvert\ddot T\rvert < 0.25\,\sigma\).

!!! exercise "11.3"
    Estimate the memory needed to store the trajectory for a 200 × 100 × 8 grid and 1000 frames. Suggest a
    checkpointing scheme that trades memory for recomputation.
