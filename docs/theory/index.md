# Theory overview

## Physical problem

A hot metal plate of length \(L_x\), width \(L_y\) and thickness \(L_z\), initially at \(T_0\), is
quenched by a water jet impinging on the face \(\Gamma_j\) (\(z = 0\)). An infrared camera records the
temperature of the opposite face \(\Gamma_c\) (\(z = L_z\)), which is coated for high emissivity.
The plate may move along \(x\) with velocity \(v\) relative to the jet and camera.

The jet produces several boiling regimes on \(\Gamma_j\) at the same time: nucleate and transition boiling
around the stagnation point, a rewetting front moving outwards, and film boiling further away. The
resulting heat flux \(q(x, y, t)\) spans one to two orders of magnitude within millimetres and
fractions of a second. It is the unknown of the problem.

## Direct problem

In the plate frame (see [moving-plate formulation](moving-plate.md)):

\[
\rho(T)\,c_p(T)\,\frac{\partial T}{\partial t} = \nabla\cdot\big(k(T)\,\nabla T\big)
\quad\text{in } \Omega\times(0, t_f],
\qquad T(\mathbf{x}, 0) = T_0,
\]

\[
-k(T)\,\nabla T\cdot\hat{\mathbf n} = q(x, y, t)\ \text{on } \Gamma_j,\qquad
-k(T)\,\nabla T\cdot\hat{\mathbf n} = q_\text{loss}(T)\ \text{on } \Gamma_c,\qquad
-k(T)\,\nabla T\cdot\hat{\mathbf n} = 0\ \text{on the edges}.
\]

A positive \(q\) extracts heat from the plate. The loss on the IR face is optional:

\[
q_\text{loss}(T) = h\,(T - T_\infty) + \varepsilon\sigma\,(T^4 - T_\infty^4).
\]

## Inverse problem

Given noisy IR temperatures \(Y(x, y, t_n)\) on \(\Gamma_c\), find \(q\) that minimises the Tikhonov functional

\[
J(q) = \frac12 \sum_{n=1}^{N_t}\sum_{i} m_{n,i}\,\big(T(\mathbf x_i, t_n; q) - Y_{n,i}\big)^2
      + \frac{\gamma}{2}\,\lVert q\rVert^2 ,
\]

where \(m_{n,i}\in\{0,1\}\) is the observation mask (camera field of view, bad pixels). This is the
3D, moving-plate generalisation of Eq. (25) of the reference paper. The problem is ill-posed:
conduction smooths high-frequency content of \(q\) exponentially before it reaches \(\Gamma_c\), so
regularisation is essential. BoilFlux uses **iterative regularisation**: conjugate-gradient iterations stopped by
Morozov's discrepancy principle, optionally combined with zeroth-order Tikhonov regularisation (\(\gamma>0\)).

Details: [Adjoint CGM](adjoint-cgm.md).

## Assumptions and limits

- Conduction only inside the solid. Latent heat of solid-state phase transformations is not modelled
  (add it through an effective \(c_p(T)\)).
- IR temperatures are assumed already calibrated. Emissivity errors pass directly into the result.
- The unknown is a free field \(q(x, y, t)\). It is not parameterised by a boiling-curve model, so
  boiling curves are *outputs*, not assumptions.
- The flux recovered near the final time is biased (end-time effect). Discard the last
  \(\approx 0.5\,L_z^2/\alpha\) seconds.
- Resolution: features smaller than about the plate thickness are attenuated through the plate. This is a physical
  limit of back-face measurement, not a numerical one.
