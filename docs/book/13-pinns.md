# 13. Physics-informed neural networks

<p class="bf-lead">A physics-informed neural network (PINN) represents the temperature field by a network and
penalises the residual of the heat equation together with the data misfit. This chapter formulates PINNs for
the forward and the inverse quench problem, shows that they are a penalty form of the optimisation problem of
Part IV (and therefore inherit its ill-posedness), and discusses loss balancing, failure modes and when a PINN
is the right tool.</p>

## 13.1 Formulation

Represent the temperature by a network \(T_\theta(\mathbf x, t)\) (Chapter 12) that takes space and time coordinates as input.
Because the network is an explicit differentiable function, automatic differentiation gives \(\partial_t T_\theta\),
\(\nabla T_\theta\) and \(\nabla\cdot(k(T_\theta)\nabla T_\theta)\) exactly at any point, without a grid. The **PDE residual** of (1.5) is

\[
\mathcal R_\theta(\mathbf x,t) = \rho c_p(T_\theta)\,\partial_t T_\theta - \nabla\cdot\big(k(T_\theta)\nabla T_\theta\big) .
\tag{13.1}
\]

A PINN chooses \(\theta\) by minimising a weighted sum of mean-square penalties evaluated at random **collocation
points** in the interior, on the boundaries, at \(t=0\), and at the measurement points:

\[
\mathcal J(\theta) = w_r\,\overline{\mathcal R_\theta^2}^{\,\Omega\times(0,t_f)}
+ w_b\,\overline{\big(k\partial_n T_\theta + q_b\big)^2}^{\,\partial\Omega}
+ w_0\,\overline{(T_\theta - T_0)^2}^{\,t=0}
+ w_d\,\overline{(T_\theta - Y)^2}^{\,\text{data}} ,
\tag{13.2}
\]

where the overline denotes the mean over the corresponding points and \(q_b\) is the prescribed outward flux.

**Nondimensionalise first.** Temperatures of order \(10^2\) K, fluxes of order \(10^7\) W/m², lengths of order
\(10^{-3}\) m and times of order 1 s give residual terms that differ by many orders of magnitude, and the optimiser of
Section 12.4 fails. Use the variables of (2.1): \(\zeta = z/L\), \(\mathrm{Fo} = \alpha_0 t/L^2\), \(\vartheta = (T-T_0)/\Delta T\) and
\(q^* = qL/(k_0\Delta T)\), with reference values \(\alpha_0, k_0\). Then every term in (13.2) is of order one.

**Hard constraints.** Conditions can be built into the network instead of penalised. For example,
\(T_\theta = T_0 + \mathrm{Fo}\,N_\theta(\mathbf x,\mathrm{Fo})\) satisfies the initial condition exactly for any network \(N_\theta\), which removes one
term from (13.2) and one source of imbalance.

## 13.2 The inverse problem with a PINN

In the inverse quench problem the flux on \(\Gamma_j\) is unknown. There are two formulations:

1. **Flux network.** Add a second network \(q_\phi(x,y,t)\) and include the flux condition on \(\Gamma_j\) as a penalty, optimising
   \(\theta\) and \(\phi\) together.
2. **Flux as a derived quantity.** Impose no condition on \(\Gamma_j\) at all. The data on \(\Gamma_c\), the conditions on \(\Gamma_c\) and the
   edges, and the residual determine \(T_\theta\); the flux is then **computed** from it,
   \(q = k(T_\theta)\,\partial_z T_\theta\big|_{z=0}\).

The second form is closer to the physics: the IR face carries two conditions (the measured temperature and the
known loss), the quenched face none. It also shows plainly that the PINN must extrapolate the field across the
plate, from the data face to the face where the flux lives, which is the ill-posed step of Chapter 8.

## 13.3 PINNs as a penalty method: what they inherit

Discretise conceptually and write the PDE constraint as \(\mathbf R(\mathbf T,\mathbf q) = \mathbf 0\). The adjoint CGM of Part IV solves the
**reduced** problem: for every \(\mathbf q\) it satisfies the constraint exactly with the forward solver and minimises over
\(\mathbf q\) alone. A PINN solves the **all-at-once penalty** problem

\[
\min_{\mathbf T,\mathbf q}\ \tfrac12\lVert\mathbf H\mathbf T - \mathbf Y\rVert^2 + \tfrac{w}{2}\lVert\mathbf R(\mathbf T,\mathbf q)\rVert^2 ,
\tag{13.3}
\]

with \(\mathbf T\) and \(\mathbf q\) restricted to what the networks can represent. Two consequences:

- **As \(w\to\infty\)** the penalty enforces the constraint and (13.3) tends to the reduced least-squares problem, whose
  solution is the unregularised \(\mathbf q_\text{LS}\) of (8.7). A PINN that fits the PDE and the data perfectly reproduces the
  noise amplification of Chapter 8. **The physics constraint does not regularise the inverse problem**; it is the
  same constraint the adjoint method enforces exactly.
- **What does regularise a PINN** is implicit: the limited capacity of the network, the spectral bias of Section 12.5,
  and early stopping of the optimiser, the counterparts of the iterative regularisation of Section 9.4. They are
  harder to control than a discrepancy stop. An explicit penalty on the flux, for example
  \(\gamma\,\overline{(\partial_t q)^2}\), and stopping when the data misfit reaches the noise level (9.7) restore control.

A finite penalty weight \(w\) means the PDE is satisfied only approximately. The reconstructed flux then differs
from one that is exactly consistent with conduction in the plate, and the energy balance (1.18) need not hold. It
should be checked after training.

## 13.4 Balancing the loss terms

The terms of (13.2) have gradients of very different sizes, and the largest dominates training. A widely used
remedy (Wang, Teng and Perdikaris, 2021) adapts the weights during training so that each term's gradient has a
size comparable with the residual term's:

\[
\hat w_i = \frac{\max_\theta\lvert\nabla_\theta\mathcal J_r\rvert}{\overline{\lvert\nabla_\theta\mathcal J_i\rvert}},\qquad
w_i \leftarrow (1-\beta)\,w_i + \beta\,\hat w_i ,
\tag{13.4}
\]

with a moving-average factor \(\beta\approx0.1\). For time-dependent problems, **causal weighting** multiplies the residual at
time slice \(n\) by \(\exp\!\big(-\epsilon\sum_{j<n}\mathcal J_r^{(j)}\big)\), so the network does not try to fit late times before early times are
right, which mirrors the causal structure of the time-stepping solver.

## 13.5 Failure modes in quenching

| Difficulty | Cause | Remedy |
|---|---|---|
| sharp wetting fronts and the rewetting surge are smeared | spectral bias (Section 12.5) | Fourier features, domain decomposition, adaptive collocation near fronts |
| training stalls with a large residual | loss imbalance, unscaled variables | nondimensionalisation, adaptive weights (13.4), hard constraints |
| oscillating flux near the end of the record | the end-time effect (Section 10.7) appears here too: late flux is unobservable | discard the last \(0.5\,t_d\) as for the CGM |
| one costly training per experiment | the network represents one solution, not a solver | neural operators (Chapter 14) |
| result depends on the random seed | non-convex optimisation | ensembles of trainings; report the spread (Chapter 15) |

## 13.6 When to use a PINN

| | Adjoint CGM (BoilFlux) | PINN |
|---|---|---|
| PDE satisfied | exactly, by the forward solver | approximately, by penalty |
| regularisation | explicit and controlled (discrepancy stop, \(\gamma\)) | implicit (capacity, spectral bias, early stopping) unless added |
| gradient | exact adjoint, one backward solve | automatic differentiation through the network |
| geometry and grids | structured grid | mesh-free; irregular data points are natural |
| unknown properties or model parameters | requires extending the adjoint | add them to \(\theta\) |
| cost per reconstruction | seconds to minutes | minutes to hours of training |

PINNs are attractive when the data are scattered, when unknown model parameters (a contact resistance, an
unknown \(k(T)\)) must be identified together with the flux, or when the PDE is to be combined with other learned
components. For the standard IR quench problem on a structured plate, the adjoint method is faster and its
regularisation is better controlled. BoilFlux plans a PINN **baseline** on BF-Bench so that the two can be compared on
the same data and scores ([roadmap](../roadmap.md)).

## Exercises

!!! exercise "13.1"
    Write (13.1) in the dimensionless variables of Section 13.1 with \(k = k_0\kappa(\vartheta)\) and
    \(\rho c_p = (\rho c_p)_0\,c(\vartheta)\).

!!! exercise "13.2"
    For the linear 1D problem (8.6), show that minimising \(\tfrac12\lVert\mathbf X\mathbf q-\mathbf y\rVert^2\) over \(\mathbf q\) in the span of the first \(k\)
    right singular vectors is truncated SVD. What does this suggest about network capacity as a regulariser?

!!! exercise "13.3"
    Explain why the hard constraint \(T_\theta = T_0 + \mathrm{Fo}\,N_\theta\) does not restrict which solutions the network can
    represent, for any \(T_0\) that is constant in space.
