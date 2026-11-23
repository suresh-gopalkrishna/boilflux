# Moving-plate formulation

In moving-plate quench tests (run-out tables, plate and strip cooling, rigs that drive the specimen
past a stationary nozzle) the jet and the IR camera are fixed in the laboratory and the plate
travels with velocity \(v\) along its length.

## Two possible frames

| | Lab (Eulerian) frame | Plate (Lagrangian) frame (**QuenchFlux**) |
|---|---|---|
| Governing equation | \(\rho c_p (\partial_t T + v\,\partial_x T) = \nabla\cdot(k\nabla T)\) | \(\rho c_p\,\partial_t T = \nabla\cdot(k\nabla T)\) |
| Jet position | fixed | moves: \(x_s(t) = X_\text{jet} - X_0 - v t\) |
| Domain | a fixed window of a plate that flows through it, needing inflow/outflow conditions | the whole physical plate, with real edges |
| IR data | fixed pixels | resampled onto plate coordinates |
| Numerical issues | advection-dominated for large Péclet numbers, upwinding, artificial inflow | none beyond the stationary case |

QuenchFlux works in the **plate frame**. The physics is identical, but the discretisation stays
symmetric positive definite and the adjoint stays exact. Finite plates, plates entering and leaving
the camera view, and start/stop motions are all handled without special cases.

## Kinematics

The lab position of a plate point is

\[
X_\text{lab}(x, t) = X_0 + v\,t + x .
\]

[`PlateMotion`](../user-guide/api.md#kinematics) implements this mapping, the stagnation-point
trajectory, and the camera mask

\[
m(x, t) = \mathbb 1\big[X_\text{lab}(x, t) \in [X_\text{view,min}, X_\text{view,max}]\big].
\]

## From camera frames to plate data

`lab_frames_to_plate` resamples each IR frame (lab pixels \(X_\text{px}, Y_\text{px}\)) onto the
plate face grid by bilinear interpolation at \(X_\text{lab}(x_i, t_n)\). It returns the data and the
validity mask. Cells outside the field of view get \(m = 0\) and do not enter \(J\). Their flux is
still estimated from neighbouring information through conduction, with correspondingly larger
uncertainty.

!!! tip "Time-varying speed"
    For acceleration or stop-and-go motion, replace `PlateMotion.offset(t)` with the measured plate
    position (for example from an encoder). Everything downstream uses only this function.

## Relation to the 2D model of the 2026 paper

Setting `ny = 1`, `plate_velocity = 0` and putting the stagnation point at \(x=0\) gives the 2D
length × thickness model of Gopalkrishna et al. (2026). This is benchmark case
`S2D-AA6082-stationary`. The 3D model drops the assumption of negligible transverse gradients, which
matters near the stagnation zone and for moving plates, where the wetted region is not axisymmetric.
