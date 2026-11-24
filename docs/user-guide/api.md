# Python API

All public objects are importable from the top-level package: `import boilflux as bf`.

## Geometry and materials

`bf.PlateGrid(Lx, Ly, Lz, nx, ny, nz)`
:   Uniform cell-centred grid. Properties: `dx, dy, dz, x, y, z, shape, face_shape, n_cells`.
    `refined(f)` returns a refined grid and `coarsen_face(field, f)` block-averages face fields.

`bf.Material(name, T_table, k_table, rho_table, cp_table)`
:   Methods `k(T)`, `rho(T)`, `cp(T)`, `rho_cp(T)`, `diffusivity(T)`. Constructors
    `Material.constant(...)` and `Material.from_csv(path)`. Built-ins via `bf.get_material(name)`.

## Kinematics

`bf.PlateMotion(velocity=0.0, x0=0.0)`
:   `offset(t)`, `lab_to_plate(X, t)`, `plate_to_lab(x, t)`, `stagnation_point(X_jet, t)`,
    `view_mask(grid, t, (Xmin, Xmax))`.

`bf.lab_frames_to_plate(frames, X_px, Y_px, t, grid, motion, y_offset=0.0) -> (Y, mask)`
:   Resample lab-frame IR frames onto the plate grid.

## Forward model

`bf.ForwardModel(grid, material, dt, nt, T0, top_loss=None, solver="auto")`
:   `run(q=...)` or `run(flux_fn=...)` returns a `Trajectory` (`T`, `q`, `t`, `top_face()`,
    `bottom_cell()`, `field(n)`). `sensitivity(traj, dq)` and `adjoint(traj, weighted_residual)`
    are the linearised operators used by the inverse solver.

`bf.SurfaceLoss(h=10.0, emissivity=0.9, T_amb=25.0)`
:   Convection and radiation on the IR face.

## Data

`bf.BoilDataset(grid, material, dt, T0, Y, mask=None, sigma=None, q_true=None, motion=..., meta={})`
:   `save(path)` and `BoilDataset.load(path)`; properties `nt`, `t`.

`bf.JetQuenchScenario(...)`
:   Synthetic single-jet scenario; `generate()` returns a `BoilDataset` with `q_true`.

`bf.SpechtBoilingCurve(q_max, T_dnb, q_leid, T_leid, T_f, Q_dot)`
:   Callable boiling curve; `flux_and_derivative(Ts)`. Defaults `bf.IZ_DEFAULT`, `bf.FFZ_DEFAULT`.

## Inverse solver

`bf.AdjointCGM(data, material=None, top_loss=None, options=CGMOptions())`
:   `solve(p0=None, callback=None) -> InverseResult`; `end_effect_frames()`.

`bf.CGMOptions(max_iter=200, gamma=0.0, tau=1.02, stop_on_discrepancy=True, rel_tol=1e-7, beta="polak-ribiere")`

`bf.InverseResult`
:   `q`, `T_top`, `T_surface`, `trajectory`, `history` (`misfit`, `rms`, `alpha`), `sigma`,
    `converged`, `message`, `runtime_s`, `reliable` (time slice), `save(path)`.

## Utilities

`bf.estimate_noise_sigma(Y, mask=None)`, `bf.flux_metrics(q_est, q_true)`,
`boilflux.plotting.save_all(data, result, outdir)`, `boilflux.benchmarks.run_case(name)`.
