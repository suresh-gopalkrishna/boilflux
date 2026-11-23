# Changelog

All notable changes are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[Semantic Versioning](https://semver.org/).

## [0.1.0] - 2026-11-23

### Added
- Nonlinear 3D cell-centred finite-volume forward model of a quenched plate with
  temperature-dependent properties, implicit Euler time stepping and optional
  radiation/convection loss on the IR face. `ny = 1` reproduces the 2D model of
  Gopalkrishna et al., Applied Thermal Engineering 304 (2026) 132481.
- Moving-plate support in the plate (Lagrangian) frame: moving stagnation point,
  time-dependent camera field-of-view masks and lab-to-plate resampling of IR frames.
- Adjoint conjugate-gradient inverse solver (Polak-Ribiere, exact discrete adjoint and
  sensitivity of the frozen-coefficient problem, Tikhonov term, Morozov discrepancy
  stopping, nonlinear backtracking safeguard, end-time-effect window).
- Specht parabolic boiling-curve model and synthetic single-jet quench generator
  (impingement / free-flow zones, liquid coverage, refined-grid data to avoid inverse crime).
- Robust IR noise estimation, error metrics, plots, NPZ dataset format, CLI
  (`boilflux info | demo | bench | synth | invert`).
- BF-Bench tier-2 cases: `S2D-AA6082-stationary`, `M3D-AA6082`, `M3D-Ni`, `M3D-AA6082-window`.
- Documentation website (MkDocs Material), CI, docs deployment and release workflows.
- `boilflux.regimes`: boiling-regime analysis (MHF/DNB, rewetting/Leidenfrost point, regime maps,
  wetting-front velocity) with plots.
- `boilflux.operators` (extra `[ml]`): PyTorch FNO, BubbleML 2.0 range-read loader, autoregressive
  training with multi-step unrolling, rollout evaluation including heater heat flux.
- Physics guide (boiling curve, nucleate boiling, CHF/DNB, transition boiling, Leidenfrost, film boiling,
  quench rewetting) and modelling guide (CFD/OpenFOAM, PINNs, neural operators).
