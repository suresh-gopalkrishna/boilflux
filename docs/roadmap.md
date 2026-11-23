# Roadmap

BoilFlux is the software home of a long-term research programme on inverse boiling heat transfer.
Every method published in this programme is released here, with its benchmark scores.

## v0.1 (current)

- [x] Nonlinear 3D finite-volume forward model; the 2D model of the ATE 2026 paper as a special case
- [x] Moving-plate formulation, camera field-of-view masks, IR frame resampling
- [x] Adjoint CGM with Tikhonov regularisation and the discrepancy principle; adjoint verified against finite differences
- [x] Specht boiling-curve synthetic generator and BF-Bench tier 2 cases
- [x] CLI, documentation website, CI
- [x] Boiling-regime analysis (q_max, T_DNB, T_Leid, rewetting time, regime maps, wetting-front speed)
- [x] Physics guide: nucleate boiling, CHF/DNB, transition, Leidenfrost, film boiling, quench rewetting
- [x] Neural-operator prototype for boiling flow (FNO on BubbleML 2.0, wall-heat-flux evaluation)

## Boiling-flow track (forward models)

- [ ] Operator benchmark on BubbleML pool/flow boiling: FNO, U-Net, DeepONet, transformer; error growth and wall-flux statistics
- [ ] Physics-constrained training (mass and energy conservation, interface jump conditions) and noise-injected rollouts
- [ ] OpenFOAM jet-quench cases with conjugate heat transfer (VOF + phase change and Euler–Euler wall boiling), with case files and data generator
- [ ] Operator surrogate of the **jet-quench** wall flux \(q_w(x,t)\) trained on OpenFOAM data and validated against BoilFlux reconstructions
- [ ] PINN/PIKAN baselines for the conduction inverse problem on BF-Bench

## v0.2: baselines and validation

- [ ] Beck's sequential function specification (1D and 3D) as a baseline
- [ ] Parametric boiling-curve inversion (\(q_\text{max}, T_\text{DNB}, T_\text{Leid}\)) with Levenberg–Marquardt
- [ ] Thermocouple data with sensor-lag correction (embedded thermocouples)
- [ ] First- and second-order (space-time) Tikhonov regularisation; L-curve and GCV
- [ ] BF-Bench tier 1: analytic cases (Beck 1D, the Huang–Wang 3D plate)

## v0.3: learned priors (paper B)

- [ ] Fourier neural operator and PCA-Net priors trained on synthetic quench twins
- [ ] Hybrid H1 (CGM warm-started from the operator) and H2 (Tikhonov centred on the operator prediction)
- [ ] Label-free consistency test \(\rho = \text{RMS}/\hat\sigma\) for out-of-distribution detection

## v0.4: verifier-gated agent (paper C)

- [ ] Typed tool interface for diagnostics, method selection and escalation
- [ ] Gates: discrepancy, residual whiteness, physical plausibility; VERIFIED / UNVERIFIED reports
- [ ] Optional LLM planner evaluated against the same gates

## v0.5: uncertainty quantification

- [ ] Laplace approximation using the adjoint machinery
- [ ] Credible bands on \(q(x,y,t)\) and on derived boiling-curve features
- [ ] Interfaces to hIPPYlib / CUQIpy

## v1.0: open experimental data (BF-Bench tier 4)

- [ ] First public IR quench dataset on hot metal plates (stationary and moving), with metadata, DOI and CC BY 4.0
- [ ] Leaderboard on the website; community submissions
- [ ] JOSS software paper

## Performance

- [ ] Reuse factorisations, multigrid preconditioning, optional GPU backend (CuPy / JAX)
- [ ] Unstructured FEM backend (FEniCSx) for non-rectangular specimens
