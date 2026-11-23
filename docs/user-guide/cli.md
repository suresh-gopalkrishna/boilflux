# Command line

```text
boilflux info                         version, materials, benchmark cases
boilflux demo   [--out DIR]           3D moving-plate example with figures
boilflux bench  [--case NAME|all]     run BF-Bench cases and print the score table
boilflux synth  CONFIG.yaml -o D.npz  generate a synthetic dataset
boilflux invert D.npz [-c CFG.yaml] [-o DIR]   reconstruct the heat flux
```

## Scenario file (`synth`)

```yaml
scenario:
  material: AA6082
  Lx: 0.120            # plate length [m] (direction of motion)
  Ly: 0.060
  Lz: 0.005
  nx: 48
  ny: 24
  nz: 5
  dt: 0.01             # frame interval [s]
  nt: 200
  T0: 500.0
  plate_velocity: 0.010
  x_jet_lab: 0.045
  noise_sigma: 2.0
  refine: 2            # truth computed on a 2x finer grid
  iz:  {q_max: 6.5e+6, T_dnb: 200.0, q_leid: 2.0e+6, T_leid: 450.0}
  ffz: {q_max: 3.5e+6, T_dnb: 200.0, q_leid: 1.0e+6, T_leid: 350.0}
```

## Inversion file (`invert`)

```yaml
cgm:
  max_iter: 200
  tau: 1.02            # discrepancy safety factor
  gamma: 0.0           # Tikhonov weight [K^2 per (MW/m^2)^2]
material_csv: my_alloy.csv          # optional
top_loss: {h: 10.0, emissivity: 0.92, T_amb: 25.0}   # optional
```

Both files are in `examples/configs/`.
