# Getting started

## Install

BoilFlux needs Python 3.10 or newer and depends only on NumPy, SciPy, Matplotlib and PyYAML.

=== "pip (release)"

    ```bash
    pip install boilflux
    ```

=== "latest from GitHub"

    ```bash
    pip install git+https://github.com/boilflux/boilflux.git
    ```

=== "development"

    ```bash
    git clone https://github.com/boilflux/boilflux.git
    cd boilflux
    pip install -e ".[dev,docs]"
    pytest
    ```

Check the installation:

```bash
boilflux info
```

## Your first reconstruction

The demo generates a synthetic moving-plate quench, inverts the back-face IR data and writes
figures and scores. It runs in about a minute on a laptop.

```bash
boilflux demo --out demo
```

This produces `demo/M3D-AA6082/` containing `flux_maps.png`, `centerline_xt.png`,
`histories.png`, `boiling_curves.png`, `convergence.png`, `scores.json`, and the dataset and
result as `.npz` files.

## The same in Python

```python
import boilflux as qf
from boilflux.plotting import save_all

# 1. a single round jet on an AA6082 plate moving at 10 mm/s
scenario = qf.JetQuenchScenario(material="AA6082", T0=500.0, plate_velocity=0.010,
                                nx=48, ny=24, nz=5, dt=0.01, nt=200, noise_sigma=2.0)
data = scenario.generate()

# 2. adjoint conjugate-gradient reconstruction
result = qf.AdjointCGM(data, options=qf.CGMOptions(max_iter=150)).solve()

# 3. scores and figures
print(qf.flux_metrics(result.q[result.reliable], data.q_true[result.reliable]))
save_all(data, result, "figures")
```

`result.q` has shape `(nt, ny, nx)` and gives the heat flux in W/m² on the quenched face at each
frame. `result.T_surface` gives the matching quenched-face temperature, so `q` can be plotted
against `T_surface` to obtain local boiling curves.

!!! warning "End-time effect"
    The adjoint vanishes at the last frame, so the flux recovered in roughly the last
    \(0.5\,L_z^2/\alpha\) seconds is biased towards zero. Use `result.reliable` (a time slice) to
    discard these frames. When planning experiments, record a little longer than the interval you need.

## Next steps

- Process your own camera data: [Data and IR cameras](user-guide/data.md)
- Understand the method: [Theory](theory/index.md)
- Compare methods: [BF-Bench](benchmark.md)
