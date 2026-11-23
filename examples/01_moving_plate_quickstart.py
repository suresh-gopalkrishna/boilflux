"""Quickstart: synthetic 3D moving-plate quench -> IR data -> heat-flux reconstruction.

Run:  python examples/01_moving_plate_quickstart.py
"""

import boilflux as qf
from boilflux.plotting import save_all

# 1) A single round water jet on an AA6082 plate (120 x 60 x 5 mm) moving at 10 mm/s.
scenario = qf.JetQuenchScenario(
    material="AA6082", T0=500.0, plate_velocity=0.010, x_jet_lab=0.045,
    nx=48, ny=24, nz=5, dt=0.01, nt=200, noise_sigma=2.0,
)
data = scenario.generate()          # back-face IR temperatures + ground-truth flux
data.save("moving_plate_dataset.npz")

# 2) Reconstruct the boiling heat flux q(x, y, t) on the jet side.
solver = qf.AdjointCGM(data, options=qf.CGMOptions(max_iter=150, tau=1.02))
result = solver.solve()

# 3) Score against the ground truth (frames free of the end-time effect).
print(qf.flux_metrics(result.q[result.reliable], data.q_true[result.reliable]))
save_all(data, result, "quickstart_figures")
