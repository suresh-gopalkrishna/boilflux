"""Reference benchmark cases (QF-Bench tier 2: synthetic nonlinear quench twins).

Each case is a :class:`~quenchflux.synthetic.JetQuenchScenario` with fixed
parameters and seed, so scores are reproducible across machines. Scores are
computed on the frames unaffected by the end-time effect.
"""

from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path

from .boiling import SpechtBoilingCurve
from .inverse import AdjointCGM, CGMOptions
from .metrics import flux_metrics
from .synthetic import JetQuenchScenario

__all__ = ["CASES", "get_case", "run_case"]


CASES: dict[str, dict] = {
    "S2D-AA6082-stationary": dict(
        description="2D (length x thickness) stationary full-jet quench of AA6082, "
                    "the synthetic configuration of the ATE 2026 paper (shortened plate).",
        scenario=dict(material="AA6082", Lx=0.100, Ly=0.010, Lz=0.005, nx=80, ny=1, nz=6,
                      dt=0.01, nt=300, T0=500.0, plate_velocity=0.0, x_jet_lab=0.0,
                      r_iz=0.009, r_ffz=0.011, r_wet=10.0, wet_edge=0.004, noise_sigma=2.0),
    ),
    "M3D-AA6082": dict(
        description="3D moving plate (10 mm/s), single round jet, AA6082, T0 = 500 degC.",
        scenario=dict(material="AA6082", plate_velocity=0.01, T0=500.0),
    ),
    "M3D-Ni": dict(
        description="3D moving plate (5 mm/s), single round jet, pure nickel, T0 = 800 degC; "
                    "slow rewetting, film boiling dominates the first seconds.",
        scenario=dict(
            material="Nickel", plate_velocity=0.005, x_jet_lab=0.07, T0=800.0, nt=500,
            iz=asdict(SpechtBoilingCurve(q_max=12e6, T_dnb=250.0, q_leid=2.5e6, T_leid=600.0)),
            ffz=asdict(SpechtBoilingCurve(q_max=5e6, T_dnb=250.0, q_leid=1.0e6, T_leid=450.0)),
        ),
    ),
    "M3D-AA6082-window": dict(
        description="As M3D-AA6082 but plate faster (20 mm/s) and the IR camera sees only "
                    "the lab window 10-100 mm: plate regions leave the field of view.",
        scenario=dict(material="AA6082", plate_velocity=0.02, T0=500.0,
                      camera_view_lab=(0.010, 0.100)),
    ),
}


def get_case(name: str) -> JetQuenchScenario:
    if name not in CASES:
        raise KeyError(f"unknown case '{name}'. Available: {list(CASES)}")
    return JetQuenchScenario.from_dict(CASES[name]["scenario"])


def run_case(name: str, outdir: str | Path | None = None, options: CGMOptions | None = None,
             plots: bool = True, verbose: bool = True) -> dict:
    """Generate the case, invert it with the adjoint CGM and return the scores."""
    scen = get_case(name)
    data = scen.generate()
    opts = options or CGMOptions(verbose=verbose)
    solver = AdjointCGM(data, options=opts)
    res = solver.solve()
    sl = res.reliable
    scores = flux_metrics(res.q[sl], data.q_true[sl])
    scores.update(case=name, iterations=len(res.history["alpha"]), converged=res.converged,
                  runtime_s=res.runtime_s, sigma_K=res.sigma,
                  reliable_until_s=res.n_reliable * data.dt, method="adjoint-CGM (QuenchFlux)")
    if outdir is not None:
        out = Path(outdir) / name
        out.mkdir(parents=True, exist_ok=True)
        data.save(out / "dataset.npz")
        res.save(out / "result.npz")
        (out / "scores.json").write_text(json.dumps(scores, indent=2))
        if plots:
            from .plotting import save_all

            save_all(data, res, out)
    return scores
