import numpy as np

import boilflux as qf


def test_moving_plate_inversion_recovers_flux():
    sc = qf.JetQuenchScenario(nx=24, ny=12, nz=4, nt=60, dt=0.02, refine=2, seed=1)
    data = sc.generate()
    res = qf.AdjointCGM(data, options=qf.CGMOptions(max_iter=60, verbose=False)).solve()
    assert res.converged
    sl = res.reliable
    m = qf.flux_metrics(res.q[sl], data.q_true[sl])
    assert m["r2"] > 0.8
    assert m["rel_l2"] < 0.35
    rms = np.sqrt(res.history["misfit"][-1] / data.mask.sum())
    assert rms <= 1.03 * data.sigma


def test_camera_window_mask_is_respected():
    sc = qf.JetQuenchScenario(nx=20, ny=6, nz=3, nt=30, dt=0.02, refine=1,
                              plate_velocity=0.03, camera_view_lab=(0.0, 0.08))
    data = sc.generate()
    assert not data.mask.all()
    data.Y[~data.mask] = 1e9  # garbage outside the view must not matter
    res = qf.AdjointCGM(data, options=qf.CGMOptions(max_iter=10, verbose=False)).solve()
    assert np.all(np.isfinite(res.q))
    assert np.abs(res.q).max() < 5e7


def test_cli_info(capsys):
    from boilflux.cli import main

    assert main(["info"]) == 0
    assert "M3D-AA6082" in capsys.readouterr().out
