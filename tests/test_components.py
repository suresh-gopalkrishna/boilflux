import numpy as np
import pytest

import quenchflux as qf


def test_material_interpolation_and_csv(tmp_path):
    m = qf.get_material("AA6082")
    assert m.k(20.0) == pytest.approx(172.0)
    assert m.k(1000.0) == pytest.approx(m.k_table[-1])  # constant extrapolation
    p = tmp_path / "steel.csv"
    p.write_text("T_C,k,rho,cp\n20,50,7850,460\n800,25,7600,640\n")
    s = qf.Material.from_csv(p)
    assert s.k(410.0) == pytest.approx(37.5)
    with pytest.raises(KeyError):
        qf.get_material("unobtainium")


def test_boiling_curve_continuity_and_derivative():
    bc = qf.SpechtBoilingCurve()
    assert bc.q_100 == pytest.approx(1.5e6, rel=0.03)  # Table 5 of the ATE 2026 paper
    for T in (100.0, bc.T_dnb, bc.T_leid):
        assert bc(np.array([T - 1e-6]))[0] == pytest.approx(bc(np.array([T + 1e-6]))[0], rel=1e-4)
    T = np.linspace(30, 600, 301)
    q, dq = bc.flux_and_derivative(T)
    h = 1e-4
    fd = (bc(T + h) - bc(T - h)) / (2 * h)
    keep = np.all(np.abs(T[:, None] - np.array([100.0, bc.T_dnb, bc.T_leid])) > 1e-3, axis=1)
    np.testing.assert_allclose(dq[keep], fd[keep], rtol=1e-5, atol=1.0)
    assert bc(np.array([bc.T_dnb]))[0] == pytest.approx(bc.q_max)


def test_plate_motion_and_camera_mapping():
    g = qf.PlateGrid(0.10, 0.02, 0.005, nx=50, ny=4, nz=2)
    mo = qf.PlateMotion(velocity=0.02, x0=0.0)
    assert mo.stagnation_point(0.05, 1.0) == pytest.approx(0.03)
    t = np.array([0.0, 1.0])
    # synthetic lab-frame frames of a plate whose temperature equals its plate x [mm]
    Xpx = np.linspace(-0.05, 0.2, 501)
    Ypx = np.linspace(-0.005, 0.025, 31)
    frames = np.stack([np.broadcast_to(mo.lab_to_plate(Xpx, tt) * 1e3, (31, 501)) for tt in t])
    data, mask = qf.lab_frames_to_plate(frames, Xpx, Ypx, t, g, mo)
    assert mask.all()
    np.testing.assert_allclose(data[1, 0], g.x * 1e3, atol=1e-9)
    vm = mo.view_mask(g, t, (0.0, 0.06))
    assert vm[0].sum() > vm[1].sum()  # part of the plate has left the view


def test_noise_estimate():
    rng = np.random.default_rng(3)
    t = np.linspace(0, 1, 400)[:, None, None]
    Y = 500 - 200 * t**2 + 1.5 * rng.standard_normal((400, 5, 7))
    assert qf.estimate_noise_sigma(Y) == pytest.approx(1.5, rel=0.08)


def test_dataset_roundtrip(tmp_path):
    sc = qf.JetQuenchScenario(nx=8, ny=4, nz=3, nt=5, refine=1)
    ds = sc.generate()
    p = ds.save(tmp_path / "d.npz")
    ds2 = qf.QuenchDataset.load(p)
    np.testing.assert_array_equal(ds.Y, ds2.Y)
    np.testing.assert_array_equal(ds.q_true, ds2.q_true)
    assert ds2.grid == ds.grid
    assert ds2.motion.velocity == sc.plate_velocity


def test_metrics_perfect():
    q = np.random.default_rng(0).random((4, 3, 2)) * 1e6
    m = qf.flux_metrics(q, q)
    assert m["r2"] == 1.0 and m["rmse_MW_m2"] == 0.0
