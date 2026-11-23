import numpy as np
import pytest

import quenchflux as qf


def small_model(material=None, nt=20, dt=0.01, solver="direct", loss=None, ny=4):
    g = qf.PlateGrid(0.02, 0.01, 0.004, nx=8, ny=ny, nz=4)
    mat = material or qf.Material.constant(k=50.0, rho=8000.0, cp=500.0)
    return qf.ForwardModel(g, mat, dt, nt, 300.0, solver=solver, top_loss=loss)


def test_energy_conservation_constant_flux():
    fm = small_model()
    g = fm.grid
    q = np.full((fm.nt, *g.face_shape), 2.0e6)
    traj = fm.run(q=q)
    rc = 8000.0 * 500.0
    dE = rc * g.cell_volume * (traj.T[-1] - traj.T[0]).sum()
    removed = 2.0e6 * g.Lx * g.Ly * fm.dt * fm.nt
    assert dE == pytest.approx(-removed, rel=1e-9)


def test_energy_conservation_nonlinear_material():
    fm = small_model(material=qf.get_material("nickel"))
    g = fm.grid
    traj = fm.run(q=np.zeros((fm.nt, *g.face_shape)))
    np.testing.assert_allclose(traj.T, 300.0, atol=1e-8)  # insulated: no change


@pytest.mark.parametrize("solver", ["direct", "cg"])
def test_solvers_agree(solver):
    fm = small_model(solver=solver, material=qf.get_material("aa6082"))
    rng = np.random.default_rng(1)
    q = 1e6 * rng.random((fm.nt, *fm.grid.face_shape))
    ref = small_model(solver="direct", material=qf.get_material("aa6082")).run(q=q).T
    np.testing.assert_allclose(fm.run(q=q).T, ref, rtol=0, atol=1e-7)


def test_one_dimensional_symmetry():
    # uniform flux must give a temperature field independent of x and y
    fm = small_model()
    q = np.full((fm.nt, *fm.grid.face_shape), 1.0e6)
    f = fm.run(q=q).field(fm.nt)
    assert np.ptp(f, axis=(1, 2)).max() < 1e-9
    assert np.all(np.diff(f[:, 0, 0]) > 0)  # colder at the quenched face


def test_semi_infinite_surface_temperature():
    # thick plate, short time: compare with the semi-infinite analytic solution
    k, rho, cp, q0 = 50.0, 8000.0, 500.0, 1.0e6
    g = qf.PlateGrid(0.002, 0.002, 0.02, nx=1, ny=1, nz=200)
    fm = qf.ForwardModel(g, qf.Material.constant(k, rho, cp), dt=1e-3, nt=500, T0=500.0)
    traj = fm.run(q=np.full((fm.nt, 1, 1), q0))
    alpha = k / (rho * cp)
    t = traj.t[-1]
    z = g.z[5]
    from scipy.special import erfc

    exact = 500.0 - 2 * q0 / k * (np.sqrt(alpha * t / np.pi) * np.exp(-z**2 / (4 * alpha * t))
                                   - z / 2 * erfc(z / (2 * np.sqrt(alpha * t))))
    assert traj.field(fm.nt)[5, 0, 0] == pytest.approx(exact, abs=0.5)


def test_surface_loss_cools_top():
    fm = small_model(loss=qf.SurfaceLoss(h=50.0, emissivity=0.9))
    traj = fm.run(q=np.zeros((fm.nt, *fm.grid.face_shape)))
    assert traj.top_face()[-1].mean() < 300.0
