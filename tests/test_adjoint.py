"""Adjoint and sensitivity consistency: the core correctness test of the inverse solver."""

import numpy as np
import pytest

import boilflux as qf


def setup(material):
    g = qf.PlateGrid(0.02, 0.01, 0.004, nx=6, ny=3, nz=4)
    fm = qf.ForwardModel(g, material, 0.01, 12, 400.0, solver="direct")
    rng = np.random.default_rng(0)
    q = 1e6 * (1 + rng.random((fm.nt, *g.face_shape)))
    Y = fm.run(q=q * 0.8).top_face()[1:] + rng.standard_normal((fm.nt, *g.face_shape))
    return fm, q, Y, rng


def misfit(fm, q, Y):
    return 0.5 * np.sum((fm.run(q=q).top_face()[1:] - Y) ** 2)


def test_gradient_matches_finite_differences_linear():
    fm, q, Y, rng = setup(qf.Material.constant(60.0, 7800.0, 480.0))
    traj = fm.run(q=q)
    r = traj.top_face()[1:] - Y
    grad = fm.adjoint(traj, r)
    dq = 1e5 * rng.standard_normal(q.shape)
    eps = 1e-3
    fd = (misfit(fm, q + eps * dq, Y) - misfit(fm, q - eps * dq, Y)) / (2 * eps)
    assert np.sum(grad * dq) == pytest.approx(fd, rel=1e-6)


def test_sensitivity_is_transpose_of_adjoint():
    fm, q, Y, rng = setup(qf.get_material("nickel"))  # nonlinear, frozen coefficients
    traj = fm.run(q=q)
    dq = rng.standard_normal(q.shape)
    w = rng.standard_normal(q.shape)
    lhs = np.sum(w * fm.sensitivity(traj, dq))
    rhs = np.sum(fm.adjoint(traj, w) * dq)  # adjoint returns -B^T mu, i.e. (dT/dq)^T w
    assert lhs == pytest.approx(rhs, rel=1e-8)


def test_gradient_nonlinear_is_close():
    fm, q, Y, rng = setup(qf.get_material("aa6082"))
    traj = fm.run(q=q)
    grad = fm.adjoint(traj, traj.top_face()[1:] - Y)
    dq = 1e5 * rng.standard_normal(q.shape)
    eps = 1e-2
    fd = (misfit(fm, q + eps * dq, Y) - misfit(fm, q - eps * dq, Y)) / (2 * eps)
    # frozen-property gradient: small, bounded inconsistency for weak nonlinearity
    assert np.sum(grad * dq) == pytest.approx(fd, rel=5e-2)
