import json

import numpy as np
import pytest

torch = pytest.importorskip("torch")

from boilflux.operators.fno import FNO2d  # noqa: E402
from boilflux.operators.train import (  # noqa: E402
    BoilingSurrogate,
    TrainConfig,
    evaluate_rollout,
    train_surrogate,
    wall_heat_flux,
)


def test_fno_shapes_and_resolution():
    net = FNO2d(5, 3, width=8, modes=4, layers=2, padding=4)
    assert net(torch.randn(2, 5, 24, 16)).shape == (2, 3, 24, 16)
    assert net(torch.randn(1, 5, 32, 20)).shape == (1, 3, 32, 20)  # resolution-agnostic


def _fake_sim(path, twall, nt=20, ny=16, nx=12):
    t = np.arange(nt) * 0.4
    y = (np.arange(ny) + 0.5) * 0.1
    x = (np.arange(nx) + 0.5) * 0.1 - 0.6
    yy, xx = np.meshgrid(y, x, indexing="ij")
    T = twall - 30 * (1 - np.exp(-yy / 0.3))[None] + np.sin(t)[:, None, None]
    phi = 0.3 - np.hypot(xx, yy - 0.3 - 0.02 * t[:, None, None])
    np.savez(path, x=x, y=y, t=t, meta=json.dumps({"dt": 0.4}), wall_temp=float(twall), source="fake",
             temperature=T.astype(np.float32), velx=np.zeros((nt, ny, nx), np.float32),
             vely=0.1 * np.ones((nt, ny, nx), np.float32), dfun=phi.astype(np.float32))
    return path


def test_train_rollout_and_save(tmp_path):
    files = [_fake_sim(tmp_path / f"s{T}.npz", T) for T in (90, 100)]
    cfg = TrainConfig(history=2, width=8, modes=4, layers=2, epochs=2, batch_size=4)
    model, info = train_surrogate(files, cfg, log=lambda *_: None)
    assert np.isfinite(info["loss"]).all()
    ev = evaluate_rollout(model, files[0], steps=5)
    assert ev["pred"].shape == (5, 4, 16, 12)
    model.save(tmp_path / "m.pt")
    m2 = BoilingSurrogate.load(tmp_path / "m.pt")
    s = np.load(files[0])
    st = np.stack([s[k] for k in ("temperature", "velx", "vely", "dfun")], 1)[:2]
    np.testing.assert_allclose(m2.rollout(st, 90.0, 2), model.rollout(st, 90.0, 2), atol=1e-5)


def test_wall_heat_flux():
    T = np.full((3, 4, 5), 80.0)
    q = wall_heat_flux(T, 90.0, np.array([0.05, 0.15, 0.25, 0.35]))
    np.testing.assert_allclose(q, 200.0)
