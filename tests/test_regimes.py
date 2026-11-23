import numpy as np
import pytest

import quenchflux as qf
from quenchflux.regimes import Regime, extract_features, feature_maps, wetting_front


def test_features_of_specht_curve():
    bc = qf.SpechtBoilingCurve()  # q_max 6.5 MW/m2 at 200 C, Leidenfrost 450 C
    t = np.linspace(0, 5, 1001)
    Ts = 600 - 110 * t  # monotonic cooling 600 -> 50 C
    f = extract_features(bc(Ts), Ts, t, smooth=1, frac=0.02)
    assert f.q_max == pytest.approx(6.5e6, rel=1e-3)
    assert f.T_dnb == pytest.approx(200.0, abs=1.0)
    assert f.T_leid == pytest.approx(450.0, abs=8.0)
    assert f.wetted
    seq = [k for i, k in enumerate(f.regime) if i == 0 or k != f.regime[i - 1]]
    assert seq == [Regime.FILM, Regime.TRANSITION, Regime.NUCLEATE, Regime.CONVECTION]


def test_not_wetted_record():
    t = np.linspace(0, 1, 100)
    f = extract_features(np.linspace(1e6, 2e6, 100), np.full(100, 700.0), t)
    assert not f.wetted


def test_maps_and_front():
    t = np.linspace(0, 2, 201)
    x = np.linspace(0, 0.03, 11)
    bc = qf.SpechtBoilingCurve()
    Ts = 600 - 400 * np.clip(t[:, None] - 20 * x[None, :], 0, None)  # front moving at 0.05 m/s
    q = bc(Ts)[:, None, :]
    maps = feature_maps(q, Ts[:, None, :], t, smooth=1)
    fr = wetting_front(maps["t_wet"][0], maps["wetted"][0], x, t)
    assert fr["v_downstream"] == pytest.approx(0.05, rel=0.15)
    assert fr["v_upstream"] == 0.0
