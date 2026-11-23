"""Boiling-curve models used as ground truth in synthetic quench studies.

The default model is the parabolic boiling curve of Specht as used in
Gopalkrishna, Nallathambi & Specht, Applied Thermal Engineering 304 (2026)
132481, Eqs. (39)-(40):

* forced convection (T_s <= 100 degC):
  ``q = [14.6 (T_s+273) + 68.5 (T_f+273) - 21500] Qdot^0.33 (T_s - T_f)``
* nucleate boiling (100 < T_s <= T_DNB): parabola from q_100 up to q_max
* transition boiling (T_DNB < T_s <= T_Leid): parabola from q_max down to q_Leid
* film boiling (T_s > T_Leid): ``q = q_Leid``

Temperatures in degC, heat fluxes in W/m^2, Qdot in L/min.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np

__all__ = ["FFZ_DEFAULT", "IZ_DEFAULT", "SpechtBoilingCurve"]


@dataclass
class SpechtBoilingCurve:
    q_max: float = 6.5e6
    T_dnb: float = 200.0
    q_leid: float = 2.0e6
    T_leid: float = 450.0
    T_f: float = 20.0
    Q_dot: float = 100.0

    def __post_init__(self) -> None:
        # accept numeric strings, e.g. "6.5e6" read by YAML 1.1 parsers
        for name in ("q_max", "T_dnb", "q_leid", "T_leid", "T_f", "Q_dot"):
            setattr(self, name, float(getattr(self, name)))
        if not 100.0 < self.T_dnb < self.T_leid:
            raise ValueError("require 100 < T_dnb < T_leid")

    def _fc_coeff(self, Ts):
        return (14.6 * (Ts + 273.0) + 68.5 * (self.T_f + 273.0) - 21500.0) * self.Q_dot**0.33

    @property
    def q_100(self) -> float:
        return float(self._fc_coeff(100.0) * (100.0 - self.T_f))

    def __call__(self, Ts: np.ndarray) -> np.ndarray:
        return self.flux_and_derivative(Ts)[0]

    def flux_and_derivative(self, Ts: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Heat flux q(T_s) and dq/dT_s."""
        Ts = np.asarray(Ts, dtype=float)
        q = np.empty_like(Ts)
        dq = np.empty_like(Ts)
        q100 = self.q_100
        # forced convection
        m = Ts <= 100.0
        a = self._fc_coeff(Ts[m])
        q[m] = a * (Ts[m] - self.T_f)
        dq[m] = 14.6 * self.Q_dot**0.33 * (Ts[m] - self.T_f) + a
        # nucleate boiling
        m = (Ts > 100.0) & (Ts <= self.T_dnb)
        u = (self.T_dnb - Ts[m]) / (self.T_dnb - 100.0)
        q[m] = self.q_max - u**2 * (self.q_max - q100)
        dq[m] = 2.0 * u * (self.q_max - q100) / (self.T_dnb - 100.0)
        # transition boiling
        m = (Ts > self.T_dnb) & (Ts <= self.T_leid)
        u = (Ts[m] - self.T_dnb) / (self.T_leid - self.T_dnb)
        q[m] = self.q_max - u**2 * (self.q_max - self.q_leid)
        dq[m] = -2.0 * u * (self.q_max - self.q_leid) / (self.T_leid - self.T_dnb)
        # film boiling
        m = Ts > self.T_leid
        q[m] = self.q_leid
        dq[m] = 0.0
        return q, dq

    def to_dict(self) -> dict:
        return asdict(self)


#: Impingement-zone parameters (Table 5 of the ATE 2026 paper)
IZ_DEFAULT = SpechtBoilingCurve(q_max=6.5e6, T_dnb=200.0, q_leid=2.0e6, T_leid=450.0)
#: Free-falling-zone parameters (Table 5 of the ATE 2026 paper)
FFZ_DEFAULT = SpechtBoilingCurve(q_max=3.5e6, T_dnb=200.0, q_leid=1.0e6, T_leid=350.0)
