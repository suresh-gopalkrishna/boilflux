"""Error metrics for reconstructed heat-flux fields (benchmark scoring)."""

from __future__ import annotations

import numpy as np

__all__ = ["flux_metrics"]


def flux_metrics(q_est: np.ndarray, q_true: np.ndarray, mape_threshold: float = 0.05) -> dict:
    """Compare a reconstructed flux field with the ground truth.

    Returns RMSE [MW/m^2], relative L2 error, coefficient of determination R^2,
    MAPE [%] over points where ``|q_true| > mape_threshold * max|q_true|``, and
    the relative error of the maximum heat flux (MHF) [%].
    """
    e = np.asarray(q_est, float) - np.asarray(q_true, float)
    qt = np.asarray(q_true, float)
    rmse = float(np.sqrt(np.mean(e**2)))
    rel_l2 = float(np.linalg.norm(e) / max(np.linalg.norm(qt), 1e-300))
    ss_res = float(np.sum(e**2))
    ss_tot = float(np.sum((qt - qt.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    sel = np.abs(qt) > mape_threshold * np.abs(qt).max()
    mape = float(100.0 * np.mean(np.abs(e[sel] / qt[sel]))) if sel.any() else float("nan")
    mhf_err = float(100.0 * (np.max(q_est) - np.max(qt)) / np.max(qt))
    return {
        "rmse_MW_m2": rmse / 1e6,
        "rel_l2": rel_l2,
        "r2": r2,
        "mape_percent": mape,
        "mhf_error_percent": mhf_err,
        "q_max_true_MW_m2": float(np.max(qt)) / 1e6,
        "q_max_est_MW_m2": float(np.max(q_est)) / 1e6,
    }
