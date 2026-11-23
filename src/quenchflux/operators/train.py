"""Train and evaluate an autoregressive FNO surrogate of boiling flow.

State per frame: temperature ``T``, velocities ``u, v`` and the level set
``phi`` (``dfun``, > 0 in vapour). The operator maps the last ``history``
frames plus the heater temperature to the increment of the next frame:

    s^{n+1} = s^n + G_theta(s^{n-h+1..n}, T_wall)

The physically most important output for quenching and boiling studies is the
**wall heat flux**, evaluated here from the near-wall temperature gradient,
q_w ~ k_l (T_wall - T_1) / y_1, and compared with the simulation along
autoregressive rollouts.
"""

from __future__ import annotations

import json
import time
from collections.abc import Sequence
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import torch

from .bubbleml import load_fetched
from .fno import FNO2d

__all__ = ["TrainConfig", "BoilingSurrogate", "wall_heat_flux", "train_surrogate"]

FIELDS = ("temperature", "velx", "vely", "dfun")
PHI_CLIP = 1.0  # level set clipped to +/- 1 (only the interface neighbourhood matters)


@dataclass
class TrainConfig:
    history: int = 4
    width: int = 32
    modes: int = 12
    layers: int = 4
    epochs: int = 30
    batch_size: int = 8
    lr: float = 1e-3
    weight_decay: float = 1e-5
    rollout_unroll: int = 2  # pushforward-style multi-step training
    seed: int = 0


def _stack(sim: dict) -> np.ndarray:
    s = np.stack([sim[f] for f in FIELDS], axis=1).astype(np.float32)  # (t, 4, ny, nx)
    s[:, 3] = np.clip(s[:, 3], -PHI_CLIP, PHI_CLIP)
    return s


def wall_heat_flux(T: np.ndarray, T_wall: float, y: np.ndarray) -> np.ndarray:
    """Dimensionless wall heat flux time series (mean over the heater).

    ``T`` has shape (t, ny, nx) with row 0 adjacent to the heater at y = 0.
    Uses the one-sided gradient (T_wall - T_1) / y_1. Vapour-covered wall
    cells contribute with the same formula (their low conductivity is
    reflected by the small gradient in the data).
    """
    return ((T_wall - T[:, 0, :]) / y[0]).mean(axis=-1)


class BoilingSurrogate:
    """FNO surrogate plus the normalisation statistics it was trained with."""

    def __init__(self, cfg: TrainConfig, mean: np.ndarray, std: np.ndarray, twall_stats: tuple) -> None:
        self.cfg = cfg
        self.mean = mean.astype(np.float32)
        self.std = std.astype(np.float32)
        self.dmean = np.zeros_like(self.mean)
        self.dstd = np.ones_like(self.std)
        self.twall_stats = twall_stats
        c = len(FIELDS)
        self.net = FNO2d(cfg.history * c + 1, c, cfg.width, cfg.modes, cfg.layers)

    # --- tensors --------------------------------------------------------
    def _inp(self, hist: torch.Tensor, twall: torch.Tensor) -> torch.Tensor:
        # hist: (b, h, c, ny, nx) normalised
        b, h, c, ny, nx = hist.shape
        tw = ((twall - self.twall_stats[0]) / self.twall_stats[1]).view(b, 1, 1, 1).expand(b, 1, ny, nx)
        return torch.cat([hist.reshape(b, h * c, ny, nx), tw], dim=1)

    def _norm(self, s):
        return (s - torch.from_numpy(self.mean).view(1, -1, 1, 1)) / torch.from_numpy(self.std).view(1, -1, 1, 1)

    def _denorm(self, s):
        return s * torch.from_numpy(self.std).view(1, -1, 1, 1) + torch.from_numpy(self.mean).view(1, -1, 1, 1)

    def step(self, hist_n: torch.Tensor, twall: torch.Tensor) -> torch.Tensor:
        """One step in normalised space: returns the next normalised frame."""
        d = self.net(self._inp(hist_n, twall))
        d = d * torch.from_numpy(self.dstd).view(1, -1, 1, 1) + torch.from_numpy(self.dmean).view(1, -1, 1, 1)
        return hist_n[:, -1] + d

    @torch.no_grad()
    def rollout(self, states: np.ndarray, twall: float, steps: int) -> np.ndarray:
        """Autoregressive prediction from the first ``history`` frames of ``states``."""
        self.net.eval()
        h = self.cfg.history
        hist = self._norm(torch.from_numpy(_clip_phi(states[:h].copy()))).unsqueeze(0)
        tw = torch.tensor([twall], dtype=torch.float32)
        out = []
        for _ in range(steps):
            nxt = self.step(hist, tw)
            nxt[:, 3] = torch.clamp(nxt[:, 3], *self._phi_bounds())
            out.append(self._denorm(nxt)[0].numpy())
            hist = torch.cat([hist[:, 1:], nxt.unsqueeze(1)], dim=1)
        return np.stack(out)

    def _phi_bounds(self):
        m, s = float(self.mean[3]), float(self.std[3])
        return ((-PHI_CLIP - m) / s, (PHI_CLIP - m) / s)

    # --- persistence ----------------------------------------------------
    def save(self, path: str | Path) -> None:
        torch.save({"cfg": asdict(self.cfg), "mean": self.mean, "std": self.std, "dmean": self.dmean,
                    "dstd": self.dstd, "twall_stats": self.twall_stats, "state": self.net.state_dict()}, path)

    @classmethod
    def load(cls, path: str | Path) -> BoilingSurrogate:
        ck = torch.load(path, weights_only=False)
        m = cls(TrainConfig(**ck["cfg"]), ck["mean"], ck["std"], tuple(ck["twall_stats"]))
        m.dmean, m.dstd = ck["dmean"], ck["dstd"]
        m.net.load_state_dict(ck["state"])
        return m


def _clip_phi(s: np.ndarray) -> np.ndarray:
    s[:, 3] = np.clip(s[:, 3], -PHI_CLIP, PHI_CLIP)
    return s


def train_surrogate(files: Sequence[str | Path], cfg: TrainConfig | None = None,
                    log=print) -> tuple[BoilingSurrogate, dict]:
    """Train on fetched BubbleML windows (see :func:`fetch_simulation`)."""
    cfg = cfg or TrainConfig()
    torch.manual_seed(cfg.seed)
    rng = np.random.default_rng(cfg.seed)
    sims = [load_fetched(f) for f in files]
    states = [_stack(s) for s in sims]
    twalls = np.array([float(s["wall_temp"]) for s in sims], dtype=np.float32)
    allf = np.concatenate(states)
    mean = allf.mean(axis=(0, 2, 3))
    std = allf.std(axis=(0, 2, 3)) + 1e-6
    model = BoilingSurrogate(cfg, mean, std, (float(twalls.mean()), float(twalls.std() + 1e-6)))
    norm = [(s - mean[None, :, None, None]) / std[None, :, None, None] for s in states]
    d = np.concatenate([n[1:] - n[:-1] for n in norm])
    model.dmean = d.mean(axis=(0, 2, 3)).astype(np.float32)
    model.dstd = (d.std(axis=(0, 2, 3)) + 1e-6).astype(np.float32)
    h, k = cfg.history, cfg.rollout_unroll
    index = [(i, t) for i, n in enumerate(norm) for t in range(h, n.shape[0] - k + 1)]
    tens = [torch.from_numpy(n.astype(np.float32)) for n in norm]
    opt = torch.optim.AdamW(model.net.parameters(), lr=cfg.lr, weight_decay=cfg.weight_decay)
    nb = int(np.ceil(len(index) / cfg.batch_size))
    sched = torch.optim.lr_scheduler.OneCycleLR(opt, max_lr=cfg.lr, total_steps=cfg.epochs * nb)
    w = torch.tensor([1.0, 1.0, 1.0, 2.0]).view(1, -1, 1, 1)  # emphasise the interface
    hist_log = []
    t0 = time.perf_counter()
    for ep in range(cfg.epochs):
        model.net.train()
        perm = rng.permutation(len(index))
        tot = 0.0
        for bi in range(nb):
            sel = [index[j] for j in perm[bi * cfg.batch_size:(bi + 1) * cfg.batch_size]]
            hist = torch.stack([tens[i][t - h:t] for i, t in sel])
            tw = torch.from_numpy(twalls[[i for i, _ in sel]])
            loss = 0.0
            for u in range(k):
                tgt = torch.stack([tens[i][t + u] for i, t in sel])
                pred = model.step(hist, tw)
                loss = loss + (w * (pred - tgt) ** 2).mean() / k
                hist = torch.cat([hist[:, 1:], pred.unsqueeze(1)], dim=1)
            opt.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.net.parameters(), 1.0)
            opt.step()
            sched.step()
            tot += float(loss.detach())
        hist_log.append(tot / nb)
        log(f"[FNO] epoch {ep + 1:3d}/{cfg.epochs}  loss {tot / nb:.4e}  ({time.perf_counter() - t0:.0f} s)")
    info = {"loss": hist_log, "n_samples": len(index), "train_files": [str(f) for f in files],
            "train_time_s": time.perf_counter() - t0,
            "n_params": int(sum(p.numel() for p in model.net.parameters()))}
    return model, info


def evaluate_rollout(model: BoilingSurrogate, file: str | Path, steps: int = 60, start: int = 0) -> dict:
    """Roll out on a held-out simulation; return error metrics and arrays for plotting."""
    sim = load_fetched(file)
    s = _stack(sim)
    h = model.cfg.history
    steps = min(steps, s.shape[0] - start - h)
    pred = model.rollout(s[start:start + h], float(sim["wall_temp"]), steps)
    true = s[start + h:start + h + steps]
    rel = np.linalg.norm((pred - true).reshape(steps, 4, -1), axis=-1) / (
        np.linalg.norm((true - model.mean[None, :, None, None]).reshape(steps, 4, -1), axis=-1) + 1e-12)
    vap_p, vap_t = pred[:, 3] > 0, true[:, 3] > 0
    inter = (vap_p & vap_t).sum(axis=(1, 2))
    union = (vap_p | vap_t).sum(axis=(1, 2))
    iou = np.where(union > 0, inter / np.maximum(union, 1), 1.0)
    Tw = float(sim["wall_temp"])
    qp = wall_heat_flux(pred[:, 0], Tw, sim["y"])
    qt = wall_heat_flux(true[:, 0], Tw, sim["y"])
    t = sim["t"][start + h:start + h + steps]
    return {"rel_l2": rel, "iou_vapour": iou, "q_pred": qp, "q_true": qt, "t": t,
            "pred": pred, "true": true, "x": sim["x"], "y": sim["y"], "wall_temp": Tw,
            "q_rel_err_mean": float(np.mean(np.abs(qp - qt)) / np.mean(np.abs(qt)))}


def summary(ev: dict) -> dict:
    r = ev["rel_l2"]
    return {
        "wall_temp": ev["wall_temp"],
        "steps": int(r.shape[0]),
        "rel_l2_T_first5": float(r[:5, 0].mean()), "rel_l2_T_all": float(r[:, 0].mean()),
        "rel_l2_u_all": float(r[:, 1].mean()), "rel_l2_v_all": float(r[:, 2].mean()),
        "iou_vapour_mean": float(ev["iou_vapour"].mean()),
        "wall_heat_flux_mean_rel_err": ev["q_rel_err_mean"],
    }


def plot_rollout(ev: dict, outdir: str | Path, steps_to_show: Sequence[int] = (0, 10, 25, 50)) -> list[Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    n = ev["pred"].shape[0]
    show = [s for s in steps_to_show if s < n]
    ext = [ev["x"][0], ev["x"][-1], ev["y"][0], ev["y"][-1]]
    vmin, vmax = ev["true"][:, 0].min(), ev["true"][:, 0].max()
    fig, ax = plt.subplots(2, len(show), figsize=(2.3 * len(show), 6.2), constrained_layout=True, squeeze=False)
    im = None
    for j, s in enumerate(show):
        for i, (lab, arr) in enumerate((("simulation", ev["true"]), ("FNO", ev["pred"]))):
            a = ax[i, j]
            im = a.imshow(arr[s, 0], origin="lower", extent=ext, cmap="inferno", vmin=vmin, vmax=vmax)
            a.contour(ev["x"], ev["y"], arr[s, 3], levels=[0.0], colors="c", linewidths=0.9)
            a.set_title(f"{lab}, t = {ev['t'][s]:.1f}", fontsize=8)
            a.set_xticks([])
            a.set_yticks([])
    fig.colorbar(im, ax=ax, label="temperature", shrink=0.8)
    p1 = out / "rollout_snapshots.png"
    fig.savefig(p1, dpi=150)
    plt.close(fig)
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.2), constrained_layout=True)
    ax[0].plot(ev["t"], ev["q_true"], "k-", lw=1.4, label="simulation")
    ax[0].plot(ev["t"], ev["q_pred"], "--", color="#f26b1d", lw=1.4, label="FNO rollout")
    ax[0].set_xlabel("t (dimensionless)")
    ax[0].set_ylabel("wall heat flux (dimensionless)")
    ax[0].legend()
    ax[0].set_title("Heater heat flux during bubble growth and departure")
    ax[1].plot(ev["t"], ev["rel_l2"][:, 0], label="T")
    ax[1].plot(ev["t"], ev["rel_l2"][:, 1], label="u")
    ax[1].plot(ev["t"], ev["rel_l2"][:, 2], label="v")
    ax[1].plot(ev["t"], 1 - ev["iou_vapour"], label="1 - IoU(vapour)")
    ax[1].set_xlabel("t (dimensionless)")
    ax[1].set_ylabel("relative error")
    ax[1].set_yscale("log")
    ax[1].legend(fontsize=7)
    ax[1].set_title("Rollout error growth")
    p2 = out / "rollout_heatflux_error.png"
    fig.savefig(p2, dpi=150)
    plt.close(fig)
    return [p1, p2]


def save_report(info: dict, summaries: list[dict], path: str | Path) -> None:
    Path(path).write_text(json.dumps({"training": {k: v for k, v in info.items() if k != "loss"},
                                      "final_loss": info["loss"][-1], "evaluation": summaries}, indent=2))
