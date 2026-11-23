"""Access to the BubbleML 2.0 boiling-simulation dataset (Flash-X, 2D).

BubbleML 2.0 (Hassan et al., NeurIPS 2023 / arXiv 2507.21244) provides
interface-resolved two-phase boiling simulations: temperature, velocity,
pressure and the level-set distance function ``dfun`` (positive in vapour).
The data are hosted on Hugging Face (``hpcforge/BubbleML_2``) under the
**CC BY-NC-SA 4.0** licence: they, and models trained on them, may not be used
commercially. BoilFlux therefore never redistributes the data or trained
weights; :func:`fetch_simulation` downloads only the frames you request
(HTTP range reads), directly from the original host.
"""

from __future__ import annotations

import json
from collections.abc import Sequence
from pathlib import Path

import numpy as np

__all__ = ["HF_BASE", "simulation_url", "fetch_simulation", "load_fetched"]

HF_BASE = "https://huggingface.co/datasets/hpcforge/BubbleML_2/resolve/main"
DEFAULT_FIELDS = ("temperature", "velx", "vely", "dfun")


def simulation_url(subset: str, wall_temp: int, ext: str = "hdf5") -> str:
    """URL of one simulation, e.g. ``simulation_url("SingleBubble-Saturated-FC72-2D", 90)``."""
    return f"{HF_BASE}/{subset}/Twall_{wall_temp}.{ext}"


def _block_mean(a: np.ndarray, f: int) -> np.ndarray:
    if f == 1:
        return a
    t, ny, nx = a.shape
    ny2, nx2 = ny // f * f, nx // f * f
    return a[:, :ny2, :nx2].reshape(t, ny2 // f, f, nx2 // f, f).mean(axis=(2, 4))


def fetch_simulation(
    subset: str,
    wall_temp: int,
    out: str | Path,
    frames: tuple[int, int] = (600, 1200),
    stride: int = 2,
    downsample: int = 3,
    fields: Sequence[str] = DEFAULT_FIELDS,
) -> Path:
    """Download a window of one BubbleML simulation and save it as ``.npz``.

    Parameters
    ----------
    subset, wall_temp:
        Simulation folder and heater temperature, e.g. ``"SingleBubble-Saturated-FC72-2D", 90``.
    frames, stride:
        Frame range ``[start, stop)`` and temporal stride.
    downsample:
        Spatial block-averaging factor (the level set is averaged too, which
        keeps its zero contour to within one coarse cell).
    """
    import urllib.request

    import fsspec  # optional dependency
    import h5py

    url = simulation_url(subset, wall_temp)
    with urllib.request.urlopen(simulation_url(subset, wall_temp, "json")) as r:  # noqa: S310
        meta = json.loads(r.read().decode())
    # the HDF5 datasets are stored contiguously: read one hyperslab per field
    # (large sequential range requests), then subsample in memory
    fh = fsspec.open(url, mode="rb", block_size=16 * 2**20, cache_type="readahead").open()
    with h5py.File(fh, "r") as h:
        data = {}
        stop = min(frames[1], h[fields[0]].shape[0])
        idx = np.arange(frames[0], stop, stride)
        for name in fields:
            arr = h[name][frames[0]:stop][::stride].astype(np.float32)
            data[name] = _block_mean(arr, downsample)
        x = h["x_centers"][()]
        y = h["y_centers"][()]
    x = x[: x.size // downsample * downsample].reshape(-1, downsample).mean(1)
    y = y[: y.size // downsample * downsample].reshape(-1, downsample).mean(1)
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out, x=x, y=y, t=idx * float(meta.get("dt", 1.0)), meta=json.dumps(meta),
                        wall_temp=float(wall_temp), source=url, **data)
    return out


def load_fetched(path: str | Path) -> dict:
    """Load a file written by :func:`fetch_simulation` into a dict of arrays."""
    with np.load(path, allow_pickle=False) as f:
        d = {k: f[k] for k in f.files}
    d["meta"] = json.loads(str(d["meta"]))
    return d
