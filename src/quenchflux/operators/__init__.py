"""Neural operators for boiling flow (optional; requires ``pip install "quenchflux[ml]"``).

* :mod:`quenchflux.operators.fno`      Fourier neural operator (PyTorch)
* :mod:`quenchflux.operators.bubbleml` range-read access to BubbleML 2.0 simulations
* :mod:`quenchflux.operators.train`    autoregressive training, rollout, wall-heat-flux evaluation

Submodules are imported on demand so the core package works without PyTorch.
"""

__all__ = ["fno", "bubbleml", "train"]
