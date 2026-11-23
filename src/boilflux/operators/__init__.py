"""Neural operators for boiling flow (optional; requires ``pip install "boilflux[ml]"``).

* :mod:`boilflux.operators.fno`      Fourier neural operator (PyTorch)
* :mod:`boilflux.operators.bubbleml` range-read access to BubbleML 2.0 simulations
* :mod:`boilflux.operators.train`    autoregressive training, rollout, wall-heat-flux evaluation

Submodules are imported on demand so the core package works without PyTorch.
"""

__all__ = ["fno", "bubbleml", "train"]
