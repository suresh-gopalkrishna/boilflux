# Contributing to QuenchFlux

Thank you for helping build open tools for inverse boiling heat transfer. Contributions of code,
benchmark cases, experimental datasets, documentation and bug reports are all welcome.

## Ways to contribute

- **Report a bug or ask a question**: open an issue using the templates.
- **Add a method** (e.g. Beck's sequential function specification, Kalman filter, neural-operator
  priors, Bayesian UQ): open an issue first to agree on the interface, then a pull request.
- **Contribute a dataset**: experimental IR quench records are the most valuable contribution to this
  project. See `docs/user-guide/data.md` for the format and metadata we ask for. Data are published
  under CC BY 4.0 with a Zenodo DOI and full credit to the contributors.
- **Improve the documentation**: typos, clarifications and tutorials are very welcome.

## Development setup

```bash
git clone https://github.com/quenchflux/quenchflux.git
cd quenchflux
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev,docs]"
pytest              # run the test suite
mkdocs serve        # preview the website at http://127.0.0.1:8000
```

## Pull-request checklist

- Tests pass (`pytest`) and new code has tests. Numerical code needs a verification test
  (analytic solution, conservation law, or adjoint/finite-difference check).
- `ruff check src tests` reports no errors.
- Public functions have docstrings stating units.
- The documentation and `CHANGELOG.md` are updated when behaviour changes.

## Scientific conduct

Benchmark results must be reproducible from the repository (fixed seeds, recorded options).
Never tune a method on the test cases of the benchmark it is scored on. Report negative results
honestly; they are useful.

## Licensing of contributions

By contributing code you agree that it is released under the project's MIT licence. Data and
documentation contributions are released under CC BY 4.0. Only contribute material you have the
right to share (check your employer's and institution's IP policy for code and experimental data).

## Code of conduct

Participation is governed by our [Code of Conduct](CODE_OF_CONDUCT.md).
