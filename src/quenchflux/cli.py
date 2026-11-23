"""Command-line interface: ``quenchflux <command>``.

Commands
--------
info        Show version, materials and benchmark cases.
demo        Generate the moving-plate example, invert it and write figures.
bench       Run one or all benchmark cases and print the score table.
synth       Generate a synthetic dataset from a YAML scenario file.
invert      Invert a dataset (.npz) and write the result and figures.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def _load_yaml(path):
    import yaml

    with open(path) as f:
        return yaml.safe_load(f) or {}


def cmd_info(_args) -> int:
    from . import __version__
    from .benchmarks import CASES
    from .materials import MATERIALS

    print(f"QuenchFlux {__version__}")
    print("materials:", ", ".join(sorted({m.name for m in MATERIALS.values()})))
    print("benchmark cases:")
    for k, v in CASES.items():
        print(f"  {k:24s} {v['description']}")
    return 0


def cmd_demo(args) -> int:
    from .benchmarks import run_case

    scores = run_case("M3D-AA6082", outdir=args.out, verbose=not args.quiet)
    print(json.dumps(scores, indent=2))
    print(f"figures and data written to {Path(args.out) / 'M3D-AA6082'}")
    return 0


def cmd_bench(args) -> int:
    from .benchmarks import CASES, run_case

    names = list(CASES) if args.case == "all" else [args.case]
    rows = [run_case(n, outdir=args.out, verbose=not args.quiet, plots=not args.no_plots) for n in names]
    hdr = f"{'case':24s} {'R2':>6s} {'RMSE':>7s} {'MAPE%':>7s} {'MHF%':>7s} {'iters':>5s} {'time s':>7s}"
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        print(f"{r['case']:24s} {r['r2']:6.3f} {r['rmse_MW_m2']:7.3f} {r['mape_percent']:7.2f} "
              f"{r['mhf_error_percent']:7.2f} {r['iterations']:5d} {r['runtime_s']:7.1f}")
    if args.out:
        Path(args.out).mkdir(parents=True, exist_ok=True)
        (Path(args.out) / "scores.json").write_text(json.dumps(rows, indent=2))
    return 0


def cmd_synth(args) -> int:
    from .synthetic import JetQuenchScenario

    cfg = _load_yaml(args.config)
    scen = JetQuenchScenario.from_dict(cfg.get("scenario", cfg))
    data = scen.generate()
    path = data.save(args.output)
    print(f"wrote {path}  (nt={data.nt}, grid={data.grid.nx}x{data.grid.ny}x{data.grid.nz})")
    return 0


def cmd_invert(args) -> int:
    from .dataset import QuenchDataset
    from .forward import SurfaceLoss
    from .inverse import AdjointCGM, CGMOptions
    from .materials import Material
    from .metrics import flux_metrics

    cfg = _load_yaml(args.config) if args.config else {}
    data = QuenchDataset.load(args.dataset)
    opts = CGMOptions(**cfg.get("cgm", {}))
    mat = Material.from_csv(cfg["material_csv"]) if cfg.get("material_csv") else None
    loss = SurfaceLoss(**cfg["top_loss"]) if cfg.get("top_loss") else None
    res = AdjointCGM(data, material=mat, top_loss=loss, options=opts).solve()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    res.save(out / "result.npz")
    if data.q_true is not None:
        sl = res.reliable
        scores = flux_metrics(res.q[sl], data.q_true[sl])
        (out / "scores.json").write_text(json.dumps(scores, indent=2))
        print(json.dumps(scores, indent=2))
    if not args.no_plots:
        from .plotting import save_all

        save_all(data, res, out)
    print(f"results written to {out}")
    return 0


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        prog="quenchflux", description="Inverse heat transfer for boiling and quenching."
    )
    sub = p.add_subparsers(dest="command", required=True)
    sub.add_parser("info", help="version, materials, benchmark cases").set_defaults(func=cmd_info)

    d = sub.add_parser("demo", help="run the 3D moving-plate example")
    d.add_argument("--out", default="quenchflux_demo")
    d.add_argument("--quiet", action="store_true")
    d.set_defaults(func=cmd_demo)

    b = sub.add_parser("bench", help="run benchmark cases")
    b.add_argument("--case", default="all")
    b.add_argument("--out", default="quenchflux_bench")
    b.add_argument("--no-plots", action="store_true")
    b.add_argument("--quiet", action="store_true")
    b.set_defaults(func=cmd_bench)

    s = sub.add_parser("synth", help="generate a synthetic dataset from YAML")
    s.add_argument("config")
    s.add_argument("-o", "--output", default="dataset.npz")
    s.set_defaults(func=cmd_synth)

    i = sub.add_parser("invert", help="invert a dataset")
    i.add_argument("dataset")
    i.add_argument("-c", "--config", default=None)
    i.add_argument("-o", "--out", default="quenchflux_result")
    i.add_argument("--no-plots", action="store_true")
    i.set_defaults(func=cmd_invert)

    args = p.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
