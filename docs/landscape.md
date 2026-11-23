# Landscape: software and data for inverse boiling heat transfer

This page lists the tools and datasets BoilFlux relates to, so users can choose the right one and
see where BoilFlux fits. It is based on a literature and repository survey carried out in
November 2026. Corrections and additions are welcome through a pull request.

!!! info "Summary"
    The inverse methods themselves are mature: Beck's function specification, Tikhonov, Alifanov iterative
    regularisation, adjoint CGM, Kalman filters and Bayesian inference. Maintained open software for
    **boiling and quenching** is not. Most quench IHCP solvers are commercial and 1D, restricted, in-house,
    or described in papers without a public repository. No public inverse benchmark with
    ground truth existed for boiling or quenching.

## Open or partially open IHCP codes

| Tool | Scope | Boiling-specific | Status |
|---|---|---|---|
| Bohacek et al., sequential IHCP in OpenFOAM (OpenFOAM Journal, 2021) | 2D nonlinear, sequential + BOBYQA | applied to spray cooling | source code as a journal supplement |
| Bohacek et al., generic-surface 3D IHCP in OpenFOAM (IJTS, 2026) | up to 3D, polyhedral meshes, MPI | applied to spray/jet coolers | described as open source; no public repository found |
| ihtFoam, Nabbout & Sommerfeld (arXiv 2403.07490) | 3D adjoint CGM in OpenFOAM | jet quenching with IR | no public repository found |
| [IHCP_FEniCS](https://github.com/mberg88/IHCP_FEniCS) | 2D adjoint CGM, FEniCS 2016 | no | dormant since 2017 |
| [ROMSOC inverse heat transfer benchmark](https://github.com/ROMSOC/benchmarks-inverse-heat-transfer) | 3D continuous-casting mould, OpenFOAM/ITHACA-FV | no | v1.0.0 on Zenodo (2021) |
| [HFITS](https://github.com/ulfsri/HFITS) | IR thermography on thin planar surfaces (SoftwareX 2024) | no (fire research) | public |
| HYPERION (Nuclear Fusion, 2024) | 2D IR inversion for tokamak divertors | no | described as open source |
| PIKAN boiling, Khodakarami et al. (arXiv 2610.20613) | 3D PIKAN/PINN, IR pool boiling | yes (pool boiling) | code promised on publication |

## Commercial and restricted tools

| Tool | Scope | Availability |
|---|---|---|
| IHCP1D (Beck Engineering Consultants) | 1D nonlinear function specification | commercial |
| INTEMP | 1D/2D nonlinear quench inverse | commercial |
| NASA IHC Solve (LAR-18924-1) | 1D, two thermocouples, MATLAB | U.S. government purpose release only |
| ivf SmartQuench / SQintegra | ISO 9950 probe, HTC inversion | commercial instrument and software |
| CHTE probe software, Liščić–Nanmac probe | quench probes | consortium / commercial |
| DEFORM-HT, Abaqus, ANSYS, COMSOL | inverse via optimisation loops | general-purpose, licensed |

## General inverse-problem frameworks (building blocks)

[hIPPYlib](https://hippylib.github.io) (adjoint-based deterministic and Bayesian inversion on FEniCS),
CUQIpy (computational UQ for inverse problems), FEniCSx with
dolfin-adjoint/pyadjoint, pyGIMLi, and SciML libraries (DeepXDE, NVIDIA PhysicsNeMo, neuraloperator).
They are powerful, but none of them includes boiling curves, wetting fronts, IR-camera pipelines or quench
benchmarks. BoilFlux plans to interoperate with them (for example for UQ) rather than replace them.

## Datasets

| Dataset | Content | Inverse-ready |
|---|---|---|
| ROMSOC benchmark | continuous-casting mould, thermocouple data with ground truth | yes (not boiling) |
| Szénási & Felde (Data, 2019) | synthetic 1D quench HTC curves with cooling histories | yes (synthetic 1D) |
| [BubbleML](https://github.com/HPCForge/BubbleML), BubbleML 2.0 | simulations of pool, flow and subcooled boiling | forward; can be adapted |
| Purdue HFE-7100 IR pool-boiling data | IR + phase maps, ITO-sapphire heater | listed as in preparation |
| **Public IR quench data on hot metals** | — | **none found; BF-Bench tier 4 is planned to fill this gap** |

## Where BoilFlux fits

BoilFlux aims to be the **reference implementation and benchmark** for inverse problems in boiling and
quenching. It brings together a nonlinear multi-dimensional adjoint core, moving-plate kinematics, IR data
handling and reproducible ground-truth cases in one tested, documented and permissively licensed
package. Wrappers for the tools above and comparisons against them are part of the [roadmap](roadmap.md).
