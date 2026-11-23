"""BoilFlux: open inverse heat transfer for boiling and quenching.

Reconstruct the transient, spatially resolved boiling heat flux on the
quenched face of a (moving) hot plate from infrared thermography of the
opposite face, using a nonlinear 3D finite-volume model and the adjoint
conjugate-gradient method with Tikhonov / discrepancy-principle regularisation.

Reference method: S.B. Gopalkrishna, A.K. Nallathambi, E. Specht,
"Inverse heat conduction method to estimate heat flux during quenching of
hot metals with infrared thermography measurements", Applied Thermal
Engineering 304 (2026) 132481, doi:10.1016/j.applthermaleng.2026.132481.
"""

from .boiling import FFZ_DEFAULT, IZ_DEFAULT, SpechtBoilingCurve
from .dataset import BoilDataset
from .forward import ForwardModel, SurfaceLoss, Trajectory
from .grid import PlateGrid
from .inverse import AdjointCGM, CGMOptions, InverseResult
from .kinematics import PlateMotion, lab_frames_to_plate
from .materials import Material, get_material
from .metrics import flux_metrics
from .noise import estimate_noise_sigma
from .synthetic import JetQuenchScenario

__version__ = "0.1.0"

__all__ = [
    "FFZ_DEFAULT",
    "IZ_DEFAULT",
    "AdjointCGM",
    "CGMOptions",
    "ForwardModel",
    "InverseResult",
    "JetQuenchScenario",
    "Material",
    "PlateGrid",
    "PlateMotion",
    "BoilDataset",
    "SpechtBoilingCurve",
    "SurfaceLoss",
    "Trajectory",
    "__version__",
    "estimate_noise_sigma",
    "flux_metrics",
    "get_material",
    "lab_frames_to_plate",
]
