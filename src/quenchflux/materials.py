"""Temperature-dependent thermophysical properties.

Properties are stored as tables and evaluated by piecewise-linear
interpolation (constant extrapolation outside the table range).

The built-in tables for AA6082 and pure nickel are representative literature
values intended for synthetic studies and examples. For quantitative work on
experimental data, replace them with properties measured for your own
specimens (see :meth:`Material.from_csv`).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

__all__ = ["MATERIALS", "Material", "get_material"]


@dataclass
class Material:
    """A solid with temperature-dependent conductivity, density and heat capacity.

    Parameters
    ----------
    name:
        Human-readable material name.
    T_table:
        Temperatures in degC (strictly increasing).
    k_table, rho_table, cp_table:
        Thermal conductivity [W/(m K)], density [kg/m^3] and specific heat
        [J/(kg K)] at ``T_table``.
    """

    name: str
    T_table: np.ndarray
    k_table: np.ndarray
    rho_table: np.ndarray
    cp_table: np.ndarray
    source: str = ""
    meta: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.T_table = np.asarray(self.T_table, dtype=float)
        self.k_table = np.asarray(self.k_table, dtype=float)
        self.rho_table = np.asarray(self.rho_table, dtype=float)
        self.cp_table = np.asarray(self.cp_table, dtype=float)
        n = self.T_table.size
        for arr, nm in ((self.k_table, "k"), (self.rho_table, "rho"), (self.cp_table, "cp")):
            if arr.size != n:
                raise ValueError(f"{nm}_table must have the same length as T_table")
        if n > 1 and np.any(np.diff(self.T_table) <= 0):
            raise ValueError("T_table must be strictly increasing")

    # --- evaluation -----------------------------------------------------
    def k(self, T: np.ndarray) -> np.ndarray:
        """Thermal conductivity [W/(m K)]."""
        return np.interp(T, self.T_table, self.k_table)

    def rho(self, T: np.ndarray) -> np.ndarray:
        """Density [kg/m^3]."""
        return np.interp(T, self.T_table, self.rho_table)

    def cp(self, T: np.ndarray) -> np.ndarray:
        """Specific heat capacity [J/(kg K)]."""
        return np.interp(T, self.T_table, self.cp_table)

    def rho_cp(self, T: np.ndarray) -> np.ndarray:
        """Volumetric heat capacity [J/(m^3 K)]."""
        return self.rho(T) * self.cp(T)

    def diffusivity(self, T: np.ndarray) -> np.ndarray:
        """Thermal diffusivity [m^2/s]."""
        return self.k(T) / self.rho_cp(T)

    # --- constructors ---------------------------------------------------
    @classmethod
    def constant(cls, k: float, rho: float, cp: float, name: str = "constant") -> Material:
        """Material with temperature-independent properties."""
        return cls(name, [0.0], [k], [rho], [cp], source="user (constant)")

    @classmethod
    def from_csv(cls, path: str | Path, name: str | None = None) -> Material:
        """Load a table with columns ``T_C, k, rho, cp`` (header line required)."""
        data = np.genfromtxt(path, delimiter=",", names=True)
        cols = data.dtype.names
        T = data[cols[0]]
        return cls(name or Path(path).stem, T, data[cols[1]], data[cols[2]], data[cols[3]],
                   source=str(path))

    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "T_table": self.T_table.tolist(),
            "k_table": self.k_table.tolist(),
            "rho_table": self.rho_table.tolist(),
            "cp_table": self.cp_table.tolist(),
            "source": self.source,
        }


# Representative property tables (degC). Values compiled from standard
# handbook data; see docs/user-guide/materials.md. Replace for your specimens.
_AA6082 = Material(
    "AA6082",
    T_table=[20, 100, 200, 300, 400, 500, 600],
    k_table=[172, 177, 182, 186, 189, 191, 192],
    rho_table=[2700, 2690, 2670, 2650, 2630, 2610, 2590],
    cp_table=[894, 935, 985, 1030, 1075, 1120, 1170],
    source="representative literature values (aluminium alloy 6082-T6)",
)

_NICKEL = Material(
    "Nickel",
    T_table=[20, 100, 200, 300, 358, 400, 500, 600, 700, 800, 900],
    k_table=[90.7, 82.6, 73.3, 64.0, 60.0, 59.5, 61.5, 64.0, 66.5, 69.0, 71.8],
    rho_table=[8900, 8880, 8850, 8820, 8800, 8790, 8760, 8730, 8700, 8670, 8640],
    cp_table=[445, 470, 505, 550, 600, 535, 530, 535, 545, 555, 565],
    source="representative literature values (Ni 200 / Ni 99.5); Curie peak near 358 degC",
)

_STEEL_C45 = Material(
    "C45",
    T_table=[20, 200, 400, 600, 700, 800, 900],
    k_table=[49.8, 46.0, 40.0, 33.5, 30.0, 25.0, 26.0],
    rho_table=[7850, 7790, 7720, 7650, 7620, 7600, 7570],
    cp_table=[460, 520, 600, 730, 850, 640, 610],
    source="representative literature values (plain carbon steel C45)",
)

MATERIALS: dict[str, Material] = {
    "aa6082": _AA6082,
    "nickel": _NICKEL,
    "ni": _NICKEL,
    "c45": _STEEL_C45,
}


def get_material(name: str) -> Material:
    """Return a built-in material by (case-insensitive) name."""
    key = name.lower().replace("-", "").replace(" ", "")
    if key not in MATERIALS:
        raise KeyError(f"Unknown material '{name}'. Available: {sorted(set(MATERIALS))}")
    return MATERIALS[key]
