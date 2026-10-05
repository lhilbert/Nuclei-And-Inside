"""The only converter between microns and voxels.

Every scale parameter in this project is stated in microns. This module is the single place
that turns one into a voxel count, and it always needs the grid to do it. If a sigma appears
in voxels anywhere else, that is a bug.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict


@dataclass(frozen=True)
class Grid:
    """Voxel spacing in microns, as (z, y, x), plus the field-frame origin of voxel (0,0,0)."""
    dz_um: float
    dy_um: float
    dx_um: float
    oz_um: float = 0.0
    oy_um: float = 0.0
    ox_um: float = 0.0

    @property
    def spacing(self) -> tuple[float, float, float]:
        return (self.dz_um, self.dy_um, self.dx_um)

    @property
    def origin(self) -> tuple[float, float, float]:
        return (self.oz_um, self.oy_um, self.ox_um)

    @property
    def voxel_volume_um3(self) -> float:
        return self.dz_um * self.dy_um * self.dx_um

    @property
    def is_isotropic(self) -> bool:
        s = self.spacing
        return max(s) / min(s) < 1.01

    def um_to_vox(self, um: float, axis: int) -> float:
        return um / self.spacing[axis]

    def sigma_vox(self, sigma_um: float, sigma_z_um: float | None = None) -> tuple[float, float, float]:
        """A physical sigma as a per-axis voxel sigma. Anisotropic z is explicit, never implied."""
        sz = (sigma_z_um if sigma_z_um is not None else sigma_um) / self.dz_um
        return (sz, sigma_um / self.dy_um, sigma_um / self.dx_um)

    def index_to_um(self, zyx) -> tuple[float, float, float]:
        z, y, x = zyx
        return (self.oz_um + z * self.dz_um,
                self.oy_um + y * self.dy_um,
                self.ox_um + x * self.dx_um)

    def to_dict(self) -> dict:
        return asdict(self)

    @staticmethod
    def from_dict(d: dict) -> "Grid":
        return Grid(**{k: float(v) for k, v in d.items()})


def isotropic(spacing_um: float) -> Grid:
    return Grid(spacing_um, spacing_um, spacing_um)


def sample_um(vol, grid: Grid, pts_um, order: int = 1):
    """Interpolate `vol` at physical (z, y, x) points, both expressed in the same frame.

    Here rather than anywhere else because turning a micron coordinate into a voxel index is
    exactly what this module is for, and doing it inline somewhere else is how a crop's origin
    gets dropped.
    """
    import numpy as np
    from scipy.ndimage import map_coordinates
    p = np.asarray(pts_um, dtype=float)
    idx = np.empty_like(p)
    for a, (o, d) in enumerate(zip(grid.origin, grid.spacing)):
        idx[:, a] = (p[:, a] - o) / d
    return map_coordinates(vol, idx.T, order=order, mode="nearest")
