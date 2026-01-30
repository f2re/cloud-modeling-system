import numpy as np
from typing import Tuple
from cms.config import GridConfig

class Grid:
    """
    3D Eulerian Grid for CMS.
    Ref: IMPLEMENTATION_GUIDE.md Section 1.2 & 7.2
    """
    def __init__(self, config: GridConfig, ghost_cells: int = 3):
        self.nx = config.nx
        self.ny = config.ny
        self.nz = config.nz
        self.dx = config.dx
        self.dy = config.dy
        self.dz = config.dz
        self.ng = ghost_cells

        # Total dimensions including ghost cells
        self.gnx = self.nx + 2 * self.ng
        self.gny = self.ny + 2 * self.ng
        self.gnz = self.nz + 2 * self.ng

        # Coordinate arrays (inner domain)
        self.x = np.linspace(0, (self.nx - 1) * self.dx, self.nx)
        self.y = np.linspace(0, (self.ny - 1) * self.dy, self.ny)
        self.z = np.linspace(0, (self.nz - 1) * self.dz, self.nz)

        # Indexing slices for the inner domain
        self.inner = (
            slice(self.ng, -self.ng),
            slice(self.ng, -self.ng),
            slice(self.ng, -self.ng)
        )

    def create_field(self) -> np.ndarray:
        """Creates a zero-initialized 3D field with ghost cells."""
        return np.zeros((self.gnx, self.gny, self.gnz), dtype=np.float64)

    def get_inner(self, field: np.ndarray) -> np.ndarray:
        """Returns the inner domain of a field, excluding ghost cells."""
        return field[self.inner]

    @property
    def shape(self) -> Tuple[int, int, int]:
        """Returns the full shape of fields including ghost cells."""
        return (self.gnx, self.gny, self.gnz)
