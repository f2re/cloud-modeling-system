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

        # Cache for meshgrid arrays
        self._X = None
        self._Y = None
        self._Z = None

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

    @property
    def gx(self) -> np.ndarray:
        """1D coordinate array for the full grid (x-axis)."""
        return np.linspace(
            -(self.ng * self.dx),
            (self.nx - 1 + self.ng) * self.dx,
            self.gnx
        )

    @property
    def gy(self) -> np.ndarray:
        """1D coordinate array for the full grid (y-axis)."""
        return np.linspace(
            -(self.ng * self.dy),
            (self.ny - 1 + self.ng) * self.dy,
            self.gny
        )

    @property
    def gz(self) -> np.ndarray:
        """1D coordinate array for the full grid (z-axis)."""
        return np.linspace(
            -(self.ng * self.dz),
            (self.nz - 1 + self.ng) * self.dz,
            self.gnz
        )

    def _create_meshgrid(self):
        """Creates and caches the 3D meshgrid."""
        if self._X is None:
            self._X, self._Y, self._Z = np.meshgrid(self.gx, self.gy, self.gz, indexing='ij')

    @property
    def X(self) -> np.ndarray:
        """3D meshgrid for X coordinates (full grid)."""
        self._create_meshgrid()
        return self._X

    @property
    def Y(self) -> np.ndarray:
        """3D meshgrid for Y coordinates (full grid)."""
        self._create_meshgrid()
        return self._Y

    @property
    def Z(self) -> np.ndarray:
        """3D meshgrid for Z coordinates (full grid)."""
        self._create_meshgrid()
        return self._Z
