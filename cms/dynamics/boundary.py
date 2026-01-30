import numpy as np
from cms.core.grid import Grid

class BoundaryConditions:
    """
    Handles boundary conditions for the CMS.
    Ref: IMPLEMENTATION_GUIDE.md Section 7.2
    """
    def __init__(self, grid: Grid):
        self.grid = grid

    def apply_lateral_periodic(self, field: np.ndarray):
        """Applies periodic BCs to X and Y boundaries."""
        ng = self.grid.ng
        # X direction
        field[:ng, :, :] = field[-2*ng:-ng, :, :]
        field[-ng:, :, :] = field[ng:2*ng, :, :]
        # Y direction
        field[:, :ng, :] = field[:, -2*ng:-ng, :]
        field[:, -ng:, :] = field[:, ng:2*ng, :]

    def apply_bottom_noslip(self, u: np.ndarray, v: np.ndarray, w: np.ndarray):
        """Applies no-slip BC at the bottom (z=0)."""
        ng = self.grid.ng
        # u, v, w = 0 at the boundary
        u[:, :, :ng] = 0
        v[:, :, :ng] = 0
        w[:, :, :ng] = 0

    def apply_sponge_layer(self, field: np.ndarray, field_ref: np.ndarray, dt: float):
        """
        Applies an absorbing sponge layer at the top boundary.
        Ref: Eq 7.2
        """
        nz = self.grid.nz
        ng = self.grid.ng
        z = self.grid.z
        z_top = z[-1]
        z_s = 0.8 * z_top  # Sponge starts at 80% height
        
        # Calculate nu_sponge
        nu0 = 0.1
        for k in range(ng, self.grid.gnz):
            z_curr = (k - ng) * self.grid.dz
            if z_curr > z_s:
                nu = nu0 * np.sin(0.5 * np.pi * (z_curr - z_s) / (z_top - z_s))**2
                field[:, :, k] -= nu * (field[:, :, k] - field_ref[:, :, k]) * dt
