import numpy as np
from cms.core.grid import Grid
from cms.dynamics.advection import WENO5
from cms.config import PhysicsConfig

class NavierStokesSolver:
    """
    Compressible Navier-Stokes solver.
    Ref: IMPLEMENTATION_GUIDE.md Section 4.1
    """
    def __init__(self, grid: Grid, config: PhysicsConfig):
        self.grid = grid
        self.config = config
        self.weno = WENO5(grid)

    def compute_tendencies(self, 
                           u: np.ndarray, v: np.ndarray, w: np.ndarray, 
                           rho: np.ndarray, theta: np.ndarray, 
                           p: np.ndarray):
        """
        Computes tendencies for momentum, density, and potential temperature.
        """
        # 1. Advection tendencies using WENO5
        du_dt = self.weno.advect(u, u, v, w)
        dv_dt = self.weno.advect(v, u, v, w)
        dw_dt = self.weno.advect(w, u, v, w)
        drho_dt = self.weno.advect(rho, u, v, w)
        dtheta_dt = self.weno.advect(theta, u, v, w)

        # 2. Pressure gradient terms
        # grad(p) / rho
        dp_dx = np.gradient(p, self.grid.dx, axis=0)
        dp_dy = np.gradient(p, self.grid.dy, axis=1)
        dp_dz = np.gradient(p, self.grid.dz, axis=2)

        du_dt -= dp_dx / rho
        dv_dt -= dp_dy / rho
        dw_dt -= dp_dz / rho

        # 3. Buoyancy (for vertical momentum w)
        # g * (theta_v - theta_v0) / theta_v0
        # For now using theta as a proxy for theta_v
        theta0 = 300.0  # Placeholder for reference state
        buoyancy = self.config.g * (theta - theta0) / theta0
        dw_dt += buoyancy

        return du_dt, dv_dt, dw_dt, drho_dt, dtheta_dt
