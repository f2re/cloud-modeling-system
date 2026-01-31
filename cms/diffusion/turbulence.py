import numpy as np
from cms.core.grid import Grid
from cms.config import PhysicsConfig

class Turbulence:
    """
    Turbulence and Dispersion module.
    Ref: IMPLEMENTATION_GUIDE.md Section 2.2 & 4.2
    """
    def __init__(self, grid: Grid, config: PhysicsConfig):
        self.grid = grid
        self.config = config
        self.cs_sq = 0.18**2 # Smagorinsky constant squared

    def compute_eddy_viscosity(self, u: np.ndarray, v: np.ndarray, w: np.ndarray) -> np.ndarray:
        """
        Computes subgrid eddy viscosity (nu_t) with stability clips.
        Eq 4.2: nu_t = (Cs * Delta)^2 * |S|
        """
        # Grid spacing filter width
        delta = (self.grid.dx * self.grid.dy * self.grid.dz)**(1/3)
        
        # Strain rate tensor norm with clipping
        max_grad = 1.0  # Max gradient [1/s]
        du_dx = np.clip(np.gradient(u, self.grid.dx, axis=0), -max_grad, max_grad)
        dv_dy = np.clip(np.gradient(v, self.grid.dy, axis=1), -max_grad, max_grad)
        dw_dz = np.clip(np.gradient(w, self.grid.dz, axis=2), -max_grad, max_grad)
        
        s_sq = 2 * (du_dx**2 + dv_dy**2 + dw_dz**2)
        
        nu_t = (self.cs_sq * delta**2) * np.sqrt(s_sq)
        
        # Clip final viscosity to a realistic maximum
        nu_t_max = 1000.0  # [m^2/s]
        nu_t = np.clip(nu_t, 0.0, nu_t_max)
        
        return nu_t

    def compute_diffusion(self, field: np.ndarray, nu_t: np.ndarray) -> np.ndarray:
        """
        Computes diffusion tendency: div(nu_t * grad(field))
        """
        # Fluxes
        grad_x = np.gradient(field, self.grid.dx, axis=0)
        grad_y = np.gradient(field, self.grid.dy, axis=1)
        grad_z = np.gradient(field, self.grid.dz, axis=2)
        
        flux_x = nu_t * grad_x
        flux_y = nu_t * grad_y
        flux_z = nu_t * grad_z
        
        div_flux = (np.gradient(flux_x, self.grid.dx, axis=0) +
                    np.gradient(flux_y, self.grid.dy, axis=1) +
                    np.gradient(flux_z, self.grid.dz, axis=2))
                    
        return div_flux
