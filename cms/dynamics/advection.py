import numpy as np
from typing import Tuple
from cms.core.grid import Grid

# Try to import Numba kernels
try:
    from cms.dynamics.advection_numba import (
        weno5_reconstruct_x, 
        weno5_reconstruct_y, 
        weno5_reconstruct_z, 
        compute_advection_flux
    )
    HAS_NUMBA = True
except ImportError:
    HAS_NUMBA = False

class WENO5:
    """
    5th-order Weighted Essentially Non-Oscillatory (WENO) advection scheme.
    Ref: IMPLEMENTATION_GUIDE.md Section 6.1
    """
    def __init__(self, grid: Grid):
        self.grid = grid
        self.eps = 1e-6  # Avoid division by zero
        self.use_numba = HAS_NUMBA

    def _reconstruct_weno5(self, f: np.ndarray, axis: int) -> np.ndarray:
        """
        Performs WENO5 reconstruction on a field along a specific axis.
        Returns reconstructed values at interfaces.
        """
        # Shifted versions for stencils
        # f[i-2], f[i-1], f[i], f[i+1], f[i+2]
        f_m2 = np.roll(f, 2, axis=axis)
        f_m1 = np.roll(f, 1, axis=axis)
        f_0  = f
        f_p1 = np.roll(f, -1, axis=axis)
        f_p2 = np.roll(f, -2, axis=axis)

        # Smoothness indicators (Beta)
        b0 = (13/12) * (f_m2 - 2*f_m1 + f_0)**2 + (1/4) * (f_m2 - 4*f_m1 + 3*f_0)**2
        b1 = (13/12) * (f_m1 - 2*f_0 + f_p1)**2 + (1/4) * (f_m1 - f_p1)**2
        b2 = (13/12) * (f_0 - 2*f_p1 + f_p2)**2 + (1/4) * (3*f_0 - 4*f_p1 + f_p2)**2

        # Ideal weights
        d0, d1, d2 = 0.1, 0.6, 0.3

        # Alpha values
        a0 = d0 / (self.eps + b0)**2
        a1 = d1 / (self.eps + b1)**2
        a2 = d2 / (self.eps + b2)**2

        # Weights
        w_sum = a0 + a1 + a2
        w0, w1, w2 = a0/w_sum, a1/w_sum, a2/w_sum

        # Polynomials
        q0 = (1/3)*f_m2 - (7/6)*f_m1 + (11/6)*f_0
        q1 = -(1/6)*f_m1 + (5/6)*f_0 + (1/3)*f_p1
        q2 = (1/3)*f_0 + (5/6)*f_p1 - (1/6)*f_p2

        f_reconstructed = w0*q0 + w1*q1 + w2*q2
        return f_reconstructed

    def advect(self, q: np.ndarray, u: np.ndarray, v: np.ndarray, w: np.ndarray) -> np.ndarray:
        """
        Computes the advection tendency -div(q * U).
        """
        dq_dt = np.zeros_like(q)
        nx, ny, nz = q.shape
        
        if self.use_numba:
            # Buffer for reconstruction
            q_recons = np.empty_like(q)
            
            # X-axis
            weno5_reconstruct_x(q, q_recons, nx, ny, nz)
            compute_advection_flux(q_recons, u, dq_dt, 0, self.grid.dx, nx, ny, nz)
            
            # Y-axis
            weno5_reconstruct_y(q, q_recons, nx, ny, nz)
            compute_advection_flux(q_recons, v, dq_dt, 1, self.grid.dy, nx, ny, nz)
            
            # Z-axis
            weno5_reconstruct_z(q, q_recons, nx, ny, nz)
            compute_advection_flux(q_recons, w, dq_dt, 2, self.grid.dz, nx, ny, nz)
            
            return dq_dt

        # Fallback to pure NumPy implementation (slow)
        for axis, vel in enumerate([u, v, w]):
            # 1. Reconstruct interface values
            # q_left is used for positive velocity, q_right for negative (simplistic upwinding)
            q_star = self._reconstruct_weno5(q, axis)
            
            # 2. Compute flux derivative (assuming uniform grid)
            dx = [self.grid.dx, self.grid.dy, self.grid.dz][axis]
            
            # Simple upwind flux
            flux = vel * q_star
            
            # Grad flux (central difference of the reconstructed interface fluxes)
            # This is a placeholder for a more rigorous flux divergence
            f_p = flux
            f_m = np.roll(flux, 1, axis=axis)
            dq_dt -= (f_p - f_m) / dx
            
        return dq_dt
