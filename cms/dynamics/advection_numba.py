import numpy as np
from numba import jit

@jit(nopython=True, parallel=False, fastmath=False)
def weno5_reconstruct_x(q, out, nx, ny, nz):
    """
    WENO5 reconstruction along X-axis with adaptive epsilon.
    Ref: Borges et al. (2008), Henrick et al. (2005)
    """
    # 13/12 and 1/4 constants
    c1 = 13.0/12.0
    c2 = 0.25
    
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                # Periodic boundaries for indices
                im2 = (i - 2) % nx
                im1 = (i - 1) % nx
                i0  = i
                ip1 = (i + 1) % nx
                ip2 = (i + 2) % nx
                
                v_m2 = q[im2, j, k]
                v_m1 = q[im1, j, k]
                v_0  = q[i0,  j, k]
                v_p1 = q[ip1, j, k]
                v_p2 = q[ip2, j, k]
                
                # Smoothness indicators
                term1 = v_m2 - 2*v_m1 + v_0
                term2 = v_m2 - 4*v_m1 + 3*v_0
                b0 = c1 * term1*term1 + c2 * term2*term2
                
                term1 = v_m1 - 2*v_0 + v_p1
                term2 = v_m1 - v_p1
                b1 = c1 * term1*term1 + c2 * term2*term2
                
                term1 = v_0 - 2*v_p1 + v_p2
                term2 = 3*v_0 - 4*v_p1 + v_p2
                b2 = c1 * term1*term1 + c2 * term2*term2

                # Adaptive epsilon (Henrick et al., 2005)
                local_h = max(abs(v_m2), abs(v_m1), abs(v_0), abs(v_p1), abs(v_p2))
                eps_local = 1e-6 * (local_h**2 + 1e-40)

                # Weights with stability constraints (Borges et al., 2008)
                a0 = 0.1 / max((eps_local + b0)**2, 1e-40)
                a1 = 0.6 / max((eps_local + b1)**2, 1e-40)
                a2 = 0.3 / max((eps_local + b2)**2, 1e-40)

                # Normalize weights to prevent overflow
                a_max = max(a0, a1, a2)
                if a_max > 1e10:
                    a0 /= a_max
                    a1 /= a_max
                    a2 /= a_max

                w_sum = a0 + a1 + a2 + 1e-40
                w0 = a0 / w_sum
                w1 = a1 / w_sum
                w2 = a2 / w_sum
                
                # Candidate polynomials
                q0 = 0.3333333333333333*v_m2 - 1.1666666666666667*v_m1 + 1.8333333333333333*v_0
                q1 = -0.16666666666666666*v_m1 + 0.8333333333333334*v_0 + 0.3333333333333333*v_p1
                q2 = 0.3333333333333333*v_0 + 0.8333333333333334*v_p1 - 0.16666666666666666*v_p2
                
                out[i, j, k] = w0*q0 + w1*q1 + w2*q2

@jit(nopython=True, parallel=False, fastmath=False)
def weno5_reconstruct_y(q, out, nx, ny, nz):
    """WENO5 reconstruction along Y-axis with adaptive epsilon."""
    c1 = 13.0/12.0
    c2 = 0.25
    
    for k in range(nz):
        for i in range(nx):
            for j in range(ny):
                jm2 = (j - 2) % ny
                jm1 = (j - 1) % ny
                j0  = j
                jp1 = (j + 1) % ny
                jp2 = (j + 2) % ny
                
                v_m2 = q[i, jm2, k]
                v_m1 = q[i, jm1, k]
                v_0  = q[i, j0,  k]
                v_p1 = q[i, jp1, k]
                v_p2 = q[i, jp2, k]
                
                term1 = v_m2 - 2*v_m1 + v_0
                term2 = v_m2 - 4*v_m1 + 3*v_0
                b0 = c1 * term1*term1 + c2 * term2*term2
                
                term1 = v_m1 - 2*v_0 + v_p1
                term2 = v_m1 - v_p1
                b1 = c1 * term1*term1 + c2 * term2*term2
                
                term1 = v_0 - 2*v_p1 + v_p2
                term2 = 3*v_0 - 4*v_p1 + v_p2
                b2 = c1 * term1*term1 + c2 * term2*term2

                local_h = max(abs(v_m2), abs(v_m1), abs(v_0), abs(v_p1), abs(v_p2))
                eps_local = 1e-6 * (local_h**2 + 1e-40)

                a0 = 0.1 / max((eps_local + b0)**2, 1e-40)
                a1 = 0.6 / max((eps_local + b1)**2, 1e-40)
                a2 = 0.3 / max((eps_local + b2)**2, 1e-40)

                a_max = max(a0, a1, a2)
                if a_max > 1e10:
                    a0 /= a_max
                    a1 /= a_max
                    a2 /= a_max
                
                w_sum = a0 + a1 + a2 + 1e-40
                w0 = a0 / w_sum
                w1 = a1 / w_sum
                w2 = a2 / w_sum
                
                q0 = 0.3333333333333333*v_m2 - 1.1666666666666667*v_m1 + 1.8333333333333333*v_0
                q1 = -0.16666666666666666*v_m1 + 0.8333333333333334*v_0 + 0.3333333333333333*v_p1
                q2 = 0.3333333333333333*v_0 + 0.8333333333333334*v_p1 - 0.16666666666666666*v_p2
                
                out[i, j, k] = w0*q0 + w1*q1 + w2*q2

@jit(nopython=True, parallel=False, fastmath=False)
def weno5_reconstruct_z(q, out, nx, ny, nz):
    """WENO5 reconstruction along Z-axis with adaptive epsilon."""
    c1 = 13.0/12.0
    c2 = 0.25
    
    for j in range(ny):
        for i in range(nx):
            for k in range(nz):
                km2 = (k - 2) % nz
                km1 = (k - 1) % nz
                k0  = k
                kp1 = (k + 1) % nz
                kp2 = (k + 2) % nz
                
                v_m2 = q[i, j, km2]
                v_m1 = q[i, j, km1]
                v_0  = q[i, j, k0]
                v_p1 = q[i, j, kp1]
                v_p2 = q[i, j, kp2]
                
                term1 = v_m2 - 2*v_m1 + v_0
                term2 = v_m2 - 4*v_m1 + 3*v_0
                b0 = c1 * term1*term1 + c2 * term2*term2
                
                term1 = v_m1 - 2*v_0 + v_p1
                term2 = v_m1 - v_p1
                b1 = c1 * term1*term1 + c2 * term2*term2
                
                term1 = v_0 - 2*v_p1 + v_p2
                term2 = 3*v_0 - 4*v_p1 + v_p2
                b2 = c1 * term1*term1 + c2 * term2*term2

                local_h = max(abs(v_m2), abs(v_m1), abs(v_0), abs(v_p1), abs(v_p2))
                eps_local = 1e-6 * (local_h**2 + 1e-40)

                a0 = 0.1 / max((eps_local + b0)**2, 1e-40)
                a1 = 0.6 / max((eps_local + b1)**2, 1e-40)
                a2 = 0.3 / max((eps_local + b2)**2, 1e-40)

                a_max = max(a0, a1, a2)
                if a_max > 1e10:
                    a0 /= a_max
                    a1 /= a_max
                    a2 /= a_max
                
                w_sum = a0 + a1 + a2 + 1e-40
                w0 = a0 / w_sum
                w1 = a1 / w_sum
                w2 = a2 / w_sum
                
                q0 = 0.3333333333333333*v_m2 - 1.1666666666666667*v_m1 + 1.8333333333333333*v_0
                q1 = -0.16666666666666666*v_m1 + 0.8333333333333334*v_0 + 0.3333333333333333*v_p1
                q2 = 0.3333333333333333*v_0 + 0.8333333333333334*v_p1 - 0.16666666666666666*v_p2
                
                out[i, j, k] = w0*q0 + w1*q1 + w2*q2

@jit(nopython=True, parallel=False, fastmath=True)
def compute_advection_flux(q_recons, vel, out, axis, dx, nx, ny, nz):
    """
    Computes flux gradients -d(u*q)/dx.
    """
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                
                if axis == 0:
                    curr = i
                    prev = (i - 1) % nx
                    
                    flux_curr = vel[i, j, k] * q_recons[i, j, k]
                    flux_prev = vel[prev, j, k] * q_recons[prev, j, k]
                    
                    out[i, j, k] -= (flux_curr - flux_prev) / dx
                    
                elif axis == 1:
                    curr = j
                    prev = (j - 1) % ny
                    
                    flux_curr = vel[i, j, k] * q_recons[i, j, k]
                    flux_prev = vel[i, prev, k] * q_recons[i, prev, k]
                    
                    out[i, j, k] -= (flux_curr - flux_prev) / dx
                    
                elif axis == 2:
                    curr = k
                    prev = (k - 1) % nz
                    
                    flux_curr = vel[i, j, k] * q_recons[i, j, k]
                    flux_prev = vel[i, j, prev] * q_recons[i, j, prev]
                    
                    out[i, j, k] -= (flux_curr - flux_prev) / dx
