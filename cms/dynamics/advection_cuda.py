from numba import cuda
import math

@cuda.jit
def weno5_reconstruct_x_kernel(q, out, nx, ny, nz):
    """
    CUDA kernel for WENO5 reconstruction along X-axis.
    Mapping: i (x), j (y), k (z)
    """
    i, j, k = cuda.grid(3)
    
    if i < nx and j < ny and k < nz:
        eps = 1e-6
        c1 = 13.0/12.0
        c2 = 0.25
        
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
        
        # Weights
        a0 = 0.1 / ((eps + b0)*(eps + b0))
        a1 = 0.6 / ((eps + b1)*(eps + b1))
        a2 = 0.3 / ((eps + b2)*(eps + b2))
        
        w_sum = a0 + a1 + a2
        w0 = a0 / w_sum
        w1 = a1 / w_sum
        w2 = a2 / w_sum
        
        # Polynomials
        q0 = 0.3333333333333333*v_m2 - 1.1666666666666667*v_m1 + 1.8333333333333333*v_0
        q1 = -0.16666666666666666*v_m1 + 0.8333333333333334*v_0 + 0.3333333333333333*v_p1
        q2 = 0.3333333333333333*v_0 + 0.8333333333333334*v_p1 - 0.16666666666666666*v_p2
        
        out[i, j, k] = w0*q0 + w1*q1 + w2*q2

@cuda.jit
def weno5_reconstruct_y_kernel(q, out, nx, ny, nz):
    """
    CUDA kernel for WENO5 reconstruction along Y-axis.
    """
    i, j, k = cuda.grid(3)
    
    if i < nx and j < ny and k < nz:
        eps = 1e-6
        c1 = 13.0/12.0
        c2 = 0.25
        
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
        
        a0 = 0.1 / ((eps + b0)**2)
        a1 = 0.6 / ((eps + b1)**2)
        a2 = 0.3 / ((eps + b2)**2)
        
        w_sum = a0 + a1 + a2
        w0 = a0 / w_sum
        w1 = a1 / w_sum
        w2 = a2 / w_sum
        
        q0 = 0.3333333333333333*v_m2 - 1.1666666666666667*v_m1 + 1.8333333333333333*v_0
        q1 = -0.16666666666666666*v_m1 + 0.8333333333333334*v_0 + 0.3333333333333333*v_p1
        q2 = 0.3333333333333333*v_0 + 0.8333333333333334*v_p1 - 0.16666666666666666*v_p2
        
        out[i, j, k] = w0*q0 + w1*q1 + w2*q2

@cuda.jit
def weno5_reconstruct_z_kernel(q, out, nx, ny, nz):
    """
    CUDA kernel for WENO5 reconstruction along Z-axis.
    """
    i, j, k = cuda.grid(3)
    
    if i < nx and j < ny and k < nz:
        eps = 1e-6
        c1 = 13.0/12.0
        c2 = 0.25
        
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
        
        a0 = 0.1 / ((eps + b0)**2)
        a1 = 0.6 / ((eps + b1)**2)
        a2 = 0.3 / ((eps + b2)**2)
        
        w_sum = a0 + a1 + a2
        w0 = a0 / w_sum
        w1 = a1 / w_sum
        w2 = a2 / w_sum
        
        q0 = 0.3333333333333333*v_m2 - 1.1666666666666667*v_m1 + 1.8333333333333333*v_0
        q1 = -0.16666666666666666*v_m1 + 0.8333333333333334*v_0 + 0.3333333333333333*v_p1
        q2 = 0.3333333333333333*v_0 + 0.8333333333333334*v_p1 - 0.16666666666666666*v_p2
        
        out[i, j, k] = w0*q0 + w1*q1 + w2*q2

@cuda.jit
def compute_advection_flux_kernel(q_recons, vel, out, axis, dx, nx, ny, nz):
    """
    Computes flux gradients -d(u*q)/dx on GPU.
    """
    i, j, k = cuda.grid(3)
    
    if i < nx and j < ny and k < nz:
        if axis == 0:
            prev = (i - 1) % nx
            flux_curr = vel[i, j, k] * q_recons[i, j, k]
            flux_prev = vel[prev, j, k] * q_recons[prev, j, k]
            
        elif axis == 1:
            prev = (j - 1) % ny
            flux_curr = vel[i, j, k] * q_recons[i, j, k]
            flux_prev = vel[i, prev, k] * q_recons[i, prev, k]
            
        elif axis == 2:
            prev = (k - 1) % nz
            flux_curr = vel[i, j, k] * q_recons[i, j, k]
            flux_prev = vel[i, j, prev] * q_recons[i, j, prev]
            
        # Atomic sub not strictly needed if each thread owns output[i,j,k]
        # But out is shared accumulator, so let's verify.
        # Advect method:
        # out = 0
        # out -= flux_x
        # out -= flux_y
        # ...
        # If we run 3 kernels sequentially updating 'out', no race condition on same element.
        
        out[i, j, k] -= (flux_curr - flux_prev) / dx
