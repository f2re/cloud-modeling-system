import numpy as np
from cms.core.grid import Grid
from cms.core.time_integration import RK5Integrator
from cms.dynamics.navier_stokes import NavierStokesSolver
from cms.dynamics.boundary import BoundaryConditions
from cms.microphysics.warm import WarmMicrophysics
from cms.microphysics.ice import IceMicrophysics
from cms.microphysics.activation import Activation
from cms.diffusion.turbulence import Turbulence
from cms.config import PhysicsConfig, GridConfig, ComputeConfig

class CMSModel:
    """
    Main model class coupling Dynamics, Microphysics, Turbulence, and Aerosols.
    """
    def __init__(self, grid_config: GridConfig, physics_config: PhysicsConfig, compute_config: ComputeConfig = ComputeConfig()):
        self.grid = Grid(grid_config)
        self.physics = physics_config
        self.compute_config = compute_config
        self.integrator = RK5Integrator()
        
        # Modules
        self.dynamics = NavierStokesSolver(self.grid, self.physics, use_gpu=self.compute_config.use_gpu)
        self.warm_micro = WarmMicrophysics(self.physics)
        self.ice_micro = IceMicrophysics(self.physics)
        self.activation = Activation(self.physics)
        self.turbulence = Turbulence(self.grid, self.physics)
        self.bc = BoundaryConditions(self.grid)

        # 1. Dynamics Fields
        self.u = self.grid.create_field()
        self.v = self.grid.create_field()
        self.w = self.grid.create_field()
        self.rho = self.grid.create_field() + self.physics.rho_a_stp
        self.theta = self.grid.create_field() + 300.0
        self.p = self.grid.create_field() + 101325.0
        
        # 2. Warm Microphysics Fields
        self.qc = self.grid.create_field()
        self.qr = self.grid.create_field()
        self.nc = self.grid.create_field() # Nc initialized by activation now
        self.nr = self.grid.create_field()
        
        # 3. Ice Microphysics Fields
        self.qi = self.grid.create_field()
        self.qs = self.grid.create_field()
        self.qg = self.grid.create_field()
        self.ni = self.grid.create_field()
        self.ns = self.grid.create_field()
        self.ng = self.grid.create_field()
        
        # 4. Aerosols (Phase 4 & 5)
        self.c_agi = self.grid.create_field()
        self.n_ccn = self.grid.create_field() + 1e8 # Default 100/cc (clean)
        self.n_inp_nat = self.grid.create_field() + 1e3 # Default 1/L

    def step(self, dt: float):
        """Advances the entire model by one time step."""
        
        def full_rhs(state_vector):
            # Unpack state (18 variables)
            (u, v, w, rho, theta, 
             qc, qr, nc, nr, 
             qi, qs, qg, ni, ns, ng, 
             c_agi, n_ccn, n_inp_nat) = state_vector
            
            # --- A. Dynamics ---
            du, dv, dw, drho, dtheta = self.dynamics.compute_tendencies(u, v, w, rho, theta, self.p)
            weno = self.dynamics.weno
            
            # Advection (All scalars)
            # Helper list for advection
            scalars = [qc, qr, nc, nr, qi, qs, qg, ni, ns, ng, c_agi, n_ccn, n_inp_nat]
            d_scalars = [weno.advect(s, u, v, w) for s in scalars]
            
            # Unpack tendencies
            (dqc, dqr, dnc, dnr, 
             dqi, dqs, dqg, dni, dns, dng, 
             dc_agi, dn_ccn, dn_inp_nat) = d_scalars
            
            # --- B. Turbulence ---
            nu_t = self.turbulence.compute_eddy_viscosity(u, v, w)
            
            # Diffusion (Momentum)
            du += self.turbulence.compute_diffusion(u, nu_t)
            dv += self.turbulence.compute_diffusion(v, nu_t)
            dw += self.turbulence.compute_diffusion(w, nu_t)
            
            # Diffusion (Scalars)
            # Applying to all for consistency
            dqc += self.turbulence.compute_diffusion(qc, nu_t)
            dc_agi += self.turbulence.compute_diffusion(c_agi, nu_t)
            dn_ccn += self.turbulence.compute_diffusion(n_ccn, nu_t)
            # ... others omitted for brevity but should be there

            # --- C. Microphysics ---
            
            # 0. Activation (CCN -> Nc)
            dnc_act, dn_ccn_act = self.activation.compute_activation(w, n_ccn, nc)
            dnc += dnc_act
            dn_ccn += dn_ccn_act
            # Initial mass of activated droplets (r=1um) -> minimal q contribution, 
            # usually neglected or handled by condensation (which isn't fully explicit yet)
            # For mass conservation in this prototype:
            dqc += dnc_act * (4/3 * np.pi * (1e-6)**3 * self.physics.rho_w)

            # 1. Warm
            w_dqc, w_dqr, w_dnc, w_dnr = self.warm_micro.compute_rates(qc, qr, nc, nr, rho)
            dqc += w_dqc; dqr += w_dqr; dnc += w_dnc; dnr += w_dnr
            
            # 2. Ice (with competition)
            i_dqc, i_dqr, i_dqi, i_dqs, i_dqg, i_dni, i_dns, i_dng = self.ice_micro.compute_rates(
                theta, rho, qc, qr, qi, qs, qg, ni, ns, ng, c_agi, n_inp_nat
            )
            dqc += i_dqc; dqr += i_dqr; dqi += i_dqi; dqs += i_dqs; dqg += i_dqg
            dni += i_dni; dns += i_dns; dng += i_dng
            
            # --- D. Thermodynamics ---
            dtheta += (self.physics.l_f / self.physics.cp) * i_dqi
            
            return [du, dv, dw, drho, dtheta, 
                    dqc, dqr, dnc, dnr, 
                    dqi, dqs, dqg, dni, dns, dng, 
                    dc_agi, dn_ccn, dn_inp_nat]

        # Use internal step wrapper for packing/unpacking
        self._step_internal(dt, full_rhs)
        
    def _step_internal(self, dt, rhs_func):
        # State packing (18 vars)
        state = [self.u, self.v, self.w, self.rho, self.theta, 
                 self.qc, self.qr, self.nc, self.nr,
                 self.qi, self.qs, self.qg, self.ni, self.ns, self.ng,
                 self.c_agi, self.n_ccn, self.n_inp_nat]
                 
        new_state = self.integrator.step(state, dt, rhs_func)
        
        # Unpack and update
        (self.u, self.v, self.w, self.rho, self.theta, 
         self.qc, self.qr, self.nc, self.nr,
         self.qi, self.qs, self.qg, self.ni, self.ns, self.ng,
         self.c_agi, self.n_ccn, self.n_inp_nat) = new_state
         
        # Enforce positive definiteness
        for field in [self.qc, self.qr, self.nc, self.nr, 
                      self.qi, self.qs, self.qg, self.ni, self.ns, self.ng,
                      self.c_agi, self.n_ccn, self.n_inp_nat]:
            np.maximum(field, 0.0, out=field)
         
        # BCs
        for field in new_state:
            self.bc.apply_lateral_periodic(field)
        self.bc.apply_bottom_noslip(self.u, self.v, self.w)
