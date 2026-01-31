import numpy as np
from cms.core.grid import Grid
from cms.core.time_integration import RK5Integrator
from cms.dynamics.navier_stokes import NavierStokesSolver
from cms.dynamics.boundary import BoundaryConditions
from cms.microphysics.warm import WarmMicrophysics
from cms.microphysics.ice import IceMicrophysics
from cms.microphysics.activation import Activation
from cms.microphysics.condensation import CondensationEvaporation
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
        self.time = 0.0
        
        # Modules
        self.dynamics = NavierStokesSolver(self.grid, self.physics, use_gpu=self.compute_config.use_gpu)
        self.warm_micro = WarmMicrophysics(self.physics)
        self.ice_micro = IceMicrophysics(self.physics)
        self.activation = Activation(self.physics)
        self.condensation = CondensationEvaporation(self.physics)
        self.turbulence = Turbulence(self.grid, self.physics)
        self.bc = BoundaryConditions(self.grid)

        # 1. Dynamics Fields
        self.u = self.grid.create_field()
        self.v = self.grid.create_field()
        self.w = self.grid.create_field()
        self.rho = self.grid.create_field() + self.physics.rho_a_stp
        self.theta = self.grid.create_field() + 300.0
        self.p = self.grid.create_field() + 101325.0
        
        # 2. Water Fields (Vapor, Cloud, Rain)
        self.qv = self.grid.create_field() + 0.01 # Initial humidity 10 g/kg
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
            # Unpack state (19 variables)
            (u, v, w, rho, theta, qv,
             qc, qr, nc, nr, 
             qi, qs, qg, ni, ns, ng, 
             c_agi, n_ccn, n_inp_nat) = state_vector
            
            # Enforce physical constraints on the input state from the integrator
            np.maximum(rho, 1e-9, out=rho)
            for field in [qv, qc, qr, nc, nr, qi, qs, qg, ni, ns, ng, c_agi, n_ccn, n_inp_nat]:
                np.maximum(field, 0.0, out=field)
            
            # --- 0. Thermodynamics: Calculate absolute temperature ---
            # T = theta * (p / p0)^(Rd/cp)
            T = theta * np.power(self.p / self.physics.p0, self.physics.rd / self.physics.cp)

            # --- A. Dynamics ---
            du, dv, dw, drho, dtheta = self.dynamics.compute_tendencies(u, v, w, rho, theta, self.p)
            weno = self.dynamics.weno
            
            # Advection (All scalars)
            # Helper list for advection
            scalars = [qv, qc, qr, nc, nr, qi, qs, qg, ni, ns, ng, c_agi, n_ccn, n_inp_nat]
            d_scalars = [weno.advect(s, u, v, w) for s in scalars]
            
            # Unpack tendencies
            (dqv, dqc, dqr, dnc, dnr, 
             dqi, dqs, dqg, dni, dns, dng, 
             dc_agi, dn_ccn, dn_inp_nat) = d_scalars
            
            # --- B. Turbulence ---
            nu_t = self.turbulence.compute_eddy_viscosity(u, v, w)
            
            # Diffusion (Momentum)
            du += self.turbulence.compute_diffusion(u, nu_t)
            dv += self.turbulence.compute_diffusion(v, nu_t)
            dw += self.turbulence.compute_diffusion(w, nu_t)
            
            # Diffusion (Scalars)
            dqv += self.turbulence.compute_diffusion(qv, nu_t)
            dqc += self.turbulence.compute_diffusion(qc, nu_t)
            dqr += self.turbulence.compute_diffusion(qr, nu_t)
            dnc += self.turbulence.compute_diffusion(nc, nu_t)
            dnr += self.turbulence.compute_diffusion(nr, nu_t)
            dqi += self.turbulence.compute_diffusion(qi, nu_t)
            dqs += self.turbulence.compute_diffusion(qs, nu_t)
            dqg += self.turbulence.compute_diffusion(qg, nu_t)
            dni += self.turbulence.compute_diffusion(ni, nu_t)
            dns += self.turbulence.compute_diffusion(ns, nu_t)
            dng += self.turbulence.compute_diffusion(ng, nu_t)
            dc_agi += self.turbulence.compute_diffusion(c_agi, nu_t)
            dn_ccn += self.turbulence.compute_diffusion(n_ccn, nu_t)
            dn_inp_nat += self.turbulence.compute_diffusion(n_inp_nat, nu_t)

            # --- C. Microphysics ---

            # 1. Condensation/Evaporation
            dqc_cond, dqv_cond = self.condensation.compute_condensation_rate(qv, qc, T, self.p, rho)
            dqc += dqc_cond
            dqv += dqv_cond
            
            # 2. Activation (CCN -> Nc)
            dnc_act, dn_ccn_act = self.activation.compute_activation(w, n_ccn, nc)
            dnc += dnc_act
            dn_ccn += dn_ccn_act
            # Initial mass of activated droplets (r=1um) -> minimal q contribution, 
            # usually neglected or handled by condensation (which isn't fully explicit yet)
            # For mass conservation in this prototype:
            dqc += dnc_act * (4/3 * np.pi * (1e-6)**3 * self.physics.rho_w)

            # 3. Warm
            w_dqc, w_dqr, w_dnc, w_dnr = self.warm_micro.compute_rates(qc, qr, nc, nr, rho)
            dqc += w_dqc; dqr += w_dqr; dnc += w_dnc; dnr += w_dnr
            
            # 4. Ice (with competition)
            i_dqc, i_dqr, i_dqi, i_dqs, i_dqg, i_dni, i_dns, i_dng = self.ice_micro.compute_rates(
                T, rho, qc, qr, nc, qi, qs, qg, ni, ns, ng, c_agi, n_inp_nat
            )
            dqc += i_dqc; dqr += i_dqr; dqi += i_dqi; dqs += i_dqs; dqg += i_dqg
            dni += i_dni; dns += i_dns; dng += i_dng
            
            # --- D. Thermodynamics ---
            # Latent heat from condensation/evaporation
            dtheta -= (self.physics.l_v / self.physics.cp) * dqv_cond
            # Latent heat from ice processes
            dtheta += (self.physics.l_f / self.physics.cp) * i_dqi
            
            return [du, dv, dw, drho, dtheta,
                    dqv, dqc, dqr, dnc, dnr, 
                    dqi, dqs, dqg, dni, dns, dng, 
                    dc_agi, dn_ccn, dn_inp_nat]

        # Use internal step wrapper for packing/unpacking
        self._step_internal(dt, full_rhs)
        self.time += dt
        
    def _step_internal(self, dt, rhs_func):
        # State packing (19 vars)
        state = [self.u, self.v, self.w, self.rho, self.theta,
                 self.qv, self.qc, self.qr, self.nc, self.nr,
                 self.qi, self.qs, self.qg, self.ni, self.ns, self.ng,
                 self.c_agi, self.n_ccn, self.n_inp_nat]
                 
        new_state = self.integrator.step(state, dt, rhs_func)
        
        # Unpack and update
        (self.u, self.v, self.w, self.rho, self.theta,
         self.qv, self.qc, self.qr, self.nc, self.nr,
         self.qi, self.qs, self.qg, self.ni, self.ns, self.ng,
         self.c_agi, self.n_ccn, self.n_inp_nat) = new_state
         
        # Enforce positive definiteness
        # Density must be strictly positive
        np.maximum(self.rho, 1e-9, out=self.rho)
        
        for field in [self.qv, self.qc, self.qr, self.nc, self.nr, 
                      self.qi, self.qs, self.qg, self.ni, self.ns, self.ng,
                      self.c_agi, self.n_ccn, self.n_inp_nat]:
            np.maximum(field, 0.0, out=field)
         
        # BCs
        for field in new_state:
            self.bc.apply_lateral_periodic(field)
        self.bc.apply_bottom_noslip(self.u, self.v, self.w)
