import numpy as np
from typing import Tuple
from cms.config import PhysicsConfig

class IceMicrophysics:
    """
    Double-moment ice microphysics including AgI interaction.
    Ref: IMPLEMENTATION_GUIDE.md Section 3.5.2 - 3.5.5
    """
    def __init__(self, config: PhysicsConfig):
        self.config = config

    def compute_nucleation_agi(self, T: np.ndarray, N_agi: np.ndarray, N_inp_nat: np.ndarray) -> np.ndarray:
        """
        Primary heterogeneous nucleation by AgI with competition.
        Ref: Eq 3.5.2 & 3.6
        """
        # INF(T)
        T_c = T - 273.15
        inf = np.zeros_like(T)
        mask_nuc = T_c < -5.0
        
        # 0.0007 * exp(0.28 * (T + 15))
        inf[mask_nuc] = 0.0007 * np.exp(0.28 * (T_c[mask_nuc] + 15.0))
        
        # Saturate INF at 1.0 (100%)
        inf = np.minimum(inf, 1.0)

        # Assuming optimal size range (50-200nm) -> f_size = 1.0 for simplified bulk AgI
        f_size = 1.0
        
        N_inp_agi = N_agi * inf * f_size
        
        # Competition Effect (Eq 3.6)
        # beta = N_nat / (N_nat + 10^4)
        beta = N_inp_nat / (N_inp_nat + 1.0e4)
        
        # Total INP available to form ice
        # N_total = N_nat + N_agi * (1 - beta)
        # But this function returns the Target N_i from AgI specifically? 
        # Or total target?
        # The equation in guide says N_i,total.
        # Let's return the Total Target Ni.
        
        N_target = N_inp_nat + N_inp_agi * (1.0 - beta)
        
        return N_target

    def compute_secondary_production(self, 
                                   T: np.ndarray, 
                                   qc: np.ndarray, qs: np.ndarray, qi: np.ndarray,
                                   ni: np.ndarray, ns: np.ndarray) -> np.ndarray:
        """
        Secondary Ice Production (SIP).
        Ref: Eq 3.5.3
        """
        T_c = T - 273.15
        dni_dt = np.zeros_like(T)
        
        # 1. Hallett-Mossop (-3 to -8 C)
        # dNi/dt_HM = C_HM * R_rime * f(T)
        # Simplified R_rime ~ qc * qs assumption for prototype
        c_hm = 3.5e8
        mask_hm = (T_c > -8.0) & (T_c < -3.0)
        
        if np.any(mask_hm):
            f_t = (T_c[mask_hm] + 8.0) / 5.0
            # Rough approx of riming rate for SIP
            r_rime = 1.0 * qc[mask_hm] * qs[mask_hm] 
            dni_dt[mask_hm] += c_hm * r_rime * f_t

        # 2. Collisional Breakup (Ice-Ice)
        # dNi/dt_coll = C_coll * Ni * Ns * |dv| * phi(T)
        c_coll = 9e-5
        mask_coll = (T_c > -27.0) & (T_c < -3.0)
        
        if np.any(mask_coll):
            phi = 0.5
            dv = 1.0 # Placeholder for relative velocity diff
            dni_dt[mask_coll] += c_coll * ni[mask_coll] * ns[mask_coll] * dv * phi
            
        return dni_dt

    def compute_rates(self,
                      T: np.ndarray, rho: np.ndarray,
                      qc: np.ndarray, qr: np.ndarray, 
                      qi: np.ndarray, qs: np.ndarray, qg: np.ndarray,
                      ni: np.ndarray, ns: np.ndarray, ng: np.ndarray,
                      N_agi: np.ndarray, N_inp_nat: np.ndarray) -> Tuple[np.ndarray, ...]:
        """
        Main interface for ice physics rates.
        Returns tendencies for all hydrometeors.
        """
        dqc = np.zeros_like(qc)
        dqr = np.zeros_like(qr)
        dqi = np.zeros_like(qi)
        dqs = np.zeros_like(qs)
        dqg = np.zeros_like(qg)
        
        dni = np.zeros_like(ni)
        dns = np.zeros_like(ns)
        dng = np.zeros_like(ng)
        
        # 1. Nucleation (AgI + Natural) -> Ice Crystals
        # This is a source term for Ni and qi (via deposition)
        n_inp_target = self.compute_nucleation_agi(T, N_agi, N_inp_nat)
        # Rate to reach target (simple relaxation)
        # In a real model, we track activated nuclei. Here, we nucleate if Ni < Ninp
        nucleation_rate = np.maximum(n_inp_target - ni, 0) / 10.0 # 10s relaxation
        
        dni += nucleation_rate
        # Initial mass of nucleated crystals (e.g. 1e-12 kg)
        dqi += nucleation_rate * 1e-12
        
        # 2. Secondary Production
        sip_rate = self.compute_secondary_production(T, qc, qs, qi, ni, ns)
        dni += sip_rate
        dqi += sip_rate * 1e-12

        # 3. Melting (T > 0C)
        # Eq 3.5.5
        mask_melt = (T > 273.15)
        # Simple melting: Ice/Snow/Graupel -> Rain
        # Rate proportional to superheating
        melt_rate_i = np.zeros_like(T)
        melt_rate_i[mask_melt] = qi[mask_melt] * (T[mask_melt] - 273.15) * 0.1
        
        dqi -= melt_rate_i
        dqr += melt_rate_i # Conservation
        
        # 4. Homogeneous Freezing (T < -40C)
        mask_fz = (T < 233.15)
        fz_rate = np.zeros_like(T)
        fz_rate[mask_fz] = qc[mask_fz] / 1.0 # Instant conversion
        
        dqc -= fz_rate
        dqi += fz_rate
        
        return dqc, dqr, dqi, dqs, dqg, dni, dns, dng
