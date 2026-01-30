import numpy as np
from typing import Tuple
from cms.config import PhysicsConfig

class WarmMicrophysics:
    """
    Double-moment warm rain microphysics (Kessler/Morrison style).
    Ref: IMPLEMENTATION_GUIDE.md Section 3.5.1
    """
    def __init__(self, config: PhysicsConfig):
        self.config = config

    def compute_rates(self, 
                      qc: np.ndarray, qr: np.ndarray, 
                      nc: np.ndarray, nr: np.ndarray, 
                      rho: np.ndarray) -> Tuple[np.ndarray, ...]:
        """
        Computes transformation rates for warm microphysics.
        Returns (dqc, dqr, dnc, dnr) tendencies.
        """
        # 1. Autoconversion (cloud -> rain)
        # Eq 3.5.1: (dqr/dt)_auto = 1350 * qc^2.47 * nc^-1.79 / rho^1.47
        # Note: nc is number concentration (m^-3)
        auto_rate = np.zeros_like(qc)
        mask = qc > 1e-6
        auto_rate[mask] = 1350.0 * qc[mask]**2.47 * nc[mask]**-1.79 / rho[mask]**1.47
        
        dqr_auto = auto_rate
        dqc_auto = -auto_rate
        
        # 2. Accretion (cloud by rain)
        # Eq 3.5.1: (dqr/dt)_accr = (pi/4) * E_cr * nr * qc * Dr^2 * vr
        # Using a simplified version or full power law
        e_cr = 1.0
        
        # Safety: Clamp inputs to avoid negative roots
        qr_safe = np.maximum(qr, 0.0)
        nr_safe = np.maximum(nr, 1e-12) # Avoid division by zero
        
        # Placeholder for Dr (mean diameter) and vr (terminal velocity)
        # In a real double-moment scheme, these depend on (qr, nr)
        # Ref Section 3.4 for full parameterization
        dr = (6.0 * rho * qr_safe / (np.pi * self.config.rho_w * nr_safe))**(1/3)
        vr = 36.34 * dr**0.5  # Simplified power law
        
        accr_rate = (np.pi / 4.0) * e_cr * nr * qc * dr**2 * vr
        
        dqr_accr = accr_rate
        dqc_accr = -accr_rate

        # 3. Self-collection (rain number reduction)
        # Eq 3.5.1: (dnr/dt)_self = -5.78 * nr^2 * dr^3
        dnr_self = -5.78 * nr**2 * dr**3

        # Combine tendencies
        dqc = dqc_auto + dqc_accr
        dqr = dqr_auto + dqr_accr
        dnc = np.zeros_like(nc) # Simplified: nc constant for now
        dnr = (auto_rate / (rho * 1e-10)) + dnr_self # Simplified conversion to number

        return dqc, dqr, dnc, dnr
