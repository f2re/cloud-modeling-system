import numpy as np
from typing import Tuple
from cms.config import PhysicsConfig

class Activation:
    """
    Aerosol-Cloud Interactions: CCN Activation.
    Ref: IMPLEMENTATION_GUIDE.md Section 3.6
    """
    def __init__(self, config: PhysicsConfig):
        self.config = config
        # Constants for continental aerosol (Ref: Guide Section 3.6)
        self.k1 = 0.7
        self.k2 = 1.5
        # C1, C2 need to be tuned or standard values. 
        # Using typical values for Twomey-style power law derived Smax
        # Smax ~ w^3/4 * Nccn^-1/2 is a simplified relationship.
        # We'll use a standard parameterization wrapper here.
        
    def compute_activation(self, 
                           w: np.ndarray, 
                           n_ccn: np.ndarray, 
                           nc: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        Computes droplet activation rate.
        Returns (dnc_dt, dn_ccn_dt).
        """
        # 1. Compute Supersaturation Max (S_max)
        # S_max = C1 * w^0.75 * N_ccn^-0.5
        # We only activate in updrafts (w > 0)
        w_up = np.maximum(w, 0.0)
        
        # Avoid division by zero in N_ccn
        n_ccn_safe = np.maximum(n_ccn, 1e6) # Min 1/cm^3
        
        # Coefficients (Placeholder values consistent with Guide's proportionalities)
        c1 = 100.0 # Scaling factor for S_max equation (tuned for realistic S%)
        
        s_max = c1 * (w_up**0.75) * (n_ccn_safe**-0.5)
        
        # 2. Compute Activated Number (N_act)
        # N_act = C2 * N_ccn^k1 * S_max^k2
        # But commonly, N_act = C * S^k (Twomey). 
        # The Guide specifies: N_c_act = C2 * N_ccn^k1 * S_max^k2
        
        c2 = 1.0 # Scaling factor
        n_act = c2 * (n_ccn_safe**self.k1) * (s_max**self.k2)
        
        # Ensure N_act doesn't exceed available CCN
        n_act = np.minimum(n_act, n_ccn_safe)
        
        # 3. Compute Rate
        # Only activate if predicted N_act > current Nc
        # Rate = (N_act - Nc) / tau
        tau = 1.0 # Timescale (seconds)
        
        activation_tendency = np.maximum(n_act - nc, 0.0) / tau
        
        dnc_dt = activation_tendency
        dn_ccn_dt = -activation_tendency # Depletion of aerosol
        
        return dnc_dt, dn_ccn_dt
