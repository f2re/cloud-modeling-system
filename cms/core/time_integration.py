from typing import Callable, List
import numpy as np

class RK5Integrator:
    """
    5-stage Strong Stability Preserving Runge-Kutta (SSP-RK5) integrator.
    Ref: IMPLEMENTATION_GUIDE.md Section 6.2
    """
    def __init__(self):
        # Coefficients from IMPLEMENTATION_GUIDE.md Section 6.2
        self.alpha = [0.37, 0.38, 0.18, 0.62, 1.0]

    def step(self, 
             q, 
             dt: float, 
             rhs_func: Callable) -> List[np.ndarray]:
        """
        Advances the state q by time step dt.
        q can be a single np.ndarray or a list of np.ndarrays.
        """
        if isinstance(q, list):
            q_n = [qi.copy() for qi in q]
            q_stage = [qi.copy() for qi in q]

            for a in self.alpha:
                tendencies = rhs_func(q_stage)
                for i in range(len(q_stage)):
                    q_stage[i] = q_n[i] + a * dt * tendencies[i]
            return q_stage
        else:
            q_n = q.copy()
            q_stage = q_n.copy()

            for a in self.alpha:
                tendency = rhs_func(q_stage)
                q_stage = q_n + a * dt * tendency
            return q_stage
