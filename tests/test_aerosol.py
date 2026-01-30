import unittest
import numpy as np
from cms.model import CMSModel
from cms.config import GridConfig, PhysicsConfig

class TestAerosol(unittest.TestCase):
    def setUp(self):
        self.grid_config = GridConfig(nx=20, ny=20, nz=20)
        self.physics_config = PhysicsConfig()
        self.model = CMSModel(self.grid_config, self.physics_config)
        
    def test_ccn_activation(self):
        """
        Tests that CCN activates into droplets in an updraft.
        """
        center = self.model.grid.ng + 10
        
        # 1. Setup Environment: Updraft
        # Use uniform updraft to avoid strong divergence/advection effects
        self.model.w[:] = 5.0 
        
        # Initial state
        center = self.model.grid.ng + 10
        initial_ccn = self.model.n_ccn[center, center, center]
        initial_nc = self.model.nc[center, center, center] # Likely 0 or small init
        
        # 2. Step
        dt = 1.0
        self.model.step(dt)
        
        # 3. Check Activation
        final_ccn = self.model.n_ccn[center, center, center]
        final_nc = self.model.nc[center, center, center]
        
        # CCN should decrease
        self.assertLess(final_ccn, initial_ccn)
        
        # Nc should increase
        self.assertGreater(final_nc, initial_nc)
        
        # Conservation (roughly, ignoring advection for 1 step at center)
        # Delta_CCN + Delta_Nc ~ 0
        delta_ccn = initial_ccn - final_ccn
        delta_nc = final_nc - initial_nc
        self.assertAlmostEqual(delta_ccn, delta_nc, delta=delta_ccn*0.1)

if __name__ == '__main__':
    unittest.main()
