import unittest
import numpy as np
from cms.model import CMSModel
from cms.config import GridConfig, PhysicsConfig

class TestIceSeeding(unittest.TestCase):
    def setUp(self):
        self.grid_config = GridConfig(nx=20, ny=20, nz=20)
        self.physics_config = PhysicsConfig()
        self.model = CMSModel(self.grid_config, self.physics_config)
        
    def test_agi_nucleation(self):
        """
        Tests that AgI creates ice in supercooled conditions.
        """
        center = self.model.grid.ng + 10
        
        # 1. Setup Environment: Supercooled (-10 C)
        self.model.theta[:] = 263.15 # -10 C
        
        # 2. Add AgI (Seeding)
        self.model.c_agi[center, center, center] = 1e-9 # 1 microgram/kg
        
        # 3. Step
        dt = 1.0
        # Run for 10 steps to allow nucleation relaxation
        for _ in range(10):
            self.model.step(dt)
            
        # 4. Check for Ice Crystals (Ni)
        ni_max = np.max(self.model.ni)
        qi_max = np.max(self.model.qi)
        
        self.assertGreater(ni_max, 0.0, "AgI failed to nucleate ice crystals (Ni)")
        self.assertGreater(qi_max, 0.0, "AgI failed to generate ice mass (qi)")
        
    def test_sip_hallett_mossop(self):
        """
        Tests Secondary Ice Production (HM process).
        """
        # 1. Setup Environment: -5 C (Optimal HM zone)
        self.model.theta[:] = 268.15 
        center = self.model.grid.ng + 10
        
        # 2. Add Ingredients: Cloud Water + Snow (Riming substrate)
        self.model.qc[:] = 1e-3
        self.model.qs[:] = 1e-4
        
        # 3. Step
        dt = 1.0
        self.model.step(dt)
        
        # 4. Check for new Ice Crystals (Ni)
        # Note: HM produces *Splinters* which we count as Ni
        ni_max = np.max(self.model.ni)
        self.assertGreater(ni_max, 0.0, "Hallett-Mossop failed to produce splinters")

if __name__ == '__main__':
    unittest.main()
