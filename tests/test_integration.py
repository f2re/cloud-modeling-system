import unittest
import os
import shutil
import numpy as np
import netCDF4 as nc
from cms.model import CMSModel
from cms.config import GridConfig, PhysicsConfig
from cms.utils.io import IOManager

class TestIntegration(unittest.TestCase):
    def setUp(self):
        # Create a small grid for fast testing
        self.grid_config = GridConfig(nx=20, ny=20, nz=20)
        self.physics_config = PhysicsConfig()
        self.model = CMSModel(self.grid_config, self.physics_config)
        self.io = IOManager(output_dir="test_output")
        
    def tearDown(self):
        # Clean up test output directory
        if os.path.exists("test_output"):
            shutil.rmtree("test_output")

    def test_simulation_loop_and_output(self):
        """
        Runs the model for a few steps and verifies output generation.
        """
        # 1. Initialize a warm bubble
        center_x, center_y, center_z = 10, 10, 5
        grid = self.model.grid
        
        # Add thermal perturbation
        # Indices for inner domain + ghost cells
        ic_x = grid.ng + center_x
        ic_y = grid.ng + center_y
        ic_z = grid.ng + center_z
        
        self.model.theta[ic_x-2:ic_x+3, ic_y-2:ic_y+3, ic_z-2:ic_z+3] += 2.0
        
        # Add some cloud water to test microphysics
        self.model.qc[ic_x, ic_y, ic_z] = 1e-3 # 1 g/kg
        self.model.nc[ic_x, ic_y, ic_z] = 1e8  # 100 cm-3 (Standard value)

        # 2. Run simulation
        dt = 1.0
        steps = 5
        for n in range(steps):
            self.model.step(dt)

        # 3. Save output
        filename = "cms_test_001.nc"
        filepath = self.io.save_state(
            filename, self.model.grid, steps * dt,
            self.model.u, self.model.v, self.model.w,
            self.model.rho, self.model.theta,
            self.model.qc, self.model.qr, self.model.nc,
            self.model.qi, self.model.qs, self.model.qg,
            self.model.ni, self.model.ns, self.model.ng,
            self.model.c_agi, self.model.n_ccn, self.model.n_inp_nat
        )

        # 4. Verify file
        self.assertTrue(os.path.exists(filepath))
        
        with nc.Dataset(filepath, 'r') as ds:
            # Check dimensions
            self.assertEqual(ds.dimensions['x'].size, 20)
            self.assertEqual(ds.dimensions['z'].size, 20)
            
            # Check variables exist
            self.assertIn('w', ds.variables)
            self.assertIn('qr', ds.variables)
            
            # Check physics happened (rain formation from autoconversion)
            qr_max = np.max(ds.variables['qr'][:])
            # We initialized with qc=1g/kg, so some should convert to rain
            self.assertGreater(qr_max, 0.0)
            
            # Check dynamics happened (buoyancy generated w)
            w_max = np.max(ds.variables['w'][:])
            self.assertGreater(w_max, 0.0)

if __name__ == '__main__':
    unittest.main()
