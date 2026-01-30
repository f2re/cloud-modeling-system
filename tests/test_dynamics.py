import unittest
import numpy as np
from cms.core.grid import Grid
from cms.config import GridConfig, PhysicsConfig
from cms.dynamics.navier_stokes import NavierStokesSolver
from cms.dynamics.boundary import BoundaryConditions

class TestDynamics(unittest.TestCase):
    def setUp(self):
        self.grid_config = GridConfig(nx=20, ny=20, nz=20)
        self.physics_config = PhysicsConfig()
        self.grid = Grid(self.grid_config)
        self.solver = NavierStokesSolver(self.grid, self.physics_config)
        self.bc = BoundaryConditions(self.grid)

    def test_tendencies_initialization(self):
        """Tests that tendencies can be computed without crashing."""
        shape = self.grid.shape
        u = np.zeros(shape)
        v = np.zeros(shape)
        w = np.zeros(shape)
        rho = np.ones(shape) * self.physics_config.rho_a_stp
        theta = np.ones(shape) * 300.0
        p = np.ones(shape) * 101325.0

        # Add a small perturbation
        theta[10, 10, 5] += 1.0

        du, dv, dw, drho, dtheta = self.solver.compute_tendencies(u, v, w, rho, theta, p)
        
        # Check that buoyancy generated vertical momentum tendency
        self.assertTrue(np.any(dw != 0))
        # Check that advection didn't create NANs
        self.assertFalse(np.any(np.isnan(du)))

if __name__ == '__main__':
    unittest.main()
