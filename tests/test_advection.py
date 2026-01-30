import unittest
import numpy as np
from cms.core.grid import Grid
from cms.core.time_integration import RK5Integrator
from cms.config import GridConfig

class TestAdvection(unittest.TestCase):
    def setUp(self):
        self.config = GridConfig(nx=50, ny=10, nz=10, dx=1.0, dy=1.0, dz=1.0)
        self.grid = Grid(self.config)
        self.integrator = RK5Integrator()

    def test_linear_advection_1d(self):
        """
        Tests linear advection of a Gaussian pulse in the X direction.
        """
        u_velocity = 0.1  # m/s
        dt = 0.5         # s
        n_steps = 20
        
        # Initialize field with a Gaussian pulse
        field = self.grid.create_field()
        x_coords = np.arange(self.grid.gnx) * self.grid.dx
        center = self.grid.ng + 10
        field[:, :, :] = np.exp(-((x_coords[:, None, None] - center)**2) / 10.0)

        initial_pulse_sum = np.sum(field[self.grid.inner])

        def rhs_advection(q):
            # Simple 1st order upwind advection: -u * dq/dx
            dq_dx = np.zeros_like(q)
            # dq/dx approx (q[i] - q[i-1]) / dx
            dq_dx[1:, :, :] = (q[1:, :, :] - q[:-1, :, :]) / self.grid.dx
            return -u_velocity * dq_dx

        # Integrate in time
        curr_field = field.copy()
        for _ in range(n_steps):
            curr_field = self.integrator.step(curr_field, dt, rhs_advection)

        # Check that the pulse moved to the right
        final_pulse_sum = np.sum(curr_field[self.grid.inner])
        
        # Basic conservation check (allowing for some numerical diffusion/boundary effects)
        self.assertAlmostEqual(initial_pulse_sum, final_pulse_sum, delta=initial_pulse_sum * 0.1)
        
        # Check displacement: pulse should move approx u * dt * n_steps = 0.1 * 0.5 * 20 = 1.0 meter
        initial_peak_idx = np.argmax(field[:, 5, 5])
        final_peak_idx = np.argmax(curr_field[:, 5, 5])
        
        self.assertGreater(final_peak_idx, initial_peak_idx)

if __name__ == '__main__':
    unittest.main()
