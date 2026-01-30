import unittest
import numpy as np
from cms.utils.diagnostics import compute_radar_reflectivity
from cms.config import PhysicsConfig

class TestDiagnostics(unittest.TestCase):
    def setUp(self):
        self.config = PhysicsConfig()
        self.shape = (10, 10, 10)

    def test_reflectivity_rain_only(self):
        """Tests radar reflectivity for rain only."""
        rho = np.ones(self.shape) * 1.0
        qr = np.ones(self.shape) * 1e-3  # 1 g/kg
        nr = np.ones(self.shape) * 1e3   # 1 per Liter (1000 m^-3)
        
        dbz = compute_radar_reflectivity(rho, qr, nr, config=self.config)
        
        # Expected value calculation (approximate):
        # Z = (720 * (1*1e-3)^2) / (pi^2 * 1000^2 * 1000) * 1e18
        # Z = (720 * 1e-6) / (9.8696 * 1e6 * 1000) * 1e18
        # Z = (720 / 9.8696) * 1e-15 * 1e18 = 72.95 * 10^3 = 72950 mm^6/m^3
        # dBZ = 10 * log10(72950) = 48.63
        
        expected_dbz = 48.63
        self.assertAlmostEqual(dbz[0, 0, 0], expected_dbz, places=1)

    def test_reflectivity_empty(self):
        """Tests that empty fields return floor value."""
        rho = np.ones(self.shape) * 1.0
        qr = np.zeros(self.shape)
        nr = np.zeros(self.shape)
        
        dbz = compute_radar_reflectivity(rho, qr, nr, config=self.config)
        
        # Floor value is 10 * log10(1e-3) = -30 dBZ
        self.assertEqual(dbz[0, 0, 0], -30.0)

    def test_reflectivity_mixed(self):
        """Tests that mixed species increase reflectivity."""
        rho = np.ones(self.shape) * 1.0
        qr = np.ones(self.shape) * 1e-3
        nr = np.ones(self.shape) * 1e3
        
        dbz_rain = compute_radar_reflectivity(rho, qr, nr, config=self.config)
        
        # Add ice
        qi = np.ones(self.shape) * 1e-4
        ni = np.ones(self.shape) * 1e4
        
        dbz_mixed = compute_radar_reflectivity(rho, qr, nr, qi, ni, config=self.config)
        
        # Reflectivity should be higher
        self.assertTrue(np.all(dbz_mixed > dbz_rain))

if __name__ == '__main__':
    unittest.main()
