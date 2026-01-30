import unittest
import numpy as np
from cms.utils.validation import compute_statistics, compute_lwp, compute_iwp

class TestValidation(unittest.TestCase):
    def test_statistics(self):
        obs = np.array([1.0, 2.0, 3.0])
        model = np.array([1.1, 1.9, 3.2])
        
        stats = compute_statistics(obs, model)
        
        self.assertAlmostEqual(stats['bias'], 0.0666666, places=5)
        self.assertGreater(stats['correlation'], 0.9)
        self.assertLess(stats['rmse'], 0.2)

    def test_lwp_calculation(self):
        # Grid: 2x2x5 (z is last dim)
        qc = np.ones((2, 2, 5)) * 0.001
        qr = np.zeros((2, 2, 5))
        rho = np.ones((2, 2, 5)) * 1.0
        dz = 100.0
        
        lwp = compute_lwp(qc, qr, rho, dz)
        
        # Expected: sum(0.001 * 1.0) over 5 levels * 100 = 0.005 * 100 = 0.5 kg/m^2
        self.assertEqual(lwp[0, 0], 0.5)

if __name__ == '__main__':
    unittest.main()
