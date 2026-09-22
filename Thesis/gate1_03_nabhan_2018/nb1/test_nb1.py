"""Independent combinatorial checks for the depth-standardization calculation."""
import math
import unittest
import numpy as np
from run_nb1 import expected_detection


class DetectionTest(unittest.TestCase):
    def test_against_enumerated_combinations(self):
        for n in [4,7,12]:
            for budget in [1,3,4]:
                for k in range(n+1):
                    expected=1-(math.comb(n-k,budget)/math.comb(n,budget) if n-k>=budget else 0)
                    self.assertAlmostEqual(float(expected_detection(n,k,budget)),expected,places=12)

    def test_invalid_library_rejected(self):
        with self.assertRaises(ValueError):
            expected_detection(2,1,3)

    def test_broadcast_and_monotonicity(self):
        p=expected_detection(np.array([1000,2000])[:,None],np.array([0,1,5,20]))
        self.assertEqual(p.shape,(2,4))
        self.assertTrue(np.all(np.diff(p,axis=1)>=0))
        self.assertTrue(np.all(p[0]>=p[1]))


if __name__=="__main__":
    unittest.main()
