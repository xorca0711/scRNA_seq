"""Failure cases that guard biological units and signed evidence."""
import importlib.util
from pathlib import Path
import unittest

PATH = Path(__file__).resolve().parents[2] / 'RQ_Specified/A1_transitional_epithelial_state_distinction/scripts/a1_regulatory_fate.py'
SPEC = importlib.util.spec_from_file_location('a1_regulatory_fate', PATH)
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)


class RegulatoryFateTests(unittest.TestCase):
    def test_unequal_field_denominators_do_not_become_equal_cell_weights(self):
        self.assertAlmostEqual(M.pooled_fraction([1, 90], [10, 100]), 91/110)
        self.assertNotEqual(M.pooled_fraction([1, 90], [10, 100]), .5)

    def test_zero_denominator_is_not_zero_outcome(self):
        with self.assertRaises(AssertionError):
            M.pooled_fraction([0], [0])
        with self.assertRaises(AssertionError):
            M.pooled_fraction([11], [10])
        self.assertEqual(M.pooled_fraction([0], [10]), 0)

    def test_two_sided_three_mouse_resolution(self):
        result = M.exact_permutation([-3, -2, -1], [1, 2, 3])
        self.assertEqual(result['permutations'], 20)
        self.assertEqual(result['extreme'], 2)
        self.assertEqual(result['p_two_sided'], .1)
        self.assertEqual(M.exact_permutation([0, 0, 0], [0, 0, 0])['p_two_sided'], 1)

    def test_regional_difference_is_not_an_overall_marker_change(self):
        # WT de-novo minus in-situ = -10; mutant = +10.
        self.assertEqual(M.region_interaction([-10]*3, [10]*3), 20)

    def test_shared_significance_does_not_imply_shared_direction(self):
        self.assertEqual(M.sign_class(2, -3), 'opposite')
        self.assertEqual(M.sign_class(-2, -3), 'both_down')
        self.assertEqual(M.sign_class(2, 3), 'both_up')


if __name__ == '__main__':
    unittest.main()
