"""Failure-focused regression tests for cohort, parser and denominator choices."""
import gzip
import unittest

import numpy as np
import pandas as pd

import reproduce


class ReproductionTests(unittest.TestCase):
    def test_headerless_first_gene_is_not_lost(self):
        line = "FirstGene\t-\t-\tFirstGene\tFirstGene\tTSS1\tchr1:1-2\t-\t-\t1.0\t0\t2\tOK\n"
        frame = reproduce.read_tracking(gzip.compress(line.encode()))
        self.assertEqual(len(frame), 1)
        self.assertEqual(frame.iloc[0]["gene_short_name"], "FirstGene")
        self.assertEqual(frame.iloc[0]["FPKM"], 1.0)

    def test_bulk_excluded_even_if_strong_expression(self):
        matrix = pd.DataFrame({"Wnt5a": [0.0, 1000.0]}, index=["cell", reproduce.BULK])
        metadata = pd.DataFrame({"accession": ["cell", reproduce.BULK], "is_bulk": [False, True]})
        cells, bulk = reproduce.split_samples(matrix, metadata)
        self.assertEqual(list(cells.index), ["cell"])
        self.assertEqual(cells["Wnt5a"].sum(), 0.0)
        self.assertEqual(list(bulk.index), [reproduce.BULK])

    def test_detection_boundary(self):
        values = np.array([0, 0.0001, 0.1, 0.9999, 1, 5])
        np.testing.assert_equal(reproduce.detection(values, 0), [False, True, True, True, True, True])
        np.testing.assert_equal(reproduce.detection(values, 1), [False, False, False, False, True, True])

    def test_conditional_denominator_is_wnt_positive_not_all_cells(self):
        joint, denominator, fraction = reproduce.conditional_fraction(
            np.array([True, False, True, True]), np.array([True, True, False, False]))
        self.assertEqual((joint, denominator, fraction), (1, 2, 0.5))
        self.assertEqual(reproduce.conditional_fraction(np.array([True]), np.array([False])), (0, 0, None))


if __name__ == "__main__":
    unittest.main()
