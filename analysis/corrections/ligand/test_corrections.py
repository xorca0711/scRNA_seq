"""Targeted checks for the scientifically consequential corrected definitions."""
import unittest
import gzip
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

from common import allowed_pairs, eligible_donors, expected_detection_at_budget
from c37_annotation_depth import extract_depth


class CorrectionTests(unittest.TestCase):
    def test_mural_cells_neither_supply_target_floor_nor_become_senders(self):
        obs = pd.DataFrame({"donor": ["d"] * 150,
            "category": ["Epithelial"] * 50 + ["Stromal"] * 100,
            "celltype": ["ATII"] * 50 + ["Fibroblast"] * 20 + ["Pericyte"] * 80})
        self.assertFalse(eligible_donors(obs).loc["d", "eligible"])
        self.assertEqual(allowed_pairs(obs).to_dict("records"), [{"source": "ATII", "target": "Fibroblast"}])

    def test_actual_combined_target_floor(self):
        obs = pd.DataFrame({"donor": ["d"] * 100,
            "category": ["Epithelial"] * 50 + ["Stromal"] * 50,
            "celltype": ["ATII"] * 50 + ["Fibroblast"] * 30 + ["Myofibroblast"] * 20})
        self.assertTrue(eligible_donors(obs).loc["d", "eligible"])

    def test_exact_molecule_detection_against_scipy_distribution(self):
        n = np.array([100, 100, 100, 1000, 5000, 8])
        k = np.array([0, 1, 99, 4, 25, 1])
        actual = expected_detection_at_budget(n, k, 50)
        np.testing.assert_allclose(actual[:-1], hypergeom.sf(0, n[:-1], k[:-1], 50), rtol=1e-8, atol=1e-9)
        self.assertTrue(np.isnan(actual[-1]))
        np.testing.assert_allclose(expected_detection_at_budget(n, k, 1), k / n)

    def test_cannot_use_more_target_molecules_than_total(self):
        with self.assertRaises(ValueError):
            expected_detection_at_budget([10], [11], 5)

    def test_raw_extraction_preserves_annotation_order_and_all_gene_totals(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            source = root / "raw.txt.gz"
            with gzip.open(source, "wt") as handle:
                handle.write("Index\tc1\tc2\tc3\nAREG\t1\t0\t3\nGENE2\t4\t9\t0\nGENE3\t0\t2\t5\n")
            meta = pd.DataFrame(index=pd.Index(["c3", "c1"], name="Index"))
            result = extract_depth(meta, source, root / "selected.csv.gz", "fixture")
            self.assertEqual(result.index.tolist(), ["c3", "c1"])
            self.assertEqual(result.total_umi.tolist(), [8, 5])
            self.assertEqual(result.n_genes.tolist(), [2, 2])
            self.assertEqual(result.areg_umi.tolist(), [3, 1])


if __name__ == "__main__":
    unittest.main()
