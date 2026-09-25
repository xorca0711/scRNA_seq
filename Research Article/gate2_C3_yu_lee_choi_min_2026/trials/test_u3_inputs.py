"""Count-integrity tests for the new acquisition/QC implementation."""
import gzip
from pathlib import Path
import tempfile
import unittest
from u3_acquire_qc import read_counts


class CountIntegrity(unittest.TestCase):
    def fixture(self, tmp, matrix, features='ENSM1\tGene1\tGene Expression\nENSM2\tGene2\tGene Expression\n', barcodes='AA\nBB\n'):
        result = {}
        for kind, content in [('matrix', matrix), ('features', features), ('barcodes', barcodes)]:
            path = Path(tmp)/(kind+'.gz')
            with gzip.open(path, 'wt', encoding='utf-8') as handle: handle.write(content)
            result[kind] = path
        return result

    def test_duplicate_coordinates_are_collapsed_without_losing_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            files = self.fixture(tmp, '%%MatrixMarket matrix coordinate integer general\n2 2 3\n1 1 2\n1 1 3\n2 2 7\n')
            x, features, barcodes, duplicates = read_counts(files)
            self.assertEqual(duplicates, 1)
            self.assertEqual(x.toarray().tolist(), [[5, 0], [0, 7]])

    def test_fractional_or_negative_counts_rejected(self):
        for value in ['0.5', '-1']:
            with self.subTest(value=value), tempfile.TemporaryDirectory() as tmp:
                files = self.fixture(tmp, '%%MatrixMarket matrix coordinate real general\n2 2 1\n1 1 '+value+'\n')
                with self.assertRaisesRegex(ValueError, 'nonnegative integer'): read_counts(files)

    def test_duplicate_barcodes_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            files = self.fixture(tmp, '%%MatrixMarket matrix coordinate integer general\n2 2 1\n1 1 1\n', barcodes='AA\nAA\n')
            with self.assertRaisesRegex(ValueError, 'Nonunique'): read_counts(files)


if __name__ == '__main__':
    unittest.main()
