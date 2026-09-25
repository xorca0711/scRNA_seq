"""Hand-calculated coordinate cases; no assay files or scientific packages needed."""
import gzip
import importlib.util
from pathlib import Path
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / '15_quantify_direct_marks.py'
spec = importlib.util.spec_from_file_location('direct_marks', SCRIPT)
direct = importlib.util.module_from_spec(spec)
spec.loader.exec_module(direct)


class IntervalTests(unittest.TestCase):
    def test_half_open_clipping_gaps_and_zero(self):
        self.assertEqual(direct.interval_signal([(0, 4, 2), (6, 12, 3)], 2, 10), (2, .75))
        self.assertEqual(direct.interval_signal([(0, 2, 99), (10, 12, 99)], 2, 10), (0, 0))
        self.assertEqual(direct.interval_signal([(2, 10, 0)], 2, 10), (0, 1))
        self.assertEqual(direct.interval_signal([], 2, 10), (0, 0))

    def test_reject_invalid_records(self):
        for records in [[(0, 6, 1), (5, 10, 2)], [(6, 8, 1), (2, 4, 2)],
                        [(2, 4, -1)], [(2, 4, float('nan'))]]:
            with self.assertRaises(ValueError):
                direct.interval_signal(records, 0, 10)
        with self.assertRaises(ValueError):
            direct.interval_signal([], 2, 2)

    def test_union_nested_adjacent_and_clipped(self):
        self.assertEqual(direct.union_length([(0, 4), (2, 3), (4, 6), (8, 12)], 2, 10), 6)
        self.assertEqual(direct.union_length([(0, 2), (10, 12)], 2, 10), 0)

    def test_strand_exact_symbols_and_ambiguous_loci(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / 'genes.gff.gz'
            with gzip.open(p, 'wt') as h:
                for chrom, left, right, strand, symbol in [
                    ('chr1', 11, 20, '+', 'PLUS'), ('chr2', 31, 40, '-', 'MINUS'),
                    ('chr1', 11, 20, '+', 'DUP'), ('chr2', 31, 40, '+', 'DUP'),
                    ('chr1', 11, 20, '+', 'WRONG_suffix')]:
                    h.write(f'{chrom}\tref\tgene\t{left}\t{right}\t.\t{strand}\t.\tgene_name={symbol}\n')
            rows = {r['gene']: r for r in direct.annotations(p, 'T2T-CHM13v2.0',
                {'test': ['PLUS', 'MINUS', 'DUP', 'WRONG']})}
            self.assertEqual((rows['PLUS']['gene_start'], rows['PLUS']['tss']), (10, 10))
            self.assertEqual((rows['MINUS']['gene_end'], rows['MINUS']['tss']), (40, 39))
            self.assertEqual(rows['DUP']['status'], 'held_missing_or_multiple_loci')
            self.assertEqual(rows['WRONG']['status'], 'held_missing_or_multiple_loci')
            with gzip.open(p, 'wt') as h:
                for left, right in [(10, 30), (20, 40)]:
                    h.write(f'0\ttx\tchr1\t-\t{left}\t{right}\t0\t0\t0\t0\t0\t0\tMINUS\n')
            row = direct.annotations(p, 'hg19', {'test': ['MINUS']})[0]
            self.assertEqual((row['gene_start'], row['gene_end'], row['tss']), (10, 40, 39))


if __name__ == '__main__':
    unittest.main()
