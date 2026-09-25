"""Guard scientifically consequential endpoint and treatment-selection errors."""
import copy
import importlib.util
import json
from pathlib import Path
import unittest

try:
    import numpy as np
    import pandas
    import matplotlib
except ImportError as exc:
    raise unittest.SkipTest('Scientific-runtime integration checks; run locally with the analysis environment') from exc

BASE = Path(__file__).resolve().parents[2] / 'RQ_Specified/A1_transitional_epithelial_state_distinction'


def module(name, filename):
    spec = importlib.util.spec_from_file_location(name, BASE / 'scripts' / filename)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


lineage = module('lineage', '06_lineage_source_reanalysis.py')
ire1 = module('ire1', '07_ire1_epithelial_analysis.py')


class LineageExecutionTests(unittest.TestCase):
    def test_zero_labelled_cells_is_not_zero_percent(self):
        self.assertTrue(np.isnan(lineage.fraction(0, 0)))
        self.assertEqual(lineage.fraction(0, 12), 0)
        self.assertEqual(lineage.fraction(3, 12), 25)
        with self.assertRaises(ValueError):
            lineage.fraction(3, 2)

    def inputs(self):
        c = json.loads((BASE / 'config/ire1_kira8.json').read_text())
        m = json.loads((BASE / 'metadata/GSE190821.json').read_text(encoding='utf-8'))
        cols = ['Bleo_017_E', 'Bleo_031_E', 'Bleo_045_E', 'Bleo_KIRA8_020_E',
                'Bleo_KIRA8_032_E', 'Bleo_KIRA8_041_E', 'Bleo_E_148', 'Bleo_E_149',
                'Bleo_KIRA8_E_154', 'Bleo_KIRA8_E_155', 'Bleo_E_208R', 'Bleo_017_I']
        return m, cols, c

    def test_antibody_control_and_whole_lung_are_excluded(self):
        m, cols, c = self.inputs()
        result = ire1.select_samples(m, list(reversed(cols)), c)
        self.assertEqual(len(result), 10)
        self.assertNotIn('208R', result.mouse.tolist())
        self.assertNotIn('Bleo_017_I', result.sample_id.tolist())
        self.assertEqual(result.groupby('group').size().to_dict(), {'KIRA8': 5, 'Vehicle': 5})

    def test_missing_column_or_wrong_treatment_fails(self):
        m, cols, c = self.inputs()
        with self.assertRaises(ValueError):
            ire1.select_samples(m, cols[1:], c)
        cols = [x.replace('Bleo_KIRA8_020_E', 'Bleo_020_E') for x in cols]
        with self.assertRaises(ValueError):
            ire1.select_samples(m, cols, c)

    def test_duplicate_library_cannot_inflate_animals(self):
        m, cols, c = self.inputs()
        bad = copy.deepcopy(m)
        bad['samples'].append(copy.deepcopy(bad['samples'][0]))
        with self.assertRaises(ValueError):
            ire1.select_samples(bad, cols, c)
