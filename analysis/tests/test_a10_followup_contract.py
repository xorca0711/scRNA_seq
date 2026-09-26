"""Portable checks of the delivered A10 evidence, without raw caches or NumPy."""
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
HERE = ROOT / 'RQ_Specified/A10_organoid_growth_outcome'
OUT = HERE / 'tables/followup_v1'


def rows(name):
    with (OUT / name).open(newline='', encoding='utf-8') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


class A10FollowupContract(unittest.TestCase):
    def test_frozen_sources_and_output_integrity(self):
        for name, script in [('diagnostic_run.json','07_audit_followup.py'),('model_run.json','08_fit_followup.py')]:
            record = json.loads((OUT / name).read_text())
            self.assertEqual(hashlib.sha256((HERE/'scripts'/script).read_bytes()).hexdigest(),record['script_sha256'])
            for filename, digest in record['outputs'].items():
                with self.subTest(filename=filename):
                    self.assertEqual(hashlib.sha256((OUT/filename).read_bytes()).hexdigest(),digest)
        record = json.loads((OUT/'model_run.json').read_text())
        self.assertEqual(hashlib.sha256((HERE/'config/a10_followup_models.json').read_bytes()).hexdigest(),record['spec_sha256'])

    def test_all_sample_ids_and_target_confounding_are_visible(self):
        samples = rows('diagnostic_geo_sample_fields.tsv')
        self.assertEqual(len(samples),886)
        self.assertEqual(len({r['gsm'] for r in samples}),886)
        self.assertEqual(len({r['library_name'] for r in samples}),886)
        record = json.loads((OUT/'diagnostic_run.json').read_text())
        self.assertEqual(record['structured_preparation_fields'],[])
        self.assertEqual(set(record['shared_target_names']),{'MECOM','RNF43','TDTOMATO','TIGIT'})
        plate3 = [r for r in rows('diagnostic_target_overlap.tsv') if r['scope']=='drop_top_and_TDTOMATO'
                  and 'plate3' in (r['plate_a'],r['plate_b'])]
        self.assertEqual(len(plate3),3)
        self.assertTrue(all(int(r['common_targets'])==0 for r in plate3))

    def test_saved_predictions_support_the_error_and_decision_records(self):
        spec = json.loads((HERE/'config/a10_followup_models.json').read_text())
        error = defaultdict(float)
        seen = set()
        for row in rows('model_predictions.tsv'):
            key = (row['setting'],row['scale'],row['evaluation'])
            identity = key + (row['library name'],)
            self.assertNotIn(identity,seen)
            seen.add(identity)
            for model in spec['models']:
                error[key+(model,)] += (float(row['outcome'])-float(row[model]))**2
        self.assertEqual(len(seen),13232)
        comparison = {r['name']:r for r in spec['comparisons']}
        summary = rows('model_comparisons.tsv')
        for row in summary:
            key = (row['setting'],row['scale'],row['evaluation'])
            item = comparison[row['comparison']]
            expected = 1-error[key+(item['candidate'],)]/error[key+(item['reference'],)]
            self.assertTrue(math.isclose(expected,float(row['relative_error_reduction']),rel_tol=1e-7,abs_tol=1e-8))
        run = json.loads((OUT/'model_run.json').read_text())
        margin = spec['decisions']['descriptive_margin_relative_error_reduction']
        for name, verdict in run['verdicts'].items():
            chosen = [r for r in summary if r['comparison']==name and r['setting'] in ('primary','drop_both')]
            within = [r for r in chosen if r['evaluation']=='within_group']
            plate = [r for r in chosen if r['evaluation']=='plate_shift']
            self.assertEqual((len(within),len(plate)),(4,4))
            passed = all(float(r['relative_error_reduction'])>=margin for r in within)
            consistent = passed and all(float(r['relative_error_reduction'])>=margin and float(r['minimum_fold_reduction'])>0 for r in plate)
            self.assertEqual(verdict,{'within_group_margin_met':passed,'consistent_plate_shift_gain':consistent})

    def test_relative_improvement_does_not_hide_absolute_failure(self):
        record = json.loads((OUT/'verification.json').read_text())
        self.assertEqual(hashlib.sha256((OUT/'model_absolute_performance.tsv').read_bytes()).hexdigest(),record['absolute_diagnostic_sha256'])
        self.assertEqual(record['independent_qr_fits'],32)
        self.assertLess(record['independent_qr_max_prediction_error'],1e-8)
        full = [r for r in rows('model_absolute_performance.tsv') if r['setting']=='primary' and r['scale']=='inherited_log2'
                and r['model']=='baseline+proliferation+remaining_growth']
        self.assertEqual(len(full),4)
        self.assertEqual(sum(float(r['r2_against_heldout_mean'])<0 for r in full),3)


if __name__ == '__main__':
    unittest.main()
