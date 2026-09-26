"""Portable A0 evidence checks: independent units, frozen rules and saved decisions."""
import csv
import gzip
import hashlib
import json
import math
import statistics
import unittest
import zipfile
from pathlib import Path
from collections import defaultdict

BASE = Path(__file__).resolve().parents[2] / 'RQ_Specified/A0_conserved_epithelial_transition_program'
OUT = BASE / 'tables/pilot_v1'


def read_json(path):
    return json.loads(path.read_text(encoding='utf-8'))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class A0PilotContract(unittest.TestCase):
    def test_frozen_rules_and_executed_code_are_bound(self):
        prep = read_json(OUT / 'preparation.json')
        run = read_json(OUT / 'discovery_run.json')
        module = read_json(OUT / 'frozen_programme.json')
        cfg_hash = digest(BASE / 'config/pilot_v1.json')
        self.assertEqual({prep['config_sha256'], run['config_sha256'], module['config_sha256']}, {cfg_hash})
        self.assertEqual(prep['script_sha256'], digest(BASE / 'scripts/11_prepare_pilot.py'))
        self.assertEqual(run['script_sha256'], digest(BASE / 'scripts/12_discover_pilot.py'))
        self.assertEqual(run['preparation_sha256'], digest(OUT / 'preparation.json'))
        self.assertEqual(module['ortholog_universe_sha256'], digest(OUT / 'ortholog_universe.tsv'))
        for filename, sha in run['outputs'].items():
            self.assertEqual(digest(OUT / filename), sha, filename)
        self.assertFalse(prep['programme_effects_computed'])
        self.assertFalse(run['V1_expression_effects_read'])
        verification = read_json(OUT / 'verification.json')
        self.assertEqual(verification['status'], 'PASS')
        self.assertEqual(verification['verifier_sha256'], digest(BASE / 'scripts/14_verify_pilot.py'))
        self.assertGreater(verification['raw_count_comparisons'], 1000)
        self.assertGreater(verification['source_raw_count_comparisons'], 100000)
        self.assertLess(verification['maximum_raw_error'], 1e-9)

    def test_eligibility_counts_biological_units_and_both_endpoints(self):
        cfg = read_json(BASE / 'config/pilot_v1.json')
        prep = read_json(OUT / 'preparation.json')
        with (OUT / 'eligibility.tsv').open() as f:
            rows = list(csv.DictReader(f, delimiter='\t'))
        self.assertEqual(len({(r['role'], r['unit']) for r in rows}), len(rows))
        for role in ['D1', 'D2', 'V1']:
            cohort = [r for r in rows if r['role'] == role]
            passing = []
            for row in cohort:
                eligible = all(int(row[state]) >= cfg['cell_floor'] for state in ['start', 'intermediate', 'destination'])
                self.assertEqual(eligible, row['eligible'] == 'True')
                if eligible:
                    passing.append(row)
            self.assertGreaterEqual(len(passing), cfg['unit_floor'])
            self.assertEqual(len(passing), prep['roles'][role]['units'])
        self.assertEqual([prep['roles'][r]['units'] for r in ['D1', 'D2', 'V1']], [9, 11, 3])
        self.assertTrue(all(not r['unit'].startswith('NC-') for r in rows if r['role'] == 'D1'))
        with (OUT / 'D2_source_crosswalk.tsv').open() as f:
            crosswalk = list(csv.DictReader(f, delimiter='\t'))
        self.assertEqual(len({r['capture'] for r in crosswalk}), 39)
        self.assertEqual(len({r['donor'] for r in crosswalk}), 17)
        self.assertTrue(all(float(r['age_pcw']) == float(r['GEO_age_pcw']) for r in crosswalk))

    def test_discovery_decision_from_saved_unit_effects(self):
        cfg = read_json(BASE / 'config/pilot_v1.json')
        rules = cfg['discovery']
        module = read_json(OUT / 'frozen_programme.json')
        with gzip.open(OUT / 'discovery_gene_effects.tsv.gz', 'rt') as f:
            rows = list(csv.DictReader(f, delimiter='\t'))
        self.assertEqual(len({r['human'] for r in rows}), len(rows))
        self.assertEqual(len({r['mouse'] for r in rows}), len(rows))
        eligible = []
        for row in rows:
            passing = True
            effects = []
            for role in ['D1', 'D2']:
                passing &= float(row[f'{role}_intermediate_detection_fraction']) + 1e-12 >= rules['minimum_detected_unit_fraction']
                for endpoint in ['start', 'destination']:
                    values = [float(v) for k, v in row.items() if k.startswith(role + '_') and k.endswith('_vs_' + endpoint)]
                    median = statistics.median(values)
                    fraction = sum(v > 0 for v in values) / len(values)
                    self.assertAlmostEqual(median, float(row[f'{role}_vs_{endpoint}_median']), places=9)
                    self.assertAlmostEqual(fraction, float(row[f'{role}_vs_{endpoint}_positive_fraction']), places=9)
                    passing &= median >= rules['minimum_median_log2CPM_difference_each_contrast'] and fraction >= rules['minimum_positive_unit_fraction']
                    effects.append(median)
            self.assertEqual(passing, row['four_contrast_pass'] == 'True')
            excluded = row['human'] in cfg['label_exclusions_human'] or row['human'].startswith(('MT-', 'RPL', 'RPS')) or row['human'] == 'XIST'
            passing &= not excluded
            self.assertEqual(passing, row['eligible'] == 'True')
            if passing:
                eligible.append((min(effects), row['human'], row['mouse']))
        eligible.sort(key=lambda r: (-r[0], r[1]))
        self.assertEqual(len(eligible), module['eligible_genes'])
        selected = eligible[:rules['maximum_genes']] if len(eligible) >= rules['minimum_genes'] else []
        self.assertEqual([r[1] for r in selected], module['human_genes'])
        self.assertEqual([r[2] for r in selected], module['mouse_genes'])
        self.assertEqual(len(selected), module['selected_genes'])
        if not selected:
            self.assertEqual(module['status'], 'STOP_NO_QUALIFYING_COMMON_MODULE')
            self.assertFalse((OUT / 'transfer_run.json').exists())
            self.assertEqual([r[1] for r in eligible], module['unselected_eligible_human_genes'])

    def test_original_feasibility_evidence_is_preserved(self):
        old = read_json(BASE / 'execution_record.json')
        self.assertFalse(old['program_scoring_performed'])
        migrated = read_json(BASE / 'history/documentation_migration.json')['files']
        with zipfile.ZipFile(BASE / 'history/feasibility_documentation.zip') as archive:
            for name, record in old['artifacts'].items():
                data = archive.read(name) if name in migrated else (BASE / name).read_bytes()
                self.assertEqual(hashlib.sha256(data).hexdigest(), record['sha256'], name)

    def test_frozen_transfer_requires_both_endpoints(self):
        run = read_json(OUT / 'transfer_run.json')
        for key, path in [('script_sha256', BASE/'scripts/13_transfer_pilot.py'),
                          ('configuration_sha256', BASE/'config/pilot_v1.json'),
                          ('module_sha256', OUT/'frozen_programme.json'),
                          ('execution_plan_sha256', BASE/'TRANSFER_EXECUTION.md')]:
            self.assertEqual(run[key], digest(path))
        for filename, sha in run['outputs'].items():
            self.assertEqual(digest(OUT/filename), sha, filename)
        with gzip.open(OUT/'transfer_cell_scores.tsv.gz', 'rt') as f:
            cells = list(csv.DictReader(f, delimiter='\t'))
        self.assertEqual(len(cells), 8423)
        grouped = defaultdict(list)
        for cell in cells:
            grouped[(cell['role'],cell['unit'],cell['state'])].append(float(cell['programme_score']))
        means = {k: math.fsum(v)/len(v) for k,v in grouped.items()}
        with (OUT/'transfer_summary.tsv').open() as f:
            rows = list(csv.DictReader(f, delimiter='\t'))
        decisions = {}
        for row in rows:
            role, endpoint = row['role'], row['endpoint']
            units = sorted({k[1] for k in means if k[0] == role})
            values = [means[(role,u,'intermediate')]-means[(role,u,endpoint)] for u in units]
            med = statistics.median(values)
            positive = sum(v>0 for v in values)
            minimum_loo = min(statistics.median(values[:i]+values[i+1:]) for i in range(len(values)))
            self.assertAlmostEqual(med, float(row['median_difference']), places=12)
            self.assertEqual(positive, int(row['positive_units']))
            self.assertAlmostEqual(minimum_loo, float(row['minimum_leave_one_out_median']), places=12)
            passed = med>0 and positive/len(values)>=2/3 and minimum_loo>0
            self.assertEqual(passed, row['criterion_met']=='True')
            decisions[(role,endpoint)] = passed
        self.assertTrue(decisions[('V1','start')])
        self.assertFalse(decisions[('V1','destination')])
        self.assertFalse(run['V1_primary_rule_passed'])
        self.assertEqual(run['status'], 'STOP_PRIMARY_TRANSFER_NOT_SUPPORTED')
        self.assertIn('pruned_after_primary_transfer_failure', run['specificity'])
        verification = read_json(OUT/'transfer_verification.json')
        self.assertEqual(verification['verifier_sha256'], digest(BASE/'scripts/17_verify_transfer.py'))
        self.assertEqual(verification['transfer_run_sha256'], digest(OUT/'transfer_run.json'))
        self.assertEqual(verification['raw_score_probes'], 78)
        self.assertLess(verification['maximum_raw_score_error'], 1e-12)

    def test_figures_match_the_inspected_outputs(self):
        render = read_json(OUT/'figure_render.json')
        review = read_json(OUT/'figure_review.json')
        self.assertEqual(render['generator_sha256'], digest(BASE/'scripts/15_plot_pilot.py'))
        self.assertEqual(review['status'], 'visually_checked')
        for name, sha in render['figures'].items():
            self.assertEqual(digest(BASE/'figures/pilot_v1'/name), sha)
        for name, sha in review['png_sha256'].items():
            self.assertEqual(sha, render['figures'][name])


if __name__ == '__main__':
    unittest.main()
