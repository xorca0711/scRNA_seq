"""Portable checks for the separate repair-first exploratory run."""
import csv
import hashlib
import json
import statistics
import unittest
import zipfile
from pathlib import Path

BASE=Path(__file__).resolve().parents[2]/'RQ_Specified/A0_conserved_epithelial_transition_program'
TABLES=BASE/'tables/exploratory'


def read_json(path):
    return json.loads(path.read_text(encoding='utf-8'))


def rows(name):
    with (TABLES/name).open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class A0ExploratoryContract(unittest.TestCase):
    def test_original_run_and_archived_bytes_are_preserved(self):
        publication=read_json(BASE/'exploratory_publication_record.json')
        record=read_json(BASE/'exploratory_execution_record.json')
        self.assertEqual(digest(BASE/'exploratory_execution_record.json'),publication['original_execution_record_sha256'])
        self.assertEqual(digest(BASE/publication['archive']),publication['archive_sha256'])
        self.assertEqual(record['validation'],'PASS')
        self.assertEqual(read_json(BASE/'exploratory_validation.json')['status'],'PASS')
        with zipfile.ZipFile(BASE/publication['archive']) as archive:
            self.assertEqual(set(archive.namelist()),set(publication['archived_original_artifacts']))
            for name,info in record['artifacts'].items():
                if name in publication['archived_original_artifacts']:
                    actual=hashlib.sha256(archive.read(name)).hexdigest()
                else:actual=digest(BASE/name)
                self.assertEqual(actual,info['sha256'],name)
        self.assertFalse(publication['analytical_estimator_or_result_changes'])

    def test_primary_and_variant_remain_frozen(self):
        cfg=read_json(BASE/'exploratory_config.json')
        primary=read_json(BASE/'frozen_repair_program.json')
        variant=read_json(BASE/'frozen_control_variant.json')
        self.assertEqual(primary['config_sha256'],digest(BASE/'exploratory_config.json'))
        self.assertEqual(read_json(BASE/'developmental_expression_manifest.json')['frozen_program_sha256'],digest(BASE/'frozen_repair_program.json'))
        self.assertEqual(len(primary['genes']),len(set(primary['genes'])))
        self.assertEqual(len(primary['genes']),50)
        generic=set().union(*(set(cfg['control_modules'][n]) for n in cfg['generic_control_names']))
        self.assertEqual(variant['genes'],[g for g in primary['genes'] if g not in generic])
        self.assertEqual(len(variant['genes']),31)
        self.assertFalse(primary['transfer_outcomes_seen'])
        self.assertFalse(variant['transfer_outcomes_seen'])

    def test_repair_heldout_summary_and_depth_floors(self):
        data=rows('repair_heldout_effects.csv');result=read_json(BASE/'stage_E1_result.json')
        for endpoint in ['AT2','AT1']:
            group=[r for r in data if r['endpoint']==endpoint]
            self.assertEqual(len(group),9)
            values=[float(r['difference']) for r in group]
            self.assertEqual(sum(v>0 for v in values),result['heldout_results'][endpoint]['positive_mice'])
            self.assertAlmostEqual(statistics.median(values),result['heldout_results'][endpoint]['median_difference'],places=12)
        for r in rows('repair_depth_matched_effects.csv'):
            self.assertGreaterEqual(int(r['n_intermediate']),30)
            self.assertEqual(r['n_start'],r['n_intermediate'])
            self.assertEqual(r['n_end'],r['n_intermediate'])
            self.assertGreater(float(r['difference']),0)

    def test_transfer_arithmetic_and_biological_unit_limits(self):
        result=read_json(BASE/'stage_E3_result.json')
        for context,qualified in [('development',1),('intestine',2)]:
            data=rows(context+'_transfer_effects.csv')
            for r in data:
                self.assertAlmostEqual(float(r['difference']),float(r['intermediate_mean'])-float(r['endpoint_mean']),places=12)
                passes=min(int(r[c]) for c in ['n_start','n_intermediate','n_end'])>=30
                self.assertEqual(passes,r['all_states_at_least_30']=='True')
            primary=[r for r in data if r['module']=='candidate']
            self.assertEqual(len({r['unit'] for r in primary if r['all_states_at_least_30']=='True'}),qualified)
            for endpoint,expected in result['results'][context]['effects'].items():
                q=[r for r in primary if r['endpoint']==endpoint and r['all_states_at_least_30']=='True']
                self.assertEqual(sum(float(r['difference'])>0 for r in q),expected['positive_qualified_groups'])
        self.assertEqual(result['input_audits']['development']['independent_animals'],'not_recovered')

    def test_sensitivity_does_not_manufacture_replication(self):
        data=rows('intestine_depth_matched_sensitivity.csv')
        qualifying={r['unit'] for r in data if r['all_states_at_least_30']=='True'}
        self.assertEqual(qualifying,{'Control-Mouse3'})
        m1=next(r for r in data if r['unit']=='Control-Mouse1' and r['module']=='candidate_without_generic_controls' and r['endpoint']=='Stem')
        self.assertLess(float(m1['median_difference']),0)
        for r in data:
            self.assertLessEqual(float(r['draw_q05']),float(r['median_difference']))
            self.assertLessEqual(float(r['median_difference']),float(r['draw_q95']))
        self.assertEqual(read_json(BASE/'stage_E3a_result.json')['qualified_after_matching'],1)

    def test_status_writes_are_separate_from_other_a0_runs(self):
        text=(BASE/'scripts/verify_exploratory.py').read_text(encoding='utf-8')
        self.assertNotIn("BASE/'readiness.json'",text)
        self.assertNotIn("BASE/'decisions.json'",text)
        self.assertIn("BASE/'exploratory_readiness.json'",text)
        self.assertIn("BASE/'exploratory_decisions.json'",text)
        self.assertTrue(read_json(BASE/'exploratory_readiness.json')['exploratory_pilot_complete'])
        self.assertFalse(read_json(BASE/'readiness.json')['expression_program_learned'])
        self.assertTrue((BASE/'tables/pilot_v1/pilot_status.json').exists())


if __name__=='__main__':unittest.main()
