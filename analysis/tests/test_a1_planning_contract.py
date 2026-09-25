"""A title-derived pair and an overlapping animal pool must not pass inference."""
import copy
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[2]
SCRIPTS=ROOT/'RQ_Specified/A1_transitional_epithelial_state_distinction/scripts'
sys.path.insert(0,str(SCRIPTS))
from a1_contract import digest,readiness
spec=importlib.util.spec_from_file_location('a1_runner',SCRIPTS/'03_run_paired_counts.py')
runner=importlib.util.module_from_spec(spec);spec.loader.exec_module(runner)


class A1PlanningTests(unittest.TestCase):
    def fixture(self,root):
        path=root/'counts.tsv';path.write_text('feature\ts0\ts1\ts2\ts3\ts4\ts5\ng1\t1\t2\t3\t4\t5\t6\n')
        samples=[dict(sample_id=f's{i}',accession='GSE1',assay='bulk_ATAC',genotype='WT',
                     biological_unit_id=f'animal{i//2}',member_ids=[f'animal{i//2}'],unit_kind='animal',
                     identity_verified=True,identity_evidence='verified-source-table',counts_column=f's{i}',
                     group='positive' if i%2 else 'negative') for i in range(6)]
        contrast=dict(accession='GSE1',assay='bulk_ATAC',sample_ids=[s['sample_id'] for s in samples],
                      reference='negative',case='positive',minimum_independent_pairs=3,
                      counts_path='counts.tsv',counts_sha256=digest(path),input_scale='raw_integer_counts',
                      assay_qc_verified=True,assay_qc_evidence='verified-qc-report')
        return contrast,samples

    def test_candidate_titles_cannot_supply_biological_replication(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);c,samples=self.fixture(root)
            for s in samples:
                s['candidate_pair_id']=s['biological_unit_id']
                s['biological_unit_id']=None
                s['identity_verified']=False
            holds,_=readiness(c,samples,root)
            self.assertIn('Too few verified independent pairs',holds)

    def test_pairing_pool_overlap_and_payload_drift_fail(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);c,s=self.fixture(root)
            self.assertEqual(readiness(c,s,root)[0],[])
            overlap=copy.deepcopy(s)
            overlap[2]['member_ids']=overlap[3]['member_ids']=['animal0']
            self.assertTrue(any('overlapping' in x for x in readiness(c,overlap,root)[0]))
            bad=copy.deepcopy(s);bad[1]['group']='negative'
            self.assertTrue(any('incomplete/repeated' in x for x in readiness(c,bad,root)[0]))
            (root/'counts.tsv').write_text('tampered')
            self.assertIn('Count file missing or hash mismatch',readiness(c,s,root)[0])

    def test_counts_align_by_header_and_reject_normalized_values(self):
        with tempfile.TemporaryDirectory() as directory:
            p=Path(directory)/'counts.tsv'
            p.write_text('feature\tb\ta\ng1\t4\t8\n')
            selected=[{'counts_column':'a'},{'counts_column':'b'}]
            self.assertEqual(runner.read_counts(p,selected),[['g1',8,4]])
            p.write_text('feature\tb\ta\ng1\t4.2\t8\n')
            with self.assertRaises(ValueError):runner.read_counts(p,selected)
