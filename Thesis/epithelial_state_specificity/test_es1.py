"""Focused checks for the scientific estimand and joint depth sampler."""
import json
import unittest
from pathlib import Path
import numpy as np
from run_es1 import joint_sample

class SamplingChecks(unittest.TestCase):
    def test_full_budget_recovers_counts(self):
        counts=np.array([[2,1],[3,5],[5,4]],dtype=np.int32)
        np.testing.assert_array_equal(joint_sample(counts,np.array([10,10]),10,17),counts)

    def test_joint_budget_and_residual(self):
        counts=np.array([[8,0,4],[2,9,3],[5,1,1]],dtype=np.int32)
        total=np.array([20,12,9])
        for seed in range(10):
            sampled=joint_sample(counts,total,7,seed)
            self.assertTrue((sampled<=counts).all())
            self.assertTrue((sampled.sum(axis=0)<=7).all())
            self.assertTrue((7-sampled.sum(axis=0)<=total-counts.sum(axis=0)).all())
        np.testing.assert_array_equal(joint_sample(counts,total,7,17),joint_sample(counts,total,7,17))

    def test_hypergeometric_expectation(self):
        # A large bank of identical cells checks marginal behavior without
        # pretending they are biological replicates.
        counts=np.tile(np.array([[3],[7],[5]],dtype=np.int32),(1,20000))
        sampled=joint_sample(counts,np.full(20000,20),8,17)
        np.testing.assert_allclose(sampled.mean(axis=1),np.array([3,7,5])*8/20,atol=.035)

    def test_source_modules_and_label_exclusion(self):
        cfg=json.loads((Path(__file__).parent/'modules.json').read_text(encoding='utf-8'))
        modules={m['name']:m for m in cfg['modules']}
        for name in ['ADI','AT2','AT1']:
            self.assertEqual(len(modules[f'{name}_published_400']['genes']),400)
            self.assertFalse(set(modules[f'{name}_published_holdout']['genes']) & {'Krt8','Cldn4','Sftpc','Cebpa'})
        self.assertNotIn('Cebpa',modules['AT2_common_panel']['genes'])

if __name__=='__main__':unittest.main()
