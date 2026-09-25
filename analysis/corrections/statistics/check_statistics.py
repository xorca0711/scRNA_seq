"""Focused checks for identity, reference parity and inferential bookkeeping."""
from __future__ import annotations
import hashlib
import json
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT))
from analysis.lib.repository_paths import recorded_file
site = ROOT / '.venv-x64/Lib/site-packages'
if site.is_dir():
    sys.path.insert(0, str(site))
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control


class StatisticalCorrectionChecks(unittest.TestCase):
    def test_historical_count_parity(self):
        for name in ['mouse_myeloid','GSE136831','GSE135893']:
            record = json.loads((HERE / f'{name}_preparation.json').read_text())
            self.assertTrue(record['historical_cell_and_library_total_parity'])
            arr = np.load(HERE / 'cache' / f'{name}.npz')
            self.assertTrue(np.issubdtype(arr['counts'].dtype, np.integer))
            self.assertGreaterEqual(arr['counts'].min(), 0)
            self.assertEqual(int(arr['counts'].sum()), record['total_counts'])
            self.assertEqual(len(set(arr['genes'])), len(arr['genes']))

    def test_w1_independent_reference_residuals(self):
        old = pd.read_csv(ROOT / 'Research Article/gate1_01_niethamer_2025/trials/w1_amac_pseudobulk_de/w1_sets.csv')
        new = pd.read_csv(HERE / 'tables/w1_reference_camera.csv')
        joined = old[old.tested].merge(new[new.method=='historical_logCPM'],on=['contrast','set'],validate='one_to_one')
        self.assertEqual(len(joined), 10)
        # Independent official limma QR residuals agree with the historical
        # Python residual correlations: a strong check of input/design alignment.
        np.testing.assert_allclose(joined.inter_gene_corr,joined.Correlation,atol=1e-12,rtol=0)
        self.assertEqual(set(joined.loc[joined.contrast=='A','residual_df']), {7})
        self.assertEqual(set(joined.loc[joined.contrast=='B','residual_df']), {6})
        self.assertEqual(set(zip(joined.loc[joined.contrast=='A','n_first'],joined.loc[joined.contrast=='A','n_second'])),{(4,8)})
        self.assertEqual(set(zip(joined.loc[joined.contrast=='B','n_first'],joined.loc[joined.contrast=='B','n_second'])),{(8,3)})
        np.testing.assert_array_equal(joined['size'].to_numpy(),joined.NGenes.to_numpy())

    def test_multiplicity_uses_complete_discovery(self):
        frozen = pd.read_csv(ROOT / 'Research Article/gate1_01_niethamer_2025/trials/g2_gsea_ipf/g2_replication.csv')
        for cohort in ['GSE136831','GSE135893']:
            allsets = pd.read_csv(HERE / 'cache' / f'{cohort}_camera_allsets.csv.gz')
            self.assertTrue(allsets.PValue.between(0,1).all())
            if cohort=='GSE136831':
                self.assertGreater(len(allsets), len(frozen))
                expected = false_discovery_control(allsets.PValue.to_numpy())
            else:
                padded = np.r_[allsets.PValue.to_numpy(),np.ones(len(frozen)-len(allsets))]
                expected = false_discovery_control(padded)[:len(allsets)]
            np.testing.assert_allclose(expected,allsets.FDR_global,atol=1e-12,rtol=0)
            candidates = pd.read_csv(HERE / 'tables' / f'{cohort}_camera_candidates.csv')
            self.assertEqual(set(zip(candidates.compartment,candidates['set'])),set(zip(frozen.compartment,frozen['set'])))
            self.assertTrue(candidates.loc[~candidates.eligible,'PValue'].isna().all())

    def test_age_identifiability_and_arm_floor(self):
        a = pd.read_csv(HERE / 'tables/mouse_design_audit.csv').set_index('label')
        self.assertLess(a.loc['G1_age_time_not_identifiable','rank'],a.loc['G1_age_time_not_identifiable','columns'])
        self.assertEqual(a.loc['G1_same_round_male_only_below_arm_floor','n_second'],1)

    def test_correlation_sensitivity_keeps_the_same_input_design(self):
        for cohort in ['GSE136831','GSE135893']:
            a = pd.read_csv(HERE / 'tables' / f'{cohort}_camera_candidates.csv')
            b = pd.read_csv(HERE / 'tables' / f'{cohort}_fixed001_camera_candidates.csv')
            m = a.merge(b,on=['compartment','set'],suffixes=('_primary','_fixed'),validate='one_to_one')
            for field in ['eligible','NGenes','genes_tested','n_first','n_second','residual_df']:
                np.testing.assert_array_equal(m[field+'_primary'].to_numpy(),m[field+'_fixed'].to_numpy())
            self.assertEqual(set(b.loc[b.eligible,'Correlation']),{.01})

    def test_lodo_has_one_omission_per_donor(self):
        for cohort in ['GSE136831','GSE135893']:
            units = pd.read_csv(HERE / 'cache' / f'{cohort}_units.csv')
            lod = pd.read_csv(HERE / 'tables' / f'{cohort}_lodo_summary.csv')
            sizes = units.groupby('comp').size()
            self.assertTrue((lod.omissions.to_numpy()==lod.compartment.map(sizes).to_numpy()).all())
            self.assertTrue(lod.direction_reversals.between(0,lod.omissions).all())
            detail = pd.read_csv(HERE / 'cache' / f'{cohort}_lodo_all.csv.gz')
            self.assertEqual(len(detail[['compartment','omitted_donor']].drop_duplicates()),len(units))
            self.assertFalse(detail.duplicated(['compartment','set','omitted_donor']).any())

    def test_transcript_shares_sum_to_one(self):
        d = pd.read_csv(HERE / 'tables/subtype_transcript_contributions.csv')
        sums = d.groupby(['cohort','compartment','set'])[['mean_set_transcript_share_IPF','mean_set_transcript_share_control']].sum()
        np.testing.assert_allclose(sums.to_numpy(),1,atol=1e-10,rtol=0)

    def test_recorded_files_match_hashes(self):
        record = json.loads((HERE / 'run_record.json').read_text())
        for section in ['historical_inputs','gene_set_inputs','code','outputs']:
            for entry in record[section]:
                self.assertEqual(hashlib.sha256(recorded_file(ROOT, entry['path'], entry['sha256']).read_bytes()).hexdigest(),entry['sha256'],entry['path'])

    def test_gene_sets_match_frozen_libraries(self):
        for trial in ['g1_gsea_by_phase','g2_gsea_ipf']:
            tag = trial.split('_')[0]
            record = json.loads((ROOT / 'Research Article/gate1_01_niethamer_2025/trials' / trial / f'{tag}_run_record.json').read_text())
            for facts in record['results']['gene_set_files'].values():
                self.assertEqual(hashlib.sha256(recorded_file(ROOT, facts['file'], facts['sha256']).read_bytes()).hexdigest(),facts['sha256'])


if __name__ == '__main__':
    unittest.main(verbosity=1)
