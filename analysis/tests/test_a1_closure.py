"""Hand-calculated failure and correspondence cases, without assay downloads."""
import importlib.util
from pathlib import Path
import unittest

PATH = Path(__file__).resolve().parents[2]/'RQ_Specified/A1_transitional_epithelial_state_distinction/scripts/a1_closure.py'
SPEC = importlib.util.spec_from_file_location('a1_closure',PATH)
MODULE = importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(MODULE)


class ClosureTests(unittest.TestCase):
    def test_source_identity_uses_number_and_gate_not_order(self):
        self.assertEqual(MODULE.cd44_key('R26_7_CD44+','column'),
                         MODULE.cd44_key('HTS053_18_R267_CD44plus_S18_R1_001.fastq.gz','fastq'))
        self.assertNotEqual(MODULE.cd44_key('OG6_CD44-','column'),MODULE.cd44_key('OG7_CD44-','column'))
        self.assertNotEqual(MODULE.cd44_key('R26_1_CD44+','column'),MODULE.cd44_key('R26_1_CD44-','column'))
        with self.assertRaises(ValueError):MODULE.cd44_key('arbitrary_R261_CD44plus.fastq.gz','fastq')

    def test_abstention_is_separate_from_reclassification(self):
        rows=[{'newleiden':'0','cell type':'AT2','clusterK12':'01','clusterK12_stringent':'other'},
              {'newleiden':'1','cell type':'HPCS','clusterK12':'02','clusterK12_stringent':'02'},
              {'newleiden':'2','cell type':'bad','clusterK12':'03','clusterK12_stringent':'04'}]
        observed=MODULE.annotation_audit(rows,{'0':'AT2','1':'HPCS'})
        self.assertEqual(observed,dict(cells=3,unmapped_newleiden=1,cell_type_mismatches=1,
                                      stringent_other=1,literal_changes=2,retained_code_changes=1))
        with self.assertRaises(ValueError):MODULE.annotation_audit([], {})


if __name__=='__main__':unittest.main()
