"""Hand-calculated source weighting, partition and transcript-coordinate cases."""
import gzip
import importlib.util
from pathlib import Path
import tempfile
import unittest

PATH = Path(__file__).resolve().parents[2]/'RQ_Specified/A1_transitional_epithelial_state_distinction/scripts/a1_robustness.py'
SPEC = importlib.util.spec_from_file_location('a1_robustness_tested', PATH)
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)


class RobustnessTests(unittest.TestCase):
    def test_weighting_preserves_small_and_zero_sources(self):
        equal, pooled = M.weighted_summary([(0,1),(8,9)])
        self.assertAlmostEqual(equal,4/9)
        self.assertEqual(pooled,.8)
        self.assertEqual(M.weighted_summary([(0,1)]),(0,0))
        self.assertEqual(M.weighted_summary([(8,9)]),(8/9,8/9))

    def test_invalid_counts_do_not_become_zero(self):
        for counts in [[],[(0,0)],[(2,1)],[(-1,4)]]:
            with self.assertRaises(ValueError):M.weighted_summary(counts)

    def test_partition_agreement_does_not_require_a_label_mapping(self):
        self.assertEqual(M.adjusted_rand(['a','a','b','b'],['2','2','1','1']),1)
        self.assertAlmostEqual(M.adjusted_rand(['a','a','b','b'],['x','y','x','y']),-.5)
        self.assertEqual(M.adjusted_rand(['a']*4,['x','x','y','y']),0)
        self.assertEqual(M.adjusted_rand(['a']*4,['b']*4),1)

    def test_native_transcript_coordinates_parents_and_duplicates(self):
        loci=[dict(gene='P',gene_group='test',chrom='chr1',strand='+',gene_start=10,gene_end=100,tss=10),
              dict(gene='M',gene_group='test',chrom='chr2',strand='-',gene_start=200,gene_end=300,tss=299)]
        records=[('chr1','gene',11,100,'+','ID=gP;gene=P'),
                 ('chr2','gene',201,300,'-','ID=gM;gene=M'),
                 ('chr1','transcript',11,100,'+','ID=p1;Parent=gP;gene=P'),
                 ('chr1','transcript',21,90,'+','ID=p2;Parent=gP;gene=P'),
                 ('chr1','mRNA',21,99,'+','ID=p3;Parent=gP;gene=P'),
                 ('chr1','transcript',31,90,'+','ID=orphan;Parent=other;gene=P'),
                 ('chr1','transcript',1,90,'+','ID=outside;Parent=gP;gene=P'),
                 ('chr2','transcript',211,250,'-','ID=m1;Parent=gM;gene=M')]
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'genes.gff.gz'
            with gzip.open(p,'wt') as h:
                for chrom,kind,left,right,strand,attrs in records:
                    h.write(f'{chrom}\tref\t{kind}\t{left}\t{right}\t.\t{strand}\t.\t{attrs}\n')
            rows,audit=M.transcript_tss(p,loci)
        self.assertEqual({(r['gene'],r['tss']) for r in rows},{('P',10),('P',20),('M',249),('M',299)})
        self.assertEqual(next(r for r in rows if r['gene']=='P' and r['tss']==20)['transcript_ids'],'p2;p3')
        self.assertEqual(sum(r['is_baseline'] for r in rows),2)
        self.assertEqual({r['transcript_id'] for r in audit if not r['eligible']},{'orphan','outside'})


if __name__=='__main__':unittest.main()
