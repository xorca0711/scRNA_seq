"""Independent dense/CSC parity fixture, including empty cells and QC exclusion."""
from pathlib import Path
import tempfile,gzip,json
import numpy as np
from u5_human_count_store import build_store,aggregate_store,read_panel_store

def main():
    root=Path(__file__).resolve().parents[1];out=root/'trials/u5_human_full'
    with tempfile.TemporaryDirectory(dir=root/'cache') as temp:
        p=Path(temp);source=p/'example.txt.gz';target=p/'counts.h5'
        matrix=np.array([[0,1,2,0,3],[0,0,0,0,0],[0,4,0,0,5]],dtype=np.int64)
        with gzip.open(source,'wt') as h:
            h.write('a\tb\tc\td\te\n')
            for g,row in zip(['GENE_A','GENE_B','MT-X'],matrix):h.write(g+'\t'+'\t'.join(map(str,row))+'\n')
        panel,metadata=build_store(source,target,['MT-X','MISSING','GENE_A'])
        assert np.array_equal(panel.toarray(),np.column_stack([matrix[2],np.zeros(5),matrix[0]]))
        assert np.array_equal(metadata['full_library_size'],matrix.sum(axis=0))
        assert np.array_equal(metadata['n_genes_by_counts'],(matrix>0).sum(axis=0))
        assert np.array_equal(metadata['full_gene_max_count'],matrix.max(axis=0))
        recovered,recovered_meta=read_panel_store(target,['MT-X','MISSING','GENE_A'])
        assert np.array_equal(recovered.toarray(),panel.toarray())
        assert np.array_equal(recovered_meta['full_library_size'],metadata['full_library_size'])
        codes=np.array([-1,0,1,-1,0]);genes,counts,detected=aggregate_store(target,codes,2,block_genes=2)
        assert list(genes)==['GENE_A','GENE_B','MT-X']
        assert np.array_equal(counts,np.vstack([matrix[:,[1,4]].sum(axis=1),matrix[:,2]]))
        assert np.array_equal(detected,np.vstack([(matrix[:,[1,4]]>0).sum(axis=1),matrix[:,2]>0]))
        with gzip.open(p/'invalid.gz','wt') as h:h.write('a\nbad\t1.5\n')
        try:build_store(p/'invalid.gz',p/'bad.h5',['bad'])
        except ValueError:pass
        else:raise AssertionError('Fractional source silently converted to integer')
    report=dict(status='passed',checks=['selected panel order and zero-filled off-panel feature','full library and gene detection totals including empty cells','blocked CSC pseudobulk exact dense parity','detected-cell counts exact parity','excluded cells omitted','fractional count rejection'])
    (out/'count_store_validation.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))

if __name__=='__main__':main()
