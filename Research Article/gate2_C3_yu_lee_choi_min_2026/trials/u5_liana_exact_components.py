"""Exact sufficient input for LIANA's CellChat-like magnitude, with QA.

LIANA uses selected resource genes plus the global maximum of log-normalized
expression. A named computational auxiliary column stores each cell's observed
maximum across ALL assayed genes. It is not a gene and never enters an edge.
It also preserves nonzero-cell retention when all resource genes are zero.
No renormalization is performed after compression.
"""
from pathlib import Path
import json,sys
PAPER=Path(__file__).resolve().parents[1];ROOT=PAPER.parents[1]


def normalized_components(a):
    import numpy as np,pandas as pd,anndata as ad
    from scipy import sparse
    denominator=a.obs.full_library_size.to_numpy(dtype=np.float32)/np.float32(10000)
    nonzero=denominator>0
    a=a[nonzero].copy();denominator=denominator[nonzero]
    x=a.X.tocsr().astype(np.float32)
    x.data/=np.repeat(denominator,np.diff(x.indptr));np.log1p(x.data,out=x.data)
    maximum=np.log1p(a.obs.full_gene_max_count.to_numpy(dtype=np.float32)/denominator)
    name='__observed_full_gene_maximum_not_a_gene__'
    assert name not in a.var_names
    matrix=sparse.hstack([x,sparse.csr_matrix(maximum[:,None])],format='csr',dtype=np.float32)
    z=ad.AnnData(matrix,obs=a.obs.copy(),var=pd.DataFrame(index=list(a.var_names)+[name]))
    z.uns['compression_note']=__doc__
    z.uns['log1p']={'base':None}
    return z


def validate():
    import numpy as np,pandas as pd,anndata as ad,scanpy as sc
    from liana.method import cellchat
    cache=PAPER/'cache/u5_ipf_liana/sampled_triad.h5ad'
    ref=ad.read_h5ad(cache,backed='r')
    out=PAPER/'trials/u5_liana_robustness';out.mkdir(exist_ok=True)
    pairs=pd.read_csv(PAPER/'trials/u5_ipf_liana/fixed_donor_pairs.csv')
    metadata=ref.obs[['donor','disease']].drop_duplicates()
    donors=[metadata[metadata.disease==d].donor.iloc[0] for d in ['IPF','control']]
    checks=[]
    for donor in donors:
        full=ref[ref.obs.donor.eq(donor)].to_memory()
        for resource in ['consensus','cellchatdb']:
            table=pd.read_csv(PAPER/'trials/u5_ipf_liana'/('resource_'+resource+'.csv'))
            genes=sorted({g for col in ['ligand','receptor'] for entry in table[col] for g in entry.split('_')}&set(full.var_names))
            part=full[:,genes].copy()
            part.obs['full_library_size']=np.asarray(full.X.sum(axis=1)).ravel()
            part.obs['full_gene_max_count']=full.X.max(axis=1).toarray().ravel()
            compressed=normalized_components(part)
            baseline=full.copy();sc.pp.normalize_total(baseline,target_sum=10000);sc.pp.log1p(baseline)
            pair=pairs[pairs.donor==donor][['source','target']]
            kw=dict(groupby='celltype',resource=table[['ligand','receptor']],groupby_pairs=pair,expr_prop=.1,min_cells=50,use_raw=False,n_perms=None,verbose=False,inplace=False)
            x=cellchat(baseline,**kw);y=cellchat(compressed,**kw)
            keys=['source','target','ligand_complex','receptor_complex']
            x=x.set_index(keys).sort_index();y=y.set_index(keys).sort_index()
            assert x.index.equals(y.index)
            for col in ['ligand_props','receptor_props','ligand_trimean','receptor_trimean','mat_max','lr_probs']:
                assert np.allclose(x[col],y[col],rtol=2e-6,atol=2e-7),col
            checks.append({'donor':donor,'resource':resource,'rows':len(x),'max_absolute_score_difference':float(np.max(abs(x.lr_probs-y.lr_probs)))})
    ref.file.close()
    report={'status':'passed','scope':'full-gene normalization and native LIANA versus sufficient-component input, two donors and both resources','tolerance_rtol':2e-6,'tolerance_atol':2e-7,'checks':checks}
    (out/'compression_validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report))


if __name__=='__main__':validate()
