"""Full-cell native LIANA and one-at-a-time sensitivity checks.

All deposited IPF/control cells enter source profiles; LR uses fixed annotated
subtype pairs in the same donor. Missing edges remain missing, never zero.
"""
from pathlib import Path
import os,sys,json,time,warnings,hashlib
from datetime import datetime,timezone
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ.setdefault(k,'2')
PAPER=Path(__file__).resolve().parents[1];ROOT=PAPER.parents[1]
sys.path.insert(0,str(ROOT))
from analysis.lib.provenance import code_identity,sha256_file,write_json_atomic
from u5_liana_exact_components import normalized_components

LABELS={
 'GSE136831':{'ATII':'epithelial','Fibroblast':'fibroblast','Myofibroblast':'fibroblast','Macrophage':'myeloid','Macrophage_Alveolar':'myeloid'},
 'GSE135893':{'AT2':'epithelial','Fibroblasts':'fibroblast','HAS1 High Fibroblasts':'fibroblast','PLIN2+ Fibroblasts':'fibroblast','Myofibroblasts':'fibroblast','Macrophages':'myeloid','Proliferating Macrophages':'myeloid'}}
DIRECTIONS={('myeloid','fibroblast'),('myeloid','epithelial'),('epithelial','fibroblast'),('fibroblast','myeloid'),('fibroblast','epithelial')}
CONFIGS=[dict(name='all_cells_primary',floor=50,expr=.1,cap=None,seed=20260924),
 dict(name='cap500_seed1',floor=50,expr=.1,cap=500,seed=20260924),
 dict(name='cap500_seed2',floor=50,expr=.1,cap=500,seed=20260925),
 dict(name='floor30',floor=30,expr=.1,cap=None,seed=20260924),
 dict(name='floor100',floor=100,expr=.1,cap=None,seed=20260924),
 dict(name='expression005',floor=50,expr=.05,cap=None,seed=20260924),
 dict(name='expression020',floor=50,expr=.2,cap=None,seed=20260924)]


def profile(a,out):
    import numpy as np,pandas as pd
    rows=[]
    for (donor,disease,label),indices in a.obs.groupby(['donor','disease','label'],observed=True).indices.items():
        if disease not in ['IPF','control']:continue
        x=a.X[indices].tocsr();lib=a.obs.iloc[indices].full_library_size.to_numpy()
        norm=x.astype(float).multiply((10000/np.maximum(lib,1))[:,None]).tocsr()
        counts=np.asarray(x.sum(axis=0)).ravel();detect=np.asarray((x>0).sum(axis=0)).ravel()
        means=np.asarray(norm.mean(axis=0)).ravel()
        for j,g in enumerate(a.var_names):rows.append(dict(donor=donor,disease=disease,label=label,gene=g,cells=len(indices),detected_cells=int(detect[j]),detection_fraction=detect[j]/len(indices),count_sum=int(counts[j]),full_library_sum=int(lib.sum()),mean_normalized_10000=means[j]))
    pd.DataFrame(rows).to_csv(out/'all_label_source_recipient_profiles.csv.gz',index=False,compression='gzip')
    # Cell fractions are within recovered deposited cells, not tissue abundance.
    c=a.obs[a.obs.disease.isin(['IPF','control'])].groupby(['donor','disease','label'],observed=True).size().rename('cells').reset_index()
    c['recovered_donor_cells']=c.groupby('donor').cells.transform('sum')
    c['recovered_cell_fraction']=c.cells/c.recovered_donor_cells
    c.to_csv(out/'recovered_cell_composition.csv',index=False)


def main():
    import numpy as np,pandas as pd,anndata as ad
    from liana.method import cellchat
    warnings.filterwarnings('ignore',category=FutureWarning,module='liana')
    warnings.filterwarnings('ignore',category=ad.ImplicitModificationWarning)
    cohort=sys.argv[1];out=PAPER/'trials/u5_liana_robustness'/cohort;out.mkdir(parents=True,exist_ok=True)
    cache=PAPER/'cache/u5_liana_robustness'/cohort;cache.mkdir(parents=True,exist_ok=True)
    qa=json.loads((out.parent/'compression_validation.json').read_text());assert qa['status']=='passed'
    src=PAPER/'cache/u5_full_source_panels'/(cohort+'.h5ad')
    assert json.loads((PAPER/'trials/u5_full_source_panels'/cohort/'run_record.json').read_text())['status']=='completed'
    if (out/'run_record.json').exists():raise RuntimeError('Existing run; inspect and resume explicitly')
    resources=pd.read_csv(PAPER/'trials/u5_full_source_panels'/cohort/'resource_edges.csv')
    # Match the original frozen LR family set. WNT is kept separately as exploratory.
    resources['family']=resources.ligand.map(lambda x:'core' if any(g in ['IL1A','IL1B','AREG','HBEGF','TGFB1'] for g in x.split('_')) else 'exploratory')
    resources.to_csv(out/'frozen_resource_edges.csv',index=False)
    spec=dict(frozen_utc=datetime.now(timezone.utc).isoformat(),cohort=cohort,configs=CONFIGS,labels=LABELS[cohort],directions=sorted(DIRECTIONS),method='native LIANA CellChat-like; n_perms=None',input_sha256=sha256_file(src),compression_validation_sha256=sha256_file(out.parent/'compression_validation.json'),sampling='independent stable SHA256-derived donor/subtype seeds; cap sensitivity, not original cap replicate',inference='descriptive, missing rows not zero; resource WNT edges exploratory',source_profiles='all deposited annotated IPF/control cells, full-library normalization; detection and recovered abundance reported separately')
    write_json_atomic(out/'specification.json',spec)
    state=dict(status='running',started_utc=spec['frozen_utc'],pid=os.getpid(),code=code_identity(ROOT,__file__),completed_calls=0);write_json_atomic(out/'run_record.json',state);start=time.monotonic()
    try:
        a=ad.read_h5ad(src);profile(a,out)
        a.obs['celltype']=a.obs.label.astype(str);a.obs['compartment']=a.obs.celltype.map(LABELS[cohort])
        a=a[a.obs.compartment.notna() & a.obs.disease.isin(['IPF','control'])].copy()
        coverage=a.obs.groupby(['donor','disease','celltype','compartment'],observed=True).size().rename('cells').reset_index()
        coverage.to_csv(out/'donor_subtype_coverage.csv',index=False)
        scored_all=[];eligibility=[];pair_rows=[]
        for config in CONFIGS:
            good=coverage[coverage.cells>=config['floor']]
            for donor,groups in good.groupby('donor',sort=True,observed=True):
                pairs=[dict(source=r.celltype,target=s.celltype) for r in groups.itertuples() for s in groups.itertuples() if (r.compartment,s.compartment) in DIRECTIONS]
                if not pairs:continue
                pairs=pd.DataFrame(pairs).drop_duplicates()
                disease=groups.disease.iloc[0]
                assert groups.disease.nunique()==1
                selected=[]
                for label in groups.celltype:
                    indices=np.flatnonzero(a.obs.donor.eq(donor).to_numpy() & a.obs.celltype.eq(label).to_numpy())
                    if config['cap'] and len(indices)>config['cap']:
                        seed=int.from_bytes(hashlib.sha256(f"{config['seed']}|{donor}|{label}".encode()).digest()[:8],'little')
                        indices=np.random.default_rng(seed).choice(indices,config['cap'],replace=False)
                    selected.extend(indices)
                part=normalized_components(a[np.sort(selected)])
                for r in pairs.itertuples(index=False):pair_rows.append(dict(config=config['name'],donor=donor,disease=disease,source=r.source,target=r.target))
                for resource in ['consensus','cellchatdb']:
                    table=resources[resources.resource.eq(resource)]
                    z=cellchat(part,groupby='celltype',resource=table[['ligand','receptor']],groupby_pairs=pairs,expr_prop=config['expr'],min_cells=config['floor'],use_raw=False,n_perms=None,verbose=False,inplace=False)
                    assert not any('pval' in c.lower() for c in z.columns)
                    z=z.merge(table[['ligand','receptor','family']],left_on=['ligand_complex','receptor_complex'],right_on=['ligand','receptor'],validate='many_to_one')
                    for k,v in dict(config=config['name'],donor=donor,disease=disease,resource=resource,cohort=cohort).items():z[k]=v
                    z.to_csv(cache/f"{config['name']}_{resource}_{donor}.csv",index=False)
                    scored_all.append(z)
                    eligibility.append(dict(config=config['name'],donor=donor,disease=disease,resource=resource,cells=part.n_obs,pairs=len(pairs),retained_edges=len(z)))
                    state['completed_calls']+=1
                state.update(config=config['name'],donor=donor,elapsed_seconds=round(time.monotonic()-start,1));write_json_atomic(out/'run_record.json',state)
            print(cohort,config['name'],'done',state['completed_calls'],flush=True)
        scores=pd.concat(scored_all,ignore_index=True);scores.to_csv(out/'donor_lr_scores.csv.gz',index=False,compression='gzip')
        pd.DataFrame(eligibility).to_csv(out/'call_coverage.csv',index=False);pd.DataFrame(pair_rows).to_csv(out/'fixed_pairs_by_config.csv',index=False)
        keys=['config','donor','resource','source','target','ligand_complex','receptor_complex']
        assert not scores.duplicated(keys).any()
        assert np.isfinite(scores.lr_probs).all() and scores.lr_probs.between(0,1).all()
        for config in CONFIGS:
            z=scores[scores.config==config['name']]
            assert (z.ligand_props>=config['expr']-1e-7).all() and (z.receptor_props>=config['expr']-1e-7).all()
        state.update(status='completed',rows=len(scores),validation='unique finite bounded scores; expression eligibility; native compression parity previously checked')
    except Exception as exc:state.update(status='failed',error=type(exc).__name__+': '+str(exc));raise
    finally:state.update(updated_utc=datetime.now(timezone.utc).isoformat(),elapsed_seconds=round(time.monotonic()-start,1));write_json_atomic(out/'run_record.json',state)


if __name__=='__main__':main()
