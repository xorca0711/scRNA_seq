"""Native sample-specific mouse LR magnitude with reviewed broad lineages."""
from pathlib import Path
import os,sys,json,time,warnings
from datetime import datetime,timezone
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ.setdefault(k,'2')
PAPER=Path(__file__).resolve().parents[1];ROOT=PAPER.parents[1];sys.path.insert(0,str(ROOT))
from analysis.lib.provenance import write_json_atomic,code_identity,sha256_file
from u5_liana_exact_components import normalized_components
from u5_liana_robustness import CONFIGS,DIRECTIONS

def main():
    import numpy as np,pandas as pd,anndata as ad
    from scipy import sparse
    from liana.method import cellchat
    warnings.filterwarnings('ignore',category=FutureWarning,module='liana');warnings.filterwarnings('ignore',category=ad.ImplicitModificationWarning)
    out=PAPER/'trials/u4_mouse_niche';out.mkdir(exist_ok=True)
    cache=PAPER/'cache/u4_mouse_niche';cache.mkdir(parents=True,exist_ok=True)
    record=out/'liana_run_record.json';assert not record.exists()
    res=pd.read_csv(PAPER/'trials/u4_resources/mouse_resource_edges.csv');panel=sorted({g for col in ['ligand','receptor'] for value in res[col] for g in value.split('_')})
    annotations=pd.read_csv(PAPER/'cache/u3_lineage_annotation/cell_annotations.csv.gz',index_col=0)
    spec=dict(frozen_utc=datetime.now(timezone.utc).isoformat(),configs=CONFIGS,species='Mus musculus',source_annotation='independent broad-lineage review; myeloid is not macrophage-only; alveolar is not AT2/KAC',method='native LIANA CellChat-like, n_perms=None',directions=sorted(DIRECTIONS),resource_sha256=sha256_file(PAPER/'trials/u4_resources/mouse_resource_edges.csv'),input_annotation_sha256=sha256_file(PAPER/'cache/u3_lineage_annotation/cell_annotations.csv.gz'),inference='within-sample descriptive magnitude only, no cell-level P values')
    write_json_atomic(out/'liana_specification.json',spec);state=dict(status='running',pid=os.getpid(),started_utc=spec['frozen_utc'],code=code_identity(ROOT,__file__),completed_calls=0);write_json_atomic(record,state);start=time.monotonic()
    all_scores=[];calls=[];pairs_out=[]
    try:
        for gsm in sorted(annotations.gsm.unique()):
            a=ad.read_h5ad(PAPER/'cache/u3_qc'/(gsm+'_qc.h5ad'));a=a[a.obs.qc_pass].copy();m=annotations.loc[a.obs_names].copy()
            symbols=a.var.gene_symbol.astype(str).to_numpy();genes=sorted(set(symbols));gm={g:i for i,g in enumerate(genes)}
            project=sparse.csr_matrix((np.ones(len(symbols),dtype=np.int32),(np.arange(len(symbols)),[gm[g] for g in symbols])),shape=(len(symbols),len(genes)))
            x=(a.X@project).tocsr();total=np.asarray(x.sum(axis=1)).ravel();maxima=x.max(axis=1).toarray().ravel()
            keep=[g for g in panel if g in gm]
            part=ad.AnnData(x[:,[gm[g] for g in keep]],obs=m,var=pd.DataFrame(index=keep));del x,a
            part.obs['full_library_size']=total;part.obs['full_gene_max_count']=maxima;part.obs['celltype']=part.obs.label
            part.obs['compartment']=part.obs.compartment.replace({'alveolar':'epithelial'})
            part=part[part.obs.compartment.isin(['epithelial','fibroblast','myeloid'])].copy();part.write_h5ad(cache/(gsm+'_lr_components.h5ad'),compression='gzip')
            for config in CONFIGS:
                coverage=part.obs.groupby(['celltype','compartment'],observed=True).size().rename('cells').reset_index();good=coverage[coverage.cells>=config['floor']]
                pairs=pd.DataFrame([dict(source=r.celltype,target=s.celltype) for r in good.itertuples() for s in good.itertuples() if (r.compartment,s.compartment) in DIRECTIONS])
                if pairs.empty:continue
                pairs=pairs.drop_duplicates();selected=[];rng=np.random.default_rng(config['seed'])
                for label in good.celltype:
                    indices=np.flatnonzero(part.obs.celltype.eq(label))
                    if config['cap'] and len(indices)>config['cap']:indices=rng.choice(indices,config['cap'],replace=False)
                    selected.extend(indices)
                z=normalized_components(part[np.sort(selected)])
                props={label:pd.Series(np.asarray((z.X[np.flatnonzero(z.obs.celltype.eq(label))]>0).mean(axis=0)).ravel(),index=z.var_names) for label in good.celltype}
                for p in pairs.itertuples(index=False):pairs_out.append(dict(gsm=gsm,config=config['name'],source=p.source,target=p.target))
                for resource in ['consensus','cellchatdb']:
                    table=res[res.resource==resource];eligible_edges=[]
                    for p in pairs.itertuples(index=False):
                        for edge in table.itertuples(index=False):
                            vals=[props[p.source].get(g,-1) for g in edge.ligand.split('_')]+[props[p.target].get(g,-1) for g in edge.receptor.split('_')]
                            if min(vals)>=config['expr']:eligible_edges.append((p.source,p.target,edge.ligand,edge.receptor))
                    if eligible_edges:
                        scored=cellchat(z,groupby='celltype',resource=table[['ligand','receptor']],groupby_pairs=pairs,expr_prop=config['expr'],min_cells=config['floor'],use_raw=False,n_perms=None,verbose=False,inplace=False)
                        assert set(scored[['source','target','ligand_complex','receptor_complex']].itertuples(index=False,name=None))==set(eligible_edges)
                        scored=scored.merge(table[['ligand','receptor','human_ligand','human_receptor','family']],left_on=['ligand_complex','receptor_complex'],right_on=['ligand','receptor'],validate='many_to_one')
                        for k,v in dict(gsm=gsm,treatment=m.treatment.iloc[0],resource=resource,config=config['name']).items():scored[k]=v
                        all_scores.append(scored)
                    calls.append(dict(gsm=gsm,treatment=m.treatment.iloc[0],resource=resource,config=config['name'],cells=z.n_obs,pairs=len(pairs),retained_edges=len(eligible_edges)))
                    state['completed_calls']+=1
            state.update(gsm=gsm,elapsed_seconds=round(time.monotonic()-start,1));write_json_atomic(record,state);print(gsm,'LR completed',flush=True)
        scores=pd.concat(all_scores,ignore_index=True)
        assert not scores.duplicated(['gsm','resource','config','source','target','ligand_complex','receptor_complex']).any()
        assert np.isfinite(scores.lr_probs).all() and scores.lr_probs.between(0,1).all()
        assert not any('pval' in c.lower() for c in scores.columns)
        scores.to_csv(out/'donor_lr_scores.csv.gz',index=False,compression='gzip');pd.DataFrame(calls).to_csv(out/'call_coverage.csv',index=False);pd.DataFrame(pairs_out).to_csv(out/'fixed_pairs_by_config.csv',index=False)
        state.update(status='completed',rows=len(scores),validation='all native edges equal independent subunit/pair eligibility; unique bounded scores; no P values')
    except Exception as exc:state.update(status='failed',error=str(exc));raise
    finally:state.update(updated_utc=datetime.now(timezone.utc).isoformat(),elapsed_seconds=round(time.monotonic()-start,1));write_json_atomic(record,state)

if __name__=='__main__':main()
