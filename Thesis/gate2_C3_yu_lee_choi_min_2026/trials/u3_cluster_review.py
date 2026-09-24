"""Treatment-blind cluster diagnostics for annotation review; no effect tests."""
from __future__ import annotations
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:
    os.environ.setdefault(key,'2')
import json
import sys
import time
from datetime import datetime,timezone
from pathlib import Path
PAPER=Path(__file__).resolve().parents[1]
ROOT=PAPER.parents[1]
sys.path.insert(0,str(ROOT))
from analysis.lib.provenance import code_identity, sha256_file, write_json_atomic, archive_existing_record


def main():
    import anndata as ad
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from scipy import sparse
    start=time.monotonic()
    spec=json.loads((PAPER/'trials/u3_run_spec.json').read_text(encoding='utf-8'))
    out=PAPER/'trials/u3_cluster_review'
    out.mkdir(parents=True,exist_ok=True)
    status=out/'run_record.json'
    archive_existing_record(status)
    params={'n_top_genes':1500,'n_pcs':30,'neighbors':15,'leiden_resolution':1.0,'seed':20260924,
            'feature_selection':'seurat_log_normalized_dispersion_batch_key_GSM',
            'integration':False,'treatment_used_for_clustering':False,'biological_effect_tests':False}
    state={'started_utc':datetime.now(timezone.utc).isoformat(),'status':'running','pid':os.getpid(),
           'stage':'U3_cluster_annotation_review','parameters':params,'code':code_identity(ROOT,__file__)}
    write_json_atomic(status,state)
    try:
        objects=[]; inputs=[]
        for s in spec['samples']:
            path=PAPER/'cache/u3_qc'/(s['gsm']+'_qc.h5ad')
            d=ad.read_h5ad(path)
            objects.append(d[d.obs.qc_pass].copy())
            inputs.append({'file':str(path.relative_to(ROOT)),'sha256':sha256_file(path)})
        d=ad.concat(objects,join='inner',merge='same')
        del objects
        # Retain full counts in original per-sample caches; this object is a diagnostic transform.
        d.X=d.X.astype(np.float32)
        sc.pp.normalize_total(d,target_sum=10000)
        sc.pp.log1p(d)
        sc.pp.highly_variable_genes(d,n_top_genes=params['n_top_genes'],flavor='seurat',batch_key='gsm')
        (out/'hvg_genes.txt').write_text('\n'.join(d.var_names[d.var.highly_variable])+'\n',encoding='utf-8')
        sc.pp.pca(d,n_comps=params['n_pcs'],random_state=params['seed'],svd_solver='arpack')
        sc.pp.neighbors(d,n_neighbors=params['neighbors'],n_pcs=params['n_pcs'],random_state=params['seed'])
        sc.tl.leiden(d,resolution=params['leiden_resolution'],random_state=params['seed'],flavor='igraph',n_iterations=2,directed=False,key_added='review_cluster')
        clusters=d.obs.review_cluster.astype(str)
        levels=sorted(clusters.unique(),key=int)
        codes=pd.Categorical(clusters,categories=levels).codes
        group=sparse.csr_matrix((np.ones(d.n_obs),(codes,np.arange(d.n_obs))),shape=(len(levels),d.n_obs))
        n=np.asarray(group.sum(axis=1)).ravel()
        sums=(group@d.X).toarray()
        means=sums/n[:,None]
        other=(np.asarray(d.X.sum(axis=0))-sums)/(d.n_obs-n[:,None])
        difference=means-other
        symbols=d.var.gene_symbol.to_numpy()
        rows=[]
        for i,cluster in enumerate(levels):
            top=np.argsort(-difference[i])[:35]
            for rank,j in enumerate(top,1):
                rows.append({'cluster':cluster,'rank':rank,'feature_id':d.var_names[j],'gene_symbol':symbols[j],
                             'mean_log_expression':float(means[i,j]),'mean_log_difference_vs_other_cells':float(difference[i,j]),
                             'descriptive_no_P_value':True})
        pd.DataFrame(rows).to_csv(out/'cluster_top_genes.csv',index=False)
        marker_genes=list(dict.fromkeys(g for panel in spec['candidate_panels'].values() for g in panel))
        # Additional lineage checks are reported, not used to force labels.
        marker_genes+=['Sftpc','Sftpa1','Sftpb','Abca3','Ager','Pdpn','Hopx','Scgb1a1','Foxj1','Krt5','Trp63','Cd74','H2-Ab1','C1qa','C1qb','Pparg','Fabp4','Col13a1','Col14a1','Pi16','Acta2','Tagln','Rgs5','Mki67','Nkg7']
        markers=[]
        for gene in dict.fromkeys(marker_genes):
            positions=np.flatnonzero(symbols==gene)
            if len(positions):
                detected=(group@(d.X[:,positions].sum(axis=1)>0).astype(float)).A1 if hasattr(group@(d.X[:,positions].sum(axis=1)>0).astype(float),'A1') else np.asarray(group@(d.X[:,positions].sum(axis=1)>0).astype(float)).ravel()
                for i,cluster in enumerate(levels):
                    markers.append({'cluster':cluster,'gene_symbol':gene,'mean_log_expression':float(means[i,positions].sum()),'fraction_detected':float(detected[i]/n[i])})
        pd.DataFrame(markers).to_csv(out/'cluster_marker_summary.csv',index=False)
        pd.crosstab(clusters,d.obs.gsm).to_csv(out/'cluster_sample_counts.csv',index_label='cluster')
        pd.crosstab(clusters,d.obs.candidate_compartment).to_csv(out/'cluster_candidate_crosswalk.csv',index_label='cluster')
        d.obs[['gsm','review_cluster','candidate_compartment']].to_csv(PAPER/'cache/u3_qc/cluster_membership.csv.gz',compression='gzip',index_label='cell_id')
        d.write_h5ad(PAPER/'cache/u3_qc/cluster_diagnostics.h5ad',compression='gzip')
        state.update(status='completed_annotation_review_required',cells=d.n_obs,clusters=len(levels),inputs=inputs)
        print(json.dumps({'status':state['status'],'cells':d.n_obs,'clusters':len(levels),'elapsed_seconds':round(time.monotonic()-start,2)}),flush=True)
    except Exception as exc:
        state.update(status='failed',error=type(exc).__name__+': '+str(exc))
        raise
    finally:
        state.update(updated_utc=datetime.now(timezone.utc).isoformat(),elapsed_seconds=round(time.monotonic()-start,2))
        write_json_atomic(status,state)


if __name__=='__main__':
    main()
