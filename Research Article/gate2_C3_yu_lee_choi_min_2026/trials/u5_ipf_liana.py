"""Bounded GSE136831 triad extraction and per-donor LIANA descriptives.

Uses deposited labels. No cell-permutation or treatment-causal inference.
The 500-cell cap is a reproducible initial descriptive pass, not all-cell validation.
"""
from __future__ import annotations
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:
    os.environ.setdefault(key,'2')
import gzip
import json
from pathlib import Path
import sys
import time
from datetime import datetime,timezone
PAPER=Path(__file__).resolve().parents[1]
ROOT=PAPER.parents[1]
sys.path.insert(0,str(ROOT))
from analysis.lib.provenance import code_identity,sha256_file,write_json_atomic,archive_existing_record


def main():
    import numpy as np
    import pandas as pd
    from scipy import sparse
    import anndata as ad
    out=PAPER/'trials/u5_ipf_liana';out.mkdir(exist_ok=True)
    cache=PAPER/'cache/u5_ipf_liana';cache.mkdir(exist_ok=True)
    if not json.loads((PAPER/'analysis_contract.json').read_text())['full_run_authorized']:
        raise RuntimeError('Run not authorized')
    start=time.monotonic();record=out/'run_record.json'
    archive_existing_record(record)
    state={'status':'running','stage':'freeze_inputs','pid':os.getpid(),'started_utc':datetime.now(timezone.utc).isoformat(),
           'code':code_identity(ROOT,__file__),'completed_donors':[],'evidence':'per-donor RNA coexpression descriptives, previously inspected observational cohort'}
    write_json_atomic(record,state)
    try:
        raw=ROOT/'raw_data/GSE136831'
        paths={'metadata':raw/'GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz',
               'genes':raw/'GSE136831_AllCells.GeneIDs.txt.gz', 'matrix':raw/'GSE136831_RawCounts_Sparse.mtx.gz'}
        meta=pd.read_csv(paths['metadata'],sep='\t',usecols=['CellBarcode_Identity','Manuscript_Identity','Disease_Identity','Subject_Identity'])
        meta=meta.rename(columns={'CellBarcode_Identity':'barcode','Manuscript_Identity':'celltype','Disease_Identity':'disease','Subject_Identity':'donor'})
        meta['disease']=meta.disease.replace({'Control':'control'})
        assert meta.barcode.is_unique
        labels={'ATII':'epithelial','Fibroblast':'fibroblast','Myofibroblast':'fibroblast','Macrophage':'myeloid','Macrophage_Alveolar':'myeloid'}
        meta['compartment']=meta.celltype.map(labels)
        eligible=meta[meta.compartment.notna() & meta.disease.isin(['IPF','control'])]
        coverage=eligible.groupby(['donor','disease','celltype','compartment']).size().rename('cells').reset_index()
        coverage['eligible']=coverage.cells>=50
        coverage.to_csv(out/'donor_subtype_coverage.csv',index=False)
        good=coverage[coverage.eligible]
        directions={('myeloid','fibroblast'),('myeloid','epithelial'),('epithelial','fibroblast'),('fibroblast','myeloid'),('fibroblast','epithelial')}
        pairs=[]
        for donor,rows in good.groupby('donor',sort=True):
            for a in rows.itertuples():
                for b in rows.itertuples():
                    if (a.compartment,b.compartment) in directions:
                        pairs.append({'donor':donor,'source':a.celltype,'target':b.celltype,'disease':a.disease})
        pairs=pd.DataFrame(pairs);pairs.to_csv(out/'fixed_donor_pairs.csv',index=False)
        donors=sorted(pairs.donor.unique())
        groups=good[good.donor.isin(donors)].copy().reset_index(drop=True)
        group_map={(r.donor,r.celltype):i for i,r in groups.iterrows()}
        goc=np.array([group_map.get((r.donor,r.celltype),-1) for r in meta.itertuples()],dtype=np.int32)
        rng=np.random.default_rng(20260924);selected=[]
        for idx in range(len(groups)):
            rows=np.flatnonzero(goc==idx)
            selected.extend(rng.choice(rows,min(500,len(rows)),replace=False).tolist())
        selected=np.array(sorted(selected),dtype=np.int64)
        cell_map=np.full(len(meta),-1,dtype=np.int32);cell_map[selected]=np.arange(len(selected))
        original_genes=pd.read_csv(paths['genes'],sep='\t').iloc[:,1].astype(str).to_numpy()
        unique_genes,inverse=np.unique(original_genes,return_inverse=True)
        resource_file=ROOT/'.venv-x64/Lib/site-packages/liana/resource/omni_resource.csv'
        resources=pd.read_csv(resource_file)
        core=['IL1B','IL1A','AREG','HBEGF','TGFB1']
        extended=['CCL2','CXCL12','CXCL1','CXCL2','FGF7','FGF10','VEGFA']
        resource_paths=[]
        for name in ['consensus','cellchatdb']:
            table=resources[resources.resource.eq(name)][['source_genesymbol','target_genesymbol']].rename(columns={'source_genesymbol':'ligand','target_genesymbol':'receptor'}).drop_duplicates()
            table=table[table.ligand.map(lambda x: any(g in core+extended for g in x.split('_')))]
            table['family']=table.ligand.map(lambda x:'core' if any(g in core for g in x.split('_')) else 'exploratory')
            path=out/('resource_'+name+'.csv');table.to_csv(path,index=False);resource_paths.append(path)
        spec={'frozen_utc':datetime.now(timezone.utc).isoformat(),'cohort':'GSE136831','sampled_cells':len(selected),
              'donors':donors,'cell_floor':50,'maximum_cells_per_donor_subtype':500,'seed':20260924,
              'method':'LIANA_CellChat_like_magnitude_n_perms_None','resources':['consensus','cellchatdb'],'expression_fraction':0.1,
              'scope':'declared myeloid/fibroblast/AT2 pairs; no whole-cohort source-agnostic claim; initial descriptive pass only',
              'input_hashes':{str(p.relative_to(ROOT)):sha256_file(p) for p in list(paths.values())+[resource_file]+resource_paths}}
        write_json_atomic(out/'specification.json',spec)
        groups.to_csv(out/'selected_groups.csv',index=False)
        state.update(stage='streaming_raw_counts',planned_donors=len(donors),sampled_cells=len(selected),spec_sha256=sha256_file(out/'specification.json'))
        write_json_atomic(record,state)
        target=cache/'sampled_triad.h5ad'
        # Build one sparse subset in a single pass through the unchanged source.
        with gzip.open(paths['matrix'],'rt') as handle:
            header=handle.readline();skip=1
            if not header.startswith('%%MatrixMarket matrix coordinate'): raise ValueError('Unexpected matrix format')
            line=handle.readline();skip+=1
            while line.startswith('%'): line=handle.readline();skip+=1
            ng,nc,nnz=map(int,line.split())
        if (ng,nc)!=(len(original_genes),len(meta)): raise ValueError('Metadata/matrix dimensions differ')
        rows_list=[];cols_list=[];values_list=[];seen=0
        group_totals=np.zeros(len(groups),dtype=np.int64)
        for chunk in pd.read_csv(paths['matrix'],sep=r'\s+',skiprows=skip,header=None,names=['g','c','v'],dtype={'g':np.int32,'c':np.int32,'v':np.int64},chunksize=2_000_000):
            seen+=len(chunk);cc=chunk.c.to_numpy()-1;gg=chunk.g.to_numpy()-1;vv=chunk.v.to_numpy()
            if np.any(cc<0) or np.any(cc>=nc) or np.any(gg<0) or np.any(gg>=ng) or np.any(vv<0): raise ValueError('Invalid matrix entry')
            gx=goc[cc];all_keep=gx>=0
            np.add.at(group_totals,gx[all_keep],vv[all_keep])
            mapped=cell_map[cc];keep=mapped>=0
            rows_list.append(mapped[keep]);cols_list.append(inverse[gg[keep]].astype(np.int32));values_list.append(vv[keep].astype(np.float32))
            if seen % 100_000_000==0:
                state.update(entries_read=seen,entries_total=nnz,elapsed_seconds=round(time.monotonic()-start,2));write_json_atomic(record,state)
                print(f'LR extraction {seen}/{nnz} entries',flush=True)
        if seen!=nnz: raise ValueError('Matrix entry count mismatch')
        prior_cache=ROOT/'analysis/corrections/statistics/cache'
        units=pd.read_csv(prior_cache/'GSE136831_subtypes_units.csv')
        with np.load(prior_cache/'GSE136831_subtypes.npz') as old:
            totals=old['counts'].sum(axis=1)
        expected={(r.donor,r.label):int(totals[i]) for i,r in units.iterrows()}
        for i,r in groups.iterrows():
            if group_totals[i]!=expected[(r.donor,r.celltype)]: raise ValueError('Full-source totals do not match audited pseudobulks')
        x=sparse.coo_matrix((np.concatenate(values_list),(np.concatenate(rows_list),np.concatenate(cols_list))),shape=(len(selected),len(unique_genes))).tocsr()
        obs=meta.iloc[selected].copy();obs.index=obs.barcode.astype(str);obs.index.name='cell_id'
        data=ad.AnnData(x,obs=obs,var=pd.DataFrame(index=unique_genes));data.write_h5ad(target,compression='gzip')
        del rows_list,cols_list,values_list,x
        import scanpy as sc
        from liana.method import cellchat
        state.update(stage='per_donor_ligand_receptor',entries_read=seen,input_count_parity='passed')
        write_json_atomic(record,state)
        for donor in donors:
            part=data[data.obs.donor.eq(donor)].copy()
            sc.pp.normalize_total(part,target_sum=10000);sc.pp.log1p(part)
            pairs_d=pairs[pairs.donor.eq(donor)][['source','target']].drop_duplicates()
            for resource in ['consensus','cellchatdb']:
                table=pd.read_csv(out/('resource_'+resource+'.csv'))
                scored=cellchat(part,groupby='celltype',resource=table[['ligand','receptor']],groupby_pairs=pairs_d,
                               expr_prop=0.1,min_cells=50,use_raw=False,n_perms=None,verbose=False,inplace=False)
                scored['donor']=donor;scored['disease']=part.obs.disease.iloc[0];scored['resource']=resource
                scored.to_csv(out/(resource+'_'+donor+'.csv'),index=False)
            state['completed_donors'].append(donor);state['elapsed_seconds']=round(time.monotonic()-start,2)
            write_json_atomic(record,state)
            print(f'LR donors {len(state["completed_donors"])}/{len(donors)}',flush=True)
        state['status']='completed_descriptive_LR_requires_interpretation'
    except Exception as exc:
        state.update(status='failed',error=type(exc).__name__+': '+str(exc));raise
    finally:
        state.update(updated_utc=datetime.now(timezone.utc).isoformat(),elapsed_seconds=round(time.monotonic()-start,2));write_json_atomic(record,state)


if __name__=='__main__': main()
