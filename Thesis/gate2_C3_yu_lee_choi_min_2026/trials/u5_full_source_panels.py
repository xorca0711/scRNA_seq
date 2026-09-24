"""Stream full raw matrices into exact LR component sufficient inputs.

Retain every source cell's full library size and actual largest gene count.
The raw matrices must have unique symbols and strictly sorted unique (cell,gene)
coordinates; otherwise stop rather than approximate duplicate handling.
"""
from pathlib import Path
import gzip,json,os,sys,time
from datetime import datetime,timezone
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ.setdefault(key,'2')
PAPER=Path(__file__).resolve().parents[1];ROOT=PAPER.parents[1]
sys.path.insert(0,str(ROOT))
from analysis.lib.provenance import sha256_file,write_json_atomic,code_identity


def main():
    import numpy as np,pandas as pd,anndata as ad
    from scipy import sparse
    cohort=sys.argv[1]
    out=PAPER/'trials/u5_full_source_panels'/cohort;out.mkdir(parents=True,exist_ok=True)
    cache=PAPER/'cache/u5_full_source_panels';cache.mkdir(parents=True,exist_ok=True)
    target=cache/(cohort+'.h5ad')
    if target.exists():raise RuntimeError('Extract already exists; do not overwrite')
    sys.path.insert(0,str(ROOT/'Thesis/gate1_01_niethamer_2025/trials'))
    import g2_gsea_ipf as loaders
    class Inputs:
        def __init__(self):self.paths=[]
        def add_input(self,p):self.paths.append(Path(p))
    inputs=Inputs()
    meta,genes,raw,_=(loaders.cohort_136831 if cohort=='GSE136831' else loaders.cohort_135893)(inputs)
    meta=meta.reset_index(drop=True).rename(columns={'CellBarcode_Identity':'barcode'})
    if 'barcode' not in meta:meta['barcode']=[cohort+'_'+str(i) for i in range(len(meta))]
    meta['disease']=meta.disease.replace({'Control':'control'})
    assert len(set(genes))==len(genes),'Symbol collapsing needed before exact maximum extraction'
    omni=ROOT/'.venv-x64/Lib/site-packages/liana/resource/omni_resource.csv'
    res=pd.read_csv(omni)
    ligands=['IL1B','IL1A','AREG','HBEGF','TGFB1','CCL2','CXCL12','CXCL1','CXCL2','FGF7','FGF10','VEGFA']
    res=res[res.resource.isin(['consensus','cellchatdb'])].rename(columns={'source_genesymbol':'ligand','target_genesymbol':'receptor'})
    res=res[res.ligand.map(lambda x:any(g in ligands or g.startswith('WNT') for g in x.split('_'))) ][['resource','ligand','receptor']].drop_duplicates()
    panel=sorted({g for col in ['ligand','receptor'] for entry in res[col] for g in entry.split('_')}|{'IL1RN','IL1R2','IL1RAP','SIGIRR','NLRP3','PYCARD','CASP1','GSDMD','TNF','IL6','SPP1','CSF1','CSF1R','EGFR','TGFBR1','TGFBR2'})
    measured=sorted(set(panel)&set(genes));gm={g:i for i,g in enumerate(measured)}
    gene_map=np.array([gm.get(g,-1) for g in genes],dtype=np.int32)
    spec={'frozen_utc':datetime.now(timezone.utc).isoformat(),'cohort':cohort,'cells':len(meta),'genes':len(genes),'panel':panel,
          'missing_panel_genes':sorted(set(panel)-set(genes)),'scope':'all deposited cells, donor labels retained; disease comparisons only IPF/control',
          'exactness':'strictly unique symbols and (cell,gene) coordinates; full-library totals and observed per-cell maximum retained',
          'input_hashes':{str(p.relative_to(ROOT)):sha256_file(p) for p in inputs.paths+[omni]}}
    write_json_atomic(out/'specification.json',spec);res.to_csv(out/'resource_edges.csv',index=False)
    state={'status':'running','stage':'streaming','pid':os.getpid(),'started_utc':spec['frozen_utc'],'code':code_identity(ROOT,__file__)}
    write_json_atomic(out/'run_record.json',state);start=time.monotonic()
    try:
        with gzip.open(raw,'rt') as h:
            header=h.readline();skip=1;line=h.readline();skip+=1
            assert 'coordinate integer' in header
            while line.startswith('%'):line=h.readline();skip+=1
            ng,nc,nnz=map(int,line.split())
        assert (ng,nc)==(len(genes),len(meta))
        totals=np.zeros(nc,dtype=np.int64);maxima=np.zeros(nc,dtype=np.int64)
        rr=[];cc=[];vv=[];seen=0;last_key=-1
        for chunk in pd.read_csv(raw,sep=r'\s+',skiprows=skip,header=None,names=['g','c','v'],dtype={'g':np.int32,'c':np.int32,'v':np.int64},chunksize=2_000_000):
            g=chunk.g.to_numpy()-1;c=chunk.c.to_numpy()-1;v=chunk.v.to_numpy()
            if np.any(v<0) or np.any(g<0) or np.any(g>=ng) or np.any(c<0) or np.any(c>=nc):raise ValueError('Invalid count coordinate')
            key=c.astype(np.int64)*ng+g
            if key[0]<=last_key or np.any(key[1:]<=key[:-1]):raise ValueError('Source is not strictly unique cell/gene sorted; exact compression refused')
            last_key=int(key[-1]);np.add.at(totals,c,v);np.maximum.at(maxima,c,v)
            mapped=gene_map[g];keep=mapped>=0
            rr.append(c[keep]);cc.append(mapped[keep]);vv.append(v[keep].astype(np.float32))
            seen+=len(chunk)
            if seen%100_000_000==0:
                state.update(entries_read=seen,entries_total=nnz,elapsed_seconds=round(time.monotonic()-start,1));write_json_atomic(out/'run_record.json',state)
                print(cohort,seen,nnz,flush=True)
        assert seen==nnz
        x=sparse.coo_matrix((np.concatenate(vv),(np.concatenate(rr),np.concatenate(cc))),shape=(nc,len(measured))).tocsr()
        del rr,cc,vv
        meta['full_library_size']=totals;meta['full_gene_max_count']=maxima;meta.index=meta.barcode.astype(str);meta.index.name='cell_id'
        # Compare complete subtype library totals to already audited raw-count pseudobulks.
        units=pd.read_csv(ROOT/'analysis/corrections/statistics/cache'/(cohort+'_subtypes_units.csv'))
        with np.load(ROOT/'analysis/corrections/statistics/cache'/(cohort+'_subtypes.npz')) as prior:expected=prior['counts'].sum(axis=1)
        observed=meta.groupby(['donor','label']).full_library_size.sum()
        for i,r in units.iterrows():assert observed.loc[(r.donor,r.label)]==expected[i]
        a=ad.AnnData(x,obs=meta,var=pd.DataFrame(index=measured));a.write_h5ad(target,compression='gzip')
        state.update(status='completed',stage='exact_component_input_ready',entries_read=seen,total_counts=int(totals.sum()),cells=nc,panel_genes=len(measured),count_parity='passed',output_sha256=sha256_file(target))
    except Exception as exc:state.update(status='failed',error=str(exc));raise
    finally:
        state.update(updated_utc=datetime.now(timezone.utc).isoformat(),elapsed_seconds=round(time.monotonic()-start,1));write_json_atomic(out/'run_record.json',state)
    print(cohort,'completed exact full-source panel',nc,len(measured),flush=True)


if __name__=='__main__':main()
