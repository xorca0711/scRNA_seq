"""Stream deposited dense gene-row counts into a bounded annotation input.

Full-assay QC totals are retained. Repeated libraries stay distinct; biological
comparisons must aggregate by patient and histology after annotation.
"""
from pathlib import Path
import os,sys,gzip,json,time
from datetime import datetime,timezone
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ.setdefault(k,'2')
PAPER=Path(__file__).resolve().parents[1];ROOT=PAPER.parents[1];sys.path.insert(0,str(ROOT))
from analysis.lib.provenance import write_json_atomic,code_identity,sha256_file

def main():
    import numpy as np,pandas as pd,anndata as ad
    from scipy import sparse
    out=PAPER/'trials/u5_human_pilot';out.mkdir(exist_ok=True)
    cache=PAPER/'cache/u5_human_pilot';cache.mkdir(parents=True,exist_ok=True)
    if (cache/'annotation_input.h5ad').exists():raise RuntimeError('Pilot input exists; inspect before replacing')
    source=json.loads((PAPER/'trials/u2_extension_pilot/specification.json').read_text())
    jobs=[j for j in source['jobs'] if j['series']=='GSE308103']
    ref=ROOT/'Research Article/gate1_04_sikkema_2023_hlca/trials/s2_reference_mapping/reference'
    order=pd.read_csv(ref/'HLCA_reference_model_gene_order_ids_and_symbols.csv')
    assert order.gene_symbol.is_unique and order.gene_id.is_unique
    markers='EPCAM KRT8 KRT18 KRT19 SFTPC SFTPA1 SFTPB ABCA3 SLC34A2 NAPSA AGER HOPX PDPN CAV1 SCGB1A1 SCGB3A1 SCGB3A2 FOXJ1 PIFO KRT5 KRT14 KRT17 TP63 CLDN4 SFN KRT7 VIM PTPRC LST1 TYROBP LYZ C1QA C1QB CSF1R CD68 MRC1 FABP4 PPARG SPP1 S100A8 S100A9 FCGR3B CSF3R CD3D CD3E CD79A MS4A1 NKG7 COL1A1 COL1A2 COL3A1 DCN LUM PDGFRA COL14A1 PI16 ACTA2 TAGLN MYH11 RGS5 CSPG4 MCAM PECAM1 VWF EMCN KDR PROX1 IL1A IL1B IL1R1 IL1RAP IL1R2 IL1RN SIGIRR NLRP3 PYCARD CASP1 GSDMD MKI67 TOP2A'.split()
    pathways=pd.read_csv(PAPER/'trials/u5_ipf_spec/pathway_genes.tsv',sep='\t').gene.tolist()
    resources=pd.read_csv(PAPER/'trials/u5_full_source_panels/GSE136831/resource_edges.csv')
    rg=[g for col in ['ligand','receptor'] for x in resources[col] for g in x.split('_')]
    panel=sorted(set(order.gene_symbol)|set(markers)|set(pathways)|set(rg));gm={g:i for i,g in enumerate(panel)}
    spec=dict(frozen_utc=datetime.now(timezone.utc).isoformat(),patients=source['patients'],jobs=jobs,panel=panel,qc=dict(min_genes=500,min_counts=1000,max_mito_percent_exclusive=20),doublets='No author doublet labels deposited; broad annotation must retain ambiguous cells as unassigned; no malignancy/KAC labels inferred from atlas',model_gene_map_sha256=sha256_file(ref/'HLCA_reference_model_gene_order_ids_and_symbols.csv'),selection='First three numeric complete normal/AAH patients; all repeated libraries; no selection on biological effects')
    write_json_atomic(out/'input_specification.json',spec)
    state=dict(status='running',pid=os.getpid(),started_utc=spec['frozen_utc'],code=code_identity(ROOT,__file__),completed=[]);write_json_atomic(out/'input_run_record.json',state);start=time.monotonic()
    objects=[];qc=[];assays=[]
    try:
        for job in jobs:
            src=PAPER/'cache/extension_pilot'/job['url'].rsplit('/',1)[1]
            rr=[];cc=[];vv=[];seen=set()
            with gzip.open(src,'rt') as h:
                barcodes=h.readline().rstrip('\r\n').split('\t');n=len(barcodes)
                assert len(set(barcodes))==n
                total=np.zeros(n,dtype=np.int64);ng=np.zeros(n,dtype=np.int32);mt=np.zeros(n,dtype=np.int64);maxima=np.zeros(n,dtype=np.int64)
                for line in h:
                    gene,values=line.rstrip('\r\n').split('\t',1)
                    if gene in seen:raise ValueError('Duplicate source gene symbol '+gene)
                    seen.add(gene);v=np.fromstring(values,sep='\t',dtype=np.int64)
                    assert len(v)==n and np.all(v>=0)
                    total+=v;ng+=v>0;np.maximum(maxima,v,out=maxima)
                    if gene.startswith('MT-'):mt+=v
                    if gene in gm:
                        keep=np.flatnonzero(v);rr.append(keep);cc.append(np.full(len(keep),gm[gene],dtype=np.int32));vv.append(v[keep].astype(np.float32))
            x=sparse.coo_matrix((np.concatenate(vv),(np.concatenate(rr),np.concatenate(cc))),shape=(n,len(panel))).tocsr()
            obs=pd.DataFrame(dict(barcode=barcodes,gsm=job['gsm'],patient=str(job['patient']),histology=job['histology'],full_library_size=total,n_genes_by_counts=ng,full_gene_max_count=maxima,pct_counts_mt=100*mt/np.maximum(total,1)))
            obs.index=job['gsm']+':'+obs.barcode;obs.index.name='cell_id'
            obs['qc_pass']=(obs.n_genes_by_counts>=500)&(obs.full_library_size>=1000)&(obs.pct_counts_mt<20)
            a=ad.AnnData(x,obs=obs,var=pd.DataFrame(index=panel));a.write_h5ad(cache/(job['gsm']+'_panel.h5ad'),compression='gzip')
            objects.append(a[a.obs.qc_pass].copy())
            qc.append(dict(gsm=job['gsm'],patient=job['patient'],histology=job['histology'],source_cells=n,qc_cells=int(obs.qc_pass.sum()),assayed_genes=len(seen),model_genes_measured=int(order.gene_symbol.isin(seen).sum()),total_counts=int(total.sum()),sha256=sha256_file(src)))
            for g in panel:assays.append(dict(gsm=job['gsm'],gene=g,assayed=g in seen))
            state['completed'].append(qc[-1]);write_json_atomic(out/'input_run_record.json',state);print(job['gsm'],n,'input',int(obs.qc_pass.sum()),'QC',flush=True)
        full=ad.concat(objects,join='inner',merge='same');full.write_h5ad(cache/'annotation_input.h5ad',compression='gzip')
        pd.DataFrame(qc).to_csv(out/'library_qc.csv',index=False);pd.DataFrame(assays).to_csv(out/'panel_assay_coverage.csv',index=False)
        assert full.obs_names.is_unique and len(full)==sum(x['qc_cells'] for x in qc)
        state.update(status='completed',qc_cells=full.n_obs,panel_genes=full.n_vars)
    except Exception as exc:state.update(status='failed',error=str(exc));raise
    finally:state.update(elapsed_seconds=round(time.monotonic()-start,1),updated_utc=datetime.now(timezone.utc).isoformat());write_json_atomic(out/'input_run_record.json',state)

if __name__=='__main__':main()
