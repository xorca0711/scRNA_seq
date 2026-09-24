"""Freeze treatment-blind broad-lineage review; produce raw-count pseudobulks.

These independent annotations do not reconstruct the authors' KAC classifier.
Mixed-lineage clusters stay unassigned. A broad CSF1R+ myeloid compartment is
not relabelled as a pure macrophage population.
"""
from pathlib import Path
import os,sys,json,time
from datetime import datetime,timezone
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ.setdefault(k,'2')
PAPER=Path(__file__).resolve().parents[1];ROOT=PAPER.parents[1];sys.path.insert(0,str(ROOT))
from analysis.lib.provenance import write_json_atomic,code_identity,sha256_file

ANNOTATIONS={
0:('myeloid','CSF1R_myeloid','Lyz2/Tyrobp/Csf1r; insufficient evidence for macrophage-only identity'),
1:('myeloid','CSF1R_myeloid','Csf1r/Lyz2/Cd68; C1q and Mertk low, monocyte/macrophage unresolved'),
2:('other','APC_ambiguous','MHCII-rich; DC/macrophage distinction unresolved'),
3:('unassigned','mixed_lineage','MHCII/airway/stromal and immunoglobulin mixture'),
4:('myeloid','inflammatory_CSF1R_CSF3R_myeloid','Strong Csf1r and Csf3r; preserve as distinct inflammatory myeloid subtype'),
5:('other','pDC_candidate','Bst2/Siglech/Tcf4'),
6:('myeloid','macrophage_candidate','C1qa/C1qb/Mertk with Lyz2; only 19 cells overall'),
7:('fibroblast','fibroblast_Gsn','Pdgfra/Col1a1/Dcn; Gsn/Bgn/Mgp'),
8:('other','T_cell','Cd3d/Cd3e/Il7r'),9:('other','T_cell','Cd3d/Cd3e/Thy1'),10:('other','T_cell','Cd3e/Il7r'),
11:('other','B_cell','Cd79a/Cd79b/Igkc'),12:('other','NK_cell','Ncr1/Prf1/Nkg7'),13:('other','cytotoxic_T_cell','Cd3e/Nkg7/Ccl5'),
14:('fibroblast','fibroblast_Itga8','Pdgfra/Col1a1/Dcn; Itga8/Npnt/C7'),
15:('other','neutrophil','Csf3r/Retnlg/S100a8/S100a9'),16:('other','inflammatory_granulocyte_candidate','Csf3r high, Csf1r low; Il1b/Cxcl2/Ccl3'),
17:('other','basophil_mast_candidate','Gata2/Hdc; lineage ambiguous'),
18:('unassigned','mixed_myeloid_B','Lyz2/Tyrobp plus Cd79a/MHCII'),
19:('endothelial','vascular_endothelial','Ptprb/Plvap/Cdh5'),
20:('alveolar','mixed_AT2_AT1_alveolar','Sftpc/Sftpa1/Abca3 together with Ager/Hopx; no KAC or DATP assignment'),
21:('airway','secretory_airway','Scgb1a1/Cbr2/Wfdc2'),
22:('other','smooth_muscle','Acta2/Tagln/Myh11; not assigned as myofibroblast'),
23:('alveolar','AT1_like','Ager/Hopx/Pdpn/Aqp5; cell identity only, no lineage-fate claim'),
24:('other','mesothelial_stromal_ambiguous','Dcn/C3/Upk3b; not forced into fibroblast'),
25:('endothelial','capillary_endothelial','Car4/Cldn5/Kdr/Cdh5'),
26:('fibroblast','fibroblast_Pi16','Pdgfra/Col1a1/Dcn plus Pi16/Clec3b'),
27:('other','mural_pericyte','Pdgfrb/Gucy1a1/Cox4i2'),
28:('unassigned','mixed_fibroblast_myeloid','Pdgfra/Col1a1 together with Csf3r/S100a8/S100a9'),
29:('other','plasma_cell','Jchain/Mzb1/Igha'),
30:('airway','neuroendocrine_candidate','Calca/Resp18/Pcsk1/Krt18'),
31:('other','cardiac_muscle','Myh6/Tnnt2/Tnni3/Actc1'),
32:('unassigned','mixed_airway_cytotoxic','Scgb1a1 with Ccl5/Gzma/Prf1; 11 cells from one library')}

def main():
    import numpy as np,pandas as pd,anndata as ad
    from scipy import sparse
    out=PAPER/'trials/u3_lineage_annotation';out.mkdir(exist_ok=True)
    cache=PAPER/'cache/u3_lineage_annotation';cache.mkdir(parents=True,exist_ok=True)
    if (out/'run_record.json').exists():raise RuntimeError('Annotation run exists')
    source=PAPER/'trials/u3_cluster_review'
    ann=pd.DataFrame([dict(cluster=str(k),compartment=v[0],label=v[1],rationale=v[2]) for k,v in ANNOTATIONS.items()]);ann.to_csv(out/'frozen_cluster_annotations.csv',index=False)
    spec=dict(frozen_utc=datetime.now(timezone.utc).isoformat(),basis='treatment-blind existing cluster top genes and detection fractions; independent broad annotation',KAC='unassigned: author classifier/labels not publicly available in audited source',alveolar_denominator='clusters 20 and 23, not a KAC classifier',macrophage_specific='cluster 6 candidate only; broad CSF1R+ myeloid is not macrophage-only',exclusions='mixed clusters 3,18,28,32 stay unassigned; muscle and ambiguous stromal not fibroblasts',annotation_sha256=sha256_file(out/'frozen_cluster_annotations.csv'),inputs={p.name:sha256_file(p) for p in [source/'cluster_top_genes.csv',source/'cluster_marker_summary.csv',PAPER/'cache/u3_qc/cluster_membership.csv.gz']})
    write_json_atomic(out/'specification.json',spec)
    state=dict(status='running',pid=os.getpid(),started_utc=spec['frozen_utc'],code=code_identity(ROOT,__file__));write_json_atomic(out/'run_record.json',state);start=time.monotonic()
    try:
        membership=pd.read_csv(PAPER/'cache/u3_qc/cluster_membership.csv.gz',index_col=0,dtype={'review_cluster':str})
        membership=membership.join(ann.set_index('cluster'),on='review_cluster',validate='many_to_one');assert membership.label.notna().all()
        run=json.loads((PAPER/'trials/u3_run_spec.json').read_text());treatment={x['gsm']:x['treatment'] for x in run['samples']}
        membership['treatment']=membership.gsm.map(treatment);membership.to_csv(cache/'cell_annotations.csv.gz',index_label='cell_id',compression='gzip')
        coverage=membership.groupby(['gsm','treatment','compartment','label'],observed=True).size().rename('cells').reset_index();coverage.to_csv(out/'sample_subtype_coverage.csv',index=False)
        tables=[];units=[];panels=[];genes=None;parity=[]
        for sample in run['samples']:
            gsm=sample['gsm'];a=ad.read_h5ad(PAPER/'cache/u3_qc'/(gsm+'_qc.h5ad'));a=a[a.obs.qc_pass].copy()
            assert a.obs_names.is_unique and set(a.obs_names)==set(membership[membership.gsm==gsm].index)
            m=membership.loc[a.obs_names];symbols=a.var.gene_symbol.astype(str).to_numpy()
            unique=sorted(set(symbols))
            if genes is None:genes=unique
            assert genes==unique,'Different assayed symbols require an explicit coverage crosswalk'
            gm={g:i for i,g in enumerate(genes)};projection=sparse.csr_matrix((np.ones(len(symbols),dtype=np.int32),(np.arange(len(symbols)),[gm[g] for g in symbols])),shape=(len(symbols),len(genes)))
            x=(a.X@projection).tocsr();assert np.all(x.data>=0) and np.all(x.data==np.round(x.data)) and x.sum()==a.X.sum()
            levels=sorted(m.label.unique());gidx=pd.Categorical(m.label,categories=levels).codes
            aggregate=sparse.csr_matrix((np.ones(len(m),dtype=np.int32),(gidx,np.arange(len(m)))),shape=(len(levels),len(m)))
            pb=(aggregate@x).toarray();assert pb.sum()==x.sum()
            for i,label in enumerate(levels):
                uid='u'+str(len(units));comp=m.loc[m.label==label,'compartment'].iloc[0]
                units.append(dict(unit_id=uid,gsm=gsm,treatment=treatment[gsm],compartment=comp,label=label,cells=int((m.label==label).sum())));tables.append(pb[i].astype(np.int64))
            wanted=[g for g in genes if g in {'Il1a','Il1b','Il1r1','Il1rap','Il1rn','Il1r2','Sigirr','Nlrp3','Pycard','Casp1','Gsdmd','Areg','Hbegf','Tgfb1','Ccl2','Cxcl12','Fgf7','Fgf10','Vegfa'}]
            totals=np.asarray(x.sum(axis=1)).ravel()
            for i,label in enumerate(levels):
                rows=np.flatnonzero(gidx==i);norm=x[rows][:,[gm[g] for g in wanted]].astype(float).multiply((10000/np.maximum(totals[rows],1))[:,None]).tocsr();det=(norm>0).mean(axis=0);mean=norm.mean(axis=0)
                for j,gene in enumerate(wanted):panels.append(dict(gsm=gsm,treatment=treatment[gsm],label=label,gene=gene,cells=len(rows),detection_fraction=float(det[0,j]),mean_normalized_10000=float(mean[0,j])))
            parity.append(dict(gsm=gsm,cells=len(a),counts=int(x.sum()),collapse_and_aggregation='passed'))
        pd.DataFrame(np.array(tables).T,index=genes,columns=[u['unit_id'] for u in units]).to_csv(cache/'subtype_counts.csv.gz',compression='gzip',index_label='gene')
        pd.DataFrame(units).to_csv(out/'subtype_units.csv',index=False);pd.DataFrame(panels).to_csv(out/'source_context_donor_values.csv',index=False);pd.DataFrame(parity).to_csv(out/'count_parity.csv',index=False)
        state.update(status='completed_broad_lineage_inputs_KAC_unresolved',cells=len(membership),genes=len(genes),pseudobulk_units=len(units))
    except Exception as exc:state.update(status='failed',error=str(exc));raise
    finally:state.update(elapsed_seconds=round(time.monotonic()-start,1),updated_utc=datetime.now(timezone.utc).isoformat());write_json_atomic(out/'run_record.json',state)
    print(json.dumps({k:v for k,v in state.items() if k!='code'}))

if __name__=='__main__':main()
