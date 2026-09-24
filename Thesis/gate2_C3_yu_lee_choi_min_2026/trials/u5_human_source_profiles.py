"""All-QC human source/context profiles, retaining unassigned cells."""
from pathlib import Path
import sys,os,json,time,gc
from datetime import datetime,timezone
PAPER=Path(__file__).resolve().parents[1];ROOT=PAPER.parents[1];sys.path.insert(0,str(ROOT))
from analysis.lib.provenance import write_json_atomic,code_identity,sha256_file
from u5_human_aggregate import TRIAD

def main():
 import numpy as np,pandas as pd,anndata as ad
 out=PAPER/'trials/u5_human_sources';cache=PAPER/'cache/u5_human_sources';out.mkdir(exist_ok=True);cache.mkdir(exist_ok=True);record=out/'run_record.json';assert not record.exists()
 edges=pd.read_csv(PAPER/'trials/u5_ipf_compatibility/resource_edges.csv');panel=sorted({g for col in ['ligand','receptor'] for v in edges[col] for g in v.split('_')}|set('IL1RN IL1R2 SIGIRR IL1RAP NLRP3 PYCARD CASP1 GSDMD TNF IL6 SPP1 CSF1 CSF1R'.split()))
 jobs=json.loads((PAPER/'trials/u5_human_full/acquisition_specification.json').read_text())['jobs'];state=dict(status='running',pid=os.getpid(),started_utc=datetime.now(timezone.utc).isoformat(),code=code_identity(ROOT,__file__),completed=[]);write_json_atomic(record,state);t=time.monotonic();rows=[];coverage=[];parity=[]
 try:
  for job in jobs:
   gsm=job['gsm']
   while True:
    processing=json.loads((PAPER/'trials/u5_human_full/processing_run_record.json').read_text())
    if gsm in processing['completed_files']:break
    if processing['status']=='failed':raise RuntimeError('Upstream annotation failed')
    time.sleep(5)
   here=PAPER/'cache/u5_human_full'/gsm;a=ad.read_h5ad(here/'annotated_panel.h5ad');assayed=set(np.load(here/'unc20_pseudobulks.npz')['genes']);present=sorted(set(panel)&assayed&set(a.var_names));missing=sorted(set(panel)-set(present));coverage.append(dict(gsm=gsm,retained_measured_genes=len(present),unmeasured_or_not_retained=';'.join(missing)))
   x=a[:,present].X.tocsr();lib=a.obs.full_library_size.to_numpy();norm=x.astype(np.float64).multiply((10000/np.maximum(lib,1))[:,None]).tocsr()
   for cutoff in [.2,.3]:
    confident=a.obs.ann_finest_level_uncertainty.le(cutoff);labels=a.obs.ann_finest_level.astype(str).where(confident,'Unassigned');codes=pd.factorize(labels,sort=True);unique=codes[1];seen=0;group_counts=np.zeros(len(present),dtype=np.int64)
    for i,label in enumerate(unique):
     take=np.flatnonzero(codes[0]==i);n=len(take);seen+=n;counts=np.asarray(x[take].sum(axis=0)).ravel().astype(np.int64);detect=np.asarray((x[take]>0).sum(axis=0)).ravel();sumnorm=np.asarray(norm[take].sum(axis=0)).ravel();group_counts+=counts
     if label=='Unassigned':broad='Unassigned'
     else:
      comp=next((k for k,ls in TRIAD.items() if label in ls),None)
      if comp:broad={'AT2':'AT2-like','fibroblasts':'Fibroblasts','macrophages':'Macrophages'}[comp]
      elif a.obs.iloc[take].ann_level_2.eq('Myeloid').all():broad='Other myeloid'
      elif a.obs.iloc[take].ann_level_1.eq('Epithelial').all():broad='Other epithelial'
      else:broad='Other assigned'
     for j,g in enumerate(present):rows.append(dict(patient=job['patient'],histology=job['histology'],gsm=gsm,uncertainty=cutoff,label=label,broad=broad,gene=g,cells=n,count_sum=int(counts[j]),detected_cells=int(detect[j]),full_library_sum=int(lib[take].sum()),sum_normalized_10000=float(sumnorm[j])))
    assert seen==len(a) and np.array_equal(group_counts,np.asarray(x.sum(axis=0)).ravel());parity.append(dict(gsm=gsm,uncertainty=cutoff,all_QC_cells_assigned_or_unassigned=seen,exact_count_parity=True))
   state['completed'].append(gsm);write_json_atomic(record,state);del a,x,norm;gc.collect()
  raw=pd.DataFrame(rows);keys=['patient','histology','uncertainty','label','broad','gene'];agg=raw.groupby(keys,observed=True).agg(cells=('cells','sum'),count_sum=('count_sum','sum'),detected_cells=('detected_cells','sum'),full_library_sum=('full_library_sum','sum'),sum_normalized_10000=('sum_normalized_10000','sum'),libraries=('gsm','nunique')).reset_index();agg['detection_fraction']=agg.detected_cells/agg.cells;agg['mean_normalized_10000']=agg.sum_normalized_10000/agg.cells;agg.to_csv(cache/'all_QC_source_profiles.csv.gz',index=False,compression='gzip')
  contexts=agg[agg.gene.isin('IL1B IL1A IL1R1 IL1RAP IL1RN IL1R2 SIGIRR NLRP3 PYCARD CASP1 GSDMD'.split())];contexts.to_csv(out/'IL1_context_profiles.csv',index=False)
  source=agg[agg.gene=='IL1B'].groupby(['patient','histology','uncertainty','broad']).agg(count_sum=('count_sum','sum'),cells=('cells','sum'),detected_cells=('detected_cells','sum'),sum_normalized_10000=('sum_normalized_10000','sum')).reset_index();den=source.groupby(['patient','histology','uncertainty']).count_sum.transform('sum');source['fraction_of_observed_IL1B_counts']=source.count_sum/den.replace(0,np.nan);source['detection_fraction']=source.detected_cells/source.cells;source.to_csv(out/'IL1B_source_fractions.csv',index=False)
  pd.DataFrame(coverage).to_csv(out/'source_panel_coverage.csv',index=False);pd.DataFrame(parity).to_csv(out/'all_cell_count_parity.csv',index=False)
  state.update(status='completed',source_profile_rows=len(agg),validation='all QC cells represented exactly once per confidence cutoff, including unassigned; raw retained-gene count parity per library; no unassayed-gene zero interpretation')
 except Exception as e:state.update(status='failed',error=repr(e));raise
 finally:state.update(updated_utc=datetime.now(timezone.utc).isoformat(),elapsed_seconds=round(time.monotonic()-t,1));write_json_atomic(record,state)

if __name__=='__main__':main()
