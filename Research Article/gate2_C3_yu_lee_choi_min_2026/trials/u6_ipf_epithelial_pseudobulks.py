"""Source-labelled epithelial state inputs for the fibrosis specificity arm."""
from pathlib import Path
import sys,json,time,os
from datetime import datetime,timezone
PAPER=Path(__file__).resolve().parents[1];ROOT=PAPER.parents[1];sys.path.insert(0,str(ROOT))
from analysis.lib.provenance import write_json_atomic,sha256_file,code_identity
LABELS={'GSE136831':['ATI','ATII','Aberrant_Basaloid','Basal','Club','Ciliated'], 'GSE135893':['AT1','AT2','Transitional AT2','KRT5-/KRT17+','Basal','Proliferating Epithelial Cells','SCGB3A2+','SCGB3A2+ SCGB1A1+']}

def main():
 import numpy as np,pandas as pd,anndata as ad
 sys.path.insert(0,str(ROOT/'Research Article/gate1_01_niethamer_2025/trials'));import g2_gsea_ipf as loaders
 out=PAPER/'trials/u6_ipf_specificity';cache=PAPER/'cache/u6_ipf_specificity';out.mkdir(exist_ok=True);cache.mkdir(exist_ok=True)
 record=out/'input_run_record.json';assert not record.exists()
 spec=dict(frozen_utc=datetime.now(timezone.utc).isoformat(),labels=LABELS,selection='all source-labelled IPF/control epithelial cells in listed lineages; no score-based cell selection; counts retained before 30/50/100-cell floors',primary_contrasts={'GSE136831':['Aberrant_Basaloid','ATII'],'GSE135893':['KRT5-/KRT17+','AT2']},additional_state_contrast={'GSE135893':['Transitional AT2','AT2']},unit='donor; primary state contrasts within IPF donor; AT2 IPF/control descriptive contrast separate',identity_limit='deposited aberrant/transitional labels are not KAC, HPCS or malignancy',minimum_same_donor_pairs=3,source_modules_sha256=sha256_file(PAPER/'trials/u6_human_specificity/module_specification.json'))
 write_json_atomic(out/'specification.json',spec);state=dict(status='running',pid=os.getpid(),started_utc=spec['frozen_utc'],code=code_identity(ROOT,__file__),completed=[]);write_json_atomic(record,state);t=time.monotonic()
 class Inputs:
  def __init__(self):self.paths=[]
  def add_input(self,p):self.paths.append(Path(p))
 try:
  for cohort,labels in LABELS.items():
   inputs=Inputs();meta,genes,raw,_=(loaders.cohort_136831 if cohort=='GSE136831' else loaders.cohort_135893)(inputs);meta=meta.reset_index(drop=True);meta['disease']=meta.disease.replace({'Control':'control'})
   keep=meta.label.isin(labels)&meta.disease.isin(['IPF','control']);units=meta[keep].groupby(['donor','disease','label'],observed=True).size().rename('cells').reset_index();lookup={tuple(getattr(r,k) for k in ['donor','disease','label']):i for i,r in enumerate(units.itertuples(index=False))};codes=np.full(len(meta),-1,dtype=np.int64)
   for idx,r in meta.loc[keep,['donor','disease','label']].iterrows():codes[idx]=lookup[tuple(r)]
   state.update(cohort=cohort,stage='stream_full_raw_counts');write_json_atomic(record,state)
   counts,nnz=loaders.gu.pseudobulk_from_mtx(raw,codes,len(units));assert counts.shape==(len(genes),len(units));units['unit_id']=[cohort+'_'+str(i) for i in range(len(units))];units['full_library_sum']=counts.sum(axis=0)
   # Check all epithelial donor totals against the earlier independent stream.
   a=ad.read_h5ad(PAPER/'cache/u5_full_source_panels'/(cohort+'.h5ad'),backed='r');expected=a.obs.groupby(['donor','disease','label'],observed=True).full_library_size.sum();a.file.close()
   for r in units.itertuples():assert expected.loc[(r.donor,r.disease,r.label)]==r.full_library_sum
   units.to_csv(out/(cohort+'_units.csv'),index=False);pd.DataFrame(counts,index=genes,columns=units.unit_id).to_csv(cache/(cohort+'_counts.csv.gz'),compression='gzip',index_label='gene')
   write_json_atomic(out/(cohort+'_input_provenance.json'),dict(cohort=cohort,source_hashes={str(p.relative_to(ROOT)):sha256_file(p) for p in inputs.paths},matrix_entries=nnz,units=len(units),count_parity='all epithelial donor/label full-library totals match independent prior stream'))
   state['completed'].append(cohort);state.update(elapsed_seconds=round(time.monotonic()-t,1));write_json_atomic(record,state);print(cohort,'epithelial inputs complete',len(units),flush=True)
  state['status']='completed'
 except Exception as e:state.update(status='failed',error=repr(e));raise
 finally:state.update(updated_utc=datetime.now(timezone.utc).isoformat(),elapsed_seconds=round(time.monotonic()-t,1));write_json_atomic(record,state)

if __name__=='__main__':main()
