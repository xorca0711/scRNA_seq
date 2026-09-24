"""Paired patient RNA-compatibility contrasts and one-at-a-time sensitivities."""
from pathlib import Path
import sys,json,time,itertools
from datetime import datetime,timezone
PAPER=Path(__file__).resolve().parents[1];ROOT=PAPER.parents[1];sys.path.insert(0,str(ROOT))
from analysis.lib.provenance import write_json_atomic,code_identity
from u5_human_aggregate import CONTRASTS

def directions(ligand):
 lead=ligand.split('_')[0]
 if lead in ['IL1B','IL1A']:return [('macrophages','fibroblasts'),('macrophages','AT2'),('AT2','fibroblasts')]
 if lead in ['AREG','HBEGF']:return [('macrophages','fibroblasts'),('AT2','fibroblasts')]
 if lead=='TGFB1':return [('macrophages','fibroblasts')]
 if lead in ['CCL2','CXCL12']:return [('fibroblasts','macrophages')]
 if lead in ['FGF7','FGF10'] or lead.startswith('WNT'):return [('fibroblasts','AT2')]
 return []

def role(ligand,receptor):
 if ligand.split('_')[0] in ['IL1A','IL1B']:
  if any(g in receptor.split('_') for g in ['IL1R2','SIGIRR']):return 'inhibitory_or_decoy_context_not_activation'
  if {'IL1R1','IL1RAP'}<=set(receptor.split('_')):return 'canonical_signaling_complex_RNA_compatibility_only'
  return 'incomplete_canonical_complex_RNA_compatibility_only'
 return 'resource_RNA_compatibility_only'

def main():
 import numpy as np,pandas as pd
 out=PAPER/'trials/u5_human_niche';cache=PAPER/'cache/u5_human_niche';record=out/'compatibility_run_record.json';assert not record.exists()
 d=pd.read_csv(cache/'normalized_components.csv.gz');resources=pd.read_csv(PAPER/'trials/u5_ipf_compatibility/resource_edges.csv')
 state=dict(status='running',code=code_identity(ROOT,__file__),started_utc=datetime.now(timezone.utc).isoformat(),unit='paired patient; all source/receiver/histology combinations must be present',inference='descriptive mean paired change and leave-one-patient-out range; no cell-level P values');write_json_atomic(record,state);start=time.monotonic()
 summaries=[];eligibility=[];primary_values=[]
 try:
  for (config,floor,prior),g in d.groupby(['config','cell_floor','prior_count'],sort=False):
   lookup={key:part.set_index('gene') for key,part in g.groupby(['comp','label','patient','histology'])};labels=g[['comp','label']].drop_duplicates();patients=sorted(g.patient.unique())
   for view in ['broad','subtype']:
    labs=labels[labels.label.eq('__broad__') if view=='broad' else labels.label.ne('__broad__')]
    for edge in resources.itertuples():
     for source,target in directions(edge.ligand):
      for sl,tl in itertools.product(labs[labs.comp==source].label,labs[labs.comp==target].label):
       values={}
       for patient,hist in g[['patient','histology']].drop_duplicates().itertuples(index=False,name=None):
        a=lookup.get((source,sl,patient,hist));b=lookup.get((target,tl,patient,hist))
        if a is None or b is None:continue
        lig=edge.ligand.split('_');rec=edge.receptor.split('_')
        if not set(lig)<=set(a.index) or not set(rec)<=set(b.index):continue
        av=float(a.loc[lig,'logCPM'].min());bv=float(b.loc[rec,'logCPM'].min());values[(patient,hist)]=(av+bv)/2,av,bv
       for case,reference in CONTRASTS:
        complete=[p for p in patients if (p,case) in values and (p,reference) in values]
        keys=dict(config=config,cell_floor=floor,prior_count=prior,view=view,resource=edge.resource,family=edge.family,ligand=edge.ligand,receptor=edge.receptor,source=source,source_label=sl,target=target,target_label=tl,case=case,reference=reference,biological_role=role(edge.ligand,edge.receptor))
        eligibility.append({**keys,'complete_patients':len(complete),'eligible':len(complete)>=3})
        if len(complete)<3:continue
        differences=np.array([np.asarray(values[p,case])-values[p,reference] for p in complete]);effect=differences.mean(axis=0);loo=(differences.sum(axis=0)-differences)/(len(complete)-1)
        summaries.append({**keys,'n_patients':len(complete),'difference_case_minus_reference':effect[0],'ligand_component_difference':effect[1],'receptor_component_difference':effect[2],'loo_min':loo[:,0].min(),'loo_max':loo[:,0].max()})
        if config=='unc20_pooled' and floor==50 and prior==1 and edge.resource=='consensus':
         for p,delta in zip(complete,differences):primary_values.append({**keys,'patient':p,'case_score':values[p,case][0],'reference_score':values[p,reference][0],'paired_difference':delta[0],'ligand_difference':delta[1],'receptor_difference':delta[2]})
   print(config,floor,prior,'paired compatibility done',flush=True)
  result=pd.DataFrame(summaries);e=pd.DataFrame(eligibility);v=pd.DataFrame(primary_values)
  result.to_csv(out/'compatibility_contrasts.csv',index=False);e.to_csv(cache/'compatibility_eligibility.csv.gz',index=False,compression='gzip');v.to_csv(out/'primary_compatibility_patient_values.csv',index=False)
  primary=result[(result.config=='unc20_pooled')&(result.cell_floor==50)&(result.prior_count==1)&(result.resource=='consensus')];primary.to_csv(out/'primary_compatibility_contrasts.csv',index=False)
  # Independent group reconstruction checks every displayed primary contrast.
  keys=['view','resource','ligand','receptor','source','source_label','target','target_label','case','reference']
  index={k:x for k,x in v.groupby(keys)}
  for r in primary.itertuples():
   a=index[tuple(getattr(r,k) for k in keys)];assert a.patient.is_unique and len(a)==r.n_patients
   assert np.allclose(a.case_score-a.reference_score,a.paired_difference)
   assert np.isclose(a.paired_difference.mean(),r.difference_case_minus_reference)
   checks=[a.drop(i).paired_difference.mean() for i in a.index];assert np.allclose([min(checks),max(checks)],[r.loo_min,r.loo_max])
  state.update(status='completed',contrasts=len(result),primary_contrasts=len(primary),validation='each primary paired effect and omission range independently reconstructed')
 except Exception as e:state.update(status='failed',error=repr(e));raise
 finally:state.update(updated_utc=datetime.now(timezone.utc).isoformat(),elapsed_seconds=round(time.monotonic()-start,1));write_json_atomic(record,state)

if __name__=='__main__':main()
