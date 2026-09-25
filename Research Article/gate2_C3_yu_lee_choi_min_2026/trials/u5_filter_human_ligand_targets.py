"""Restrict human paired target fits to same-patient expression support."""
from pathlib import Path
import sys,json
from datetime import datetime,timezone
PAPER=Path(__file__).resolve().parents[1];ROOT=PAPER.parents[1];sys.path.insert(0,str(ROOT))
from analysis.lib.provenance import write_json_atomic,sha256_file,code_identity
from u5_human_aggregate import TRIAD
from u5_liana_robustness import DIRECTIONS

def main():
 import numpy as np,pandas as pd
 out=PAPER/'trials/u5_human_ligand_targets';cache=PAPER/'cache/u5_human_niche';record=out/'filter_run_record.json';assert not record.exists()
 e=pd.read_csv(out/'target_eligibility.csv');units=pd.read_csv(out/'receiver_units.csv');prior_lr=pd.read_csv(out/'prior_ligand_receptor.csv')[['from','to']].drop_duplicates()
 state=dict(status='running',code=code_identity(ROOT,__file__),started_utc=datetime.now(timezone.utc).isoformat(),target_rule='>=10 mapped q<.05 DE targets, up/down separately; paired patient design and full patient-omission refits',expression_rule='same case-histology patient, fixed source label and receiver >=50 cells and every complex subunit >=10% detection, >=3 supporting patients; omission support rechecked',score='official unsigned prior Pearson statistic, descriptive ranks; down-target fit not causal inhibition');write_json_atomic(record,state)
 try:
  if not e.eligible.any():
   state.update(status='completed_no_eligible_target_sets',primary_candidate_rows=0);write_json_atomic(record,state);return
  scores=pd.read_csv(out/'all_prior_pearson_scores.csv.gz');meta=pd.read_csv(PAPER/'trials/u5_human_niche/unc20_pooled_units.csv');z=np.load(cache/'unc20_pooled_all_labels.npz');genes=z['genes'];detected=z['detected_cells'];assert len(meta)==len(detected)
  fixed=['IL1A','IL1B','AREG','HBEGF','TGFB1','CCL2','CXCL12','CXCL1','CXCL2','FGF7','FGF10','VEGFA'];ligands=sorted((set(fixed)|{g for g in genes if g.startswith('WNT')})&set(scores.ligand));measured=set(genes)
  pos={g:i for i,g in enumerate(genes)};lookup={(int(r.patient),r.histology,r.label,r.comp):i for i,r in enumerate(meta.itertuples(index=False))};source_labels=meta[meta.view=='subtype'][['label','comp']].drop_duplicates();comp_map={'AT2':'epithelial','fibroblasts':'fibroblast','macrophages':'myeloid'}
  receptor_options={}
  for ligand in ligands:
   rec=[r for r in prior_lr.loc[prior_lr['from']==ligand,'to'].unique() if r not in ['IL1R2','SIGIRR']]
   if ligand in ['IL1A','IL1B']:rec=['IL1R1_IL1RAP'] if set(rec)&{'IL1R1','IL1RAP'} else []
   if ligand=='TGFB1':rec=['TGFBR1_TGFBR2'] if set(rec)&{'TGFBR1','TGFBR2'} else []
   receptor_options[ligand]=rec
  coverage=[dict(ligand=l,receptor=r,assayed=all(g in measured for g in l.split('_')+r.split('_')),missing_genes=';'.join(g for g in l.split('_')+r.split('_') if g not in measured)) for l,rs in receptor_options.items() for r in rs];pd.DataFrame(coverage).to_csv(out/'prior_edge_assay_coverage.csv',index=False)
  ranked=[];support=[];combined=[];stability=[]
  for key,eg in e.groupby(['comp','label','case','reference']):
   comp,label,case,reference=key;primary=eg[eg.omitted_patient=='__none__'];u=units[(units.comp==comp)&(units.label==label)&(units.case==case)&(units.reference==reference)];patients=sorted(u.patient.unique());patient_sets={}
   if not primary.eligible.any():
    for row in primary.itertuples():combined.append(dict(comp=comp,label=label,case=case,reference=reference,direction=row.direction,scope='both',target_eligible=False,eligible_ligands=0,mapped_targets=row.mapped_targets))
    continue
   for ligand in ligands:
    for source in source_labels.itertuples(index=False):
     supported=[];matching=[]
     for patient in patients:
      si=lookup.get((int(patient),case,source.label,source.comp));ri=lookup.get((int(patient),case,label,comp))
      if si is None or ri is None or min(meta.cells.iloc[si],meta.cells.iloc[ri])<50:continue
      if ligand not in pos or detected[si,pos[ligand]]/meta.cells.iloc[si]<.1:continue
      good=[r for r in receptor_options[ligand] if all(g in pos and detected[ri,pos[g]]/meta.cells.iloc[ri]>=.1 for g in r.split('_'))]
      if good:supported.append(int(patient));matching.extend(good)
     patient_sets[(ligand,source.label,source.comp)]=set(supported);focused=(comp_map.get(source.comp),comp_map[comp]) in DIRECTIONS
     support.append(dict(comp=comp,label=label,case=case,reference=reference,ligand=ligand,source_label=source.label,source_comp=source.comp,focused_direction=focused,supporting_patients=len(supported),patients=';'.join(map(str,supported)),receptors=';'.join(sorted(set(matching)))))
   for scope in ['focused_triad','source_agnostic_planned_panel']:
    for row in primary.itertuples():
     def allowed(omit):
      return [l for l in ligands if any(len(ps-({int(omit)} if omit!='__none__' else set()))>=3 for (ll,sl,sc),ps in patient_sets.items() if ll==l and (scope!='focused_triad' or (comp_map.get(sc),comp_map[comp]) in DIRECTIONS))]
     full_allowed=allowed('__none__');combined.append(dict(comp=comp,label=label,case=case,reference=reference,direction=row.direction,scope=scope,target_eligible=bool(row.eligible),eligible_ligands=len(full_allowed) if row.eligible else 0,mapped_targets=row.mapped_targets))
     if not row.eligible:continue
     subset=scores[(scores.comp==comp)&(scores.label==label)&(scores.case==case)&(scores.reference==reference)&(scores.direction==row.direction)]
     for omitted,sg in subset.groupby('omitted_patient'):
      sg=sg[sg.ligand.isin(set(allowed(omitted))&set(full_allowed))&np.isfinite(sg.pearson)].copy();sg['scope']=scope;sg['rank']=sg.pearson.rank(ascending=False,method='min');sg['ranked_ligands']=len(sg);sg['planned_omissions']=len(patients);ranked.append(sg)
  r=pd.concat(ranked,ignore_index=True) if ranked else pd.DataFrame(columns=list(scores.columns)+['scope','rank','ranked_ligands','planned_omissions']);full=r[r.omitted_patient=='__none__'];full.to_csv(out/'primary_candidate_rankings.csv',index=False);r.to_csv(out/'candidate_rankings_with_omissions.csv',index=False)
  keys=['comp','label','case','reference','direction','scope','ligand']
  assert not r.duplicated(keys+['omitted_patient']).any()
  for key,g in r.groupby(keys):
   f=g[g.omitted_patient=='__none__'];loo=g[g.omitted_patient!='__none__']
   if f.empty:continue
   row=f.iloc[0];stability.append(dict(zip(keys,key))|dict(primary_pearson=row.pearson,primary_rank=row['rank'],planned_omissions=int(row.planned_omissions),eligible_omissions=len(loo),complete_omission_coverage=len(loo)==row.planned_omissions,min_rank=loo['rank'].min(),max_rank=loo['rank'].max(),minimum_pearson=loo.pearson.min(),maximum_pearson=loo.pearson.max()))
  pd.DataFrame(stability).to_csv(out/'candidate_ranking_stability.csv',index=False);pd.DataFrame(support).to_csv(out/'fixed_source_receiver_support.csv',index=False);pd.DataFrame(combined).to_csv(out/'combined_eligibility.csv',index=False)
  state.update(status='completed',primary_candidate_rows=len(full),ranking_rows_with_omissions=len(r),assay_genes=len(genes),assay_source='all assayed genes in exact patient/histology pseudobulks, not the targeted LR panel')
 except Exception as exc:state.update(status='failed',error=repr(exc));raise
 finally:state['updated_utc']=datetime.now(timezone.utc).isoformat();write_json_atomic(record,state)

if __name__=='__main__':main()
