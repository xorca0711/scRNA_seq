"""Native within-patient/histology LR measurements with no cell P values."""
from pathlib import Path
import sys,os,json,time,hashlib,warnings,gc
from datetime import datetime,timezone
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ.setdefault(k,'2')
PAPER=Path(__file__).resolve().parents[1];ROOT=PAPER.parents[1];sys.path.insert(0,str(ROOT))
from analysis.lib.provenance import write_json_atomic,sha256_file,code_identity
from u5_human_aggregate import TRIAD
from u5_human_compatibility import role
from u5_liana_robustness import CONFIGS,DIRECTIONS
from u5_liana_exact_components import normalized_components

def main():
 import numpy as np,pandas as pd,anndata as ad
 from liana.method import cellchat
 warnings.filterwarnings('ignore',category=FutureWarning,module='liana');warnings.filterwarnings('ignore',category=ad.ImplicitModificationWarning)
 out=PAPER/'trials/u5_human_niche';out.mkdir(exist_ok=True);cache=PAPER/'cache/u5_human_niche/liana';cache.mkdir(parents=True,exist_ok=True);record=out/'liana_run_record.json';assert not record.exists()
 resource=pd.read_csv(PAPER/'trials/u5_ipf_compatibility/resource_edges.csv');panel=sorted({g for col in ['ligand','receptor'] for x in resource[col] for g in x.split('_')})
 jobs=pd.DataFrame(json.loads((PAPER/'trials/u5_human_full/acquisition_specification.json').read_text())['jobs']);labels={v:{'AT2':'epithelial','fibroblasts':'fibroblast','macrophages':'myeloid'}[k] for k,vs in TRIAD.items() for v in vs}
 configs=[{**c,'uncertainty':.2,'library_mode':'pooled'} for c in CONFIGS]+[dict(name='uncertainty030',floor=50,expr=.1,cap=None,seed=20260924,uncertainty=.3,library_mode='pooled'),dict(name='largest_library',floor=50,expr=.1,cap=None,seed=20260924,uncertainty=.2,library_mode='largest_library')]
 spec=dict(frozen_utc=datetime.now(timezone.utc).isoformat(),configs=configs,subtype_compartments=labels,resource_sha256=sha256_file(PAPER/'trials/u5_ipf_compatibility/resource_edges.csv'),normalization='full-assay totals and observed global maximum; validated exact LIANA sufficient-component input',unit='patient x deposited histology; repeated same-histology libraries pooled; only source/receiver within same tissue group',identity_limit='reference-compatible candidate labels; unknown/malignant populations not forced into atlas cell types',RNA_limit='fixed RNA background may affect scores; no secretion/activation inference')
 write_json_atomic(out/'liana_specification.json',spec);state=dict(status='running',pid=os.getpid(),started_utc=spec['frozen_utc'],code=code_identity(ROOT,__file__),completed_groups=[],completed_calls=0);write_json_atomic(record,state);start=time.monotonic();all_calls=[];all_pairs=[]
 try:
  for (patient,histology),planned in jobs.groupby(['patient','histology'],sort=True):
   while True:
    processing=json.loads((PAPER/'trials/u5_human_full/processing_run_record.json').read_text())
    if set(planned.gsm)<=set(processing['completed_files']):break
    if processing['status']=='failed':raise RuntimeError('Human annotation failed')
    time.sleep(5)
   parts=[];library_cells={}
   for gsm in planned.gsm:
    a=ad.read_h5ad(PAPER/'cache/u5_human_full'/gsm/'annotated_panel.h5ad');library_cells[gsm]=len(a)
    keep=a.obs.ann_finest_level.isin(labels)&a.obs.ann_finest_level_uncertainty.le(.3);parts.append(a[keep,[g for g in panel if g in a.var_names]].copy());del a
   a=ad.concat(parts,join='inner');del parts;assert a.obs_names.is_unique;a.obs['celltype']=a.obs.ann_finest_level.astype(str);a.obs['compartment']=a.obs.celltype.map(labels)
   largest=sorted(library_cells,key=lambda x:(-library_cells[x],x))[0];scored_groups=[]
   for config in configs:
    part=a[a.obs.ann_finest_level_uncertainty.le(config['uncertainty']) & (True if config['library_mode']=='pooled' else a.obs.gsm.eq(largest))].copy()
    coverage=part.obs.groupby(['celltype','compartment'],observed=True).size().rename('cells').reset_index();good=coverage[coverage.cells>=config['floor']]
    pairs=pd.DataFrame([dict(source=r.celltype,target=s.celltype) for r in good.itertuples() for s in good.itertuples() if (r.compartment,s.compartment) in DIRECTIONS])
    if pairs.empty:
     for res in ['consensus','cellchatdb']:all_calls.append(dict(patient=patient,histology=histology,config=config['name'],resource=res,cells=len(part),pairs=0,retained_edges=0,status='no_eligible_pair'))
     continue
    pairs=pairs.drop_duplicates();selected=[]
    for label in good.celltype:
     indices=np.flatnonzero(part.obs.celltype.eq(label));seed=int.from_bytes(hashlib.sha256(f'{patient}|{histology}|{label}|{config["seed"]}'.encode()).digest()[:8],'little');rng=np.random.default_rng(seed)
     if config['cap'] and len(indices)>config['cap']:indices=rng.choice(indices,config['cap'],replace=False)
     selected.extend(indices)
    z=normalized_components(part[np.sort(selected)]);props={}
    for label in good.celltype:
     idx=np.flatnonzero(z.obs.celltype.eq(label));props[label]=pd.Series(np.asarray((z.X[idx]>0).sum(axis=0)).ravel()/len(idx),index=z.var_names)
    for p in pairs.itertuples(index=False):all_pairs.append(dict(patient=patient,histology=histology,config=config['name'],source=p.source,target=p.target))
    for res in ['consensus','cellchatdb']:
     table=resource[resource.resource==res];eligible=[]
     for p in pairs.itertuples(index=False):
      for edge in table.itertuples(index=False):
       vals=[props[p.source].get(g,-1) for g in edge.ligand.split('_')]+[props[p.target].get(g,-1) for g in edge.receptor.split('_')]
       if min(vals)>=config['expr']:eligible.append((p.source,p.target,edge.ligand,edge.receptor))
     if eligible:
      scored=cellchat(z,groupby='celltype',resource=table[['ligand','receptor']],groupby_pairs=pairs,expr_prop=config['expr'],min_cells=config['floor'],use_raw=False,n_perms=None,verbose=False,inplace=False)
      assert set(scored[['source','target','ligand_complex','receptor_complex']].itertuples(index=False,name=None))==set(eligible)
      scored['biological_role']=[role(l,r) for l,r in zip(scored.ligand_complex,scored.receptor_complex)]
      for k,v in dict(patient=patient,histology=histology,config=config['name'],resource=res).items():scored[k]=v
      assert scored.lr_probs.between(0,1).all() and not any('pval' in c.lower() for c in scored.columns);scored_groups.append(scored)
     all_calls.append(dict(patient=patient,histology=histology,config=config['name'],resource=res,cells=z.n_obs,pairs=len(pairs),retained_edges=len(eligible),status='completed' if eligible else 'no_expression_eligible_edge'));state['completed_calls']+=1
   if scored_groups:
    result=pd.concat(scored_groups,ignore_index=True);assert not result.duplicated(['config','resource','source','target','ligand_complex','receptor_complex']).any();result.to_csv(cache/f'P{patient}_{histology}_scores.csv.gz',index=False,compression='gzip')
   pd.DataFrame(all_calls).to_csv(out/'liana_call_coverage.csv',index=False);pd.DataFrame(all_pairs).to_csv(out/'liana_fixed_pairs.csv',index=False)
   state['completed_groups'].append(f'P{patient}_{histology}');state['elapsed_seconds']=round(time.monotonic()-start,1);write_json_atomic(record,state);print('Human LR',patient,histology,'done',flush=True);del a;gc.collect()
  state.update(status='completed',validation='native edge sets equal independent complex-subunit detection eligibility; unique bounded scores; no P values')
 except Exception as e:state.update(status='failed',error=repr(e));raise
 finally:state.update(updated_utc=datetime.now(timezone.utc).isoformat(),elapsed_seconds=round(time.monotonic()-start,1));write_json_atomic(record,state)

if __name__=='__main__':main()
