"""Pool technical/repeated tissues within patient and deposited histology."""
from pathlib import Path
import json,sys,time
from datetime import datetime,timezone
PAPER=Path(__file__).resolve().parents[1];ROOT=PAPER.parents[1];sys.path.insert(0,str(ROOT))
from analysis.lib.provenance import write_json_atomic,sha256_file,code_identity

TRIAD={
 'AT2':{'AT2','AT2 proliferating'},
 'fibroblasts':{'Adventitial fibroblasts','Alveolar fibroblasts','Peribronchial fibroblasts','Myofibroblasts'},
 'macrophages':{'Alveolar macrophages','Alveolar Mph proliferating','Monocyte-derived Mph','Interstitial Mph perivascular'},
}
CONTRASTS=[('AAH','normal'),('AIS','normal'),('MIA','normal'),('LUAD','normal'),('LUAD','AAH'),('LUAD','AIS'),('LUAD','MIA')]

def main():
 import numpy as np,pandas as pd
 out=PAPER/'trials/u5_human_niche';cache=PAPER/'cache/u5_human_niche';out.mkdir(exist_ok=True);cache.mkdir(exist_ok=True)
 record=out/'aggregation_run_record.json';assert not record.exists()
 processing=json.loads((PAPER/'trials/u5_human_full/processing_run_record.json').read_text());assert processing['status']=='completed_inputs_ready' and len(processing['completed_files'])==75
 q=pd.read_csv(PAPER/'trials/u5_human_full/library_qc_and_annotation.csv');assert q.gsm.is_unique and len(q)==75 and q.patient.nunique()==23
 chosen=q.sort_values(['qc_cells','gsm'],ascending=[False,True]).drop_duplicates(['patient','histology']);selected=set(chosen.gsm)
 chosen.to_csv(out/'largest_library_sensitivity_selection.csv',index=False)
 spec=dict(frozen_utc=datetime.now(timezone.utc).isoformat(),compartments={k:sorted(v) for k,v in TRIAD.items()},contrasts=[{'case':a,'reference':b} for a,b in CONTRASTS],unit='patient; same-histology libraries pooled; precursors never pooled together',primary=dict(uncertainty=.2,library_mode='pooled',cell_floor=50),sensitivities=dict(uncertainty=.3,library_mode='largest_QC_library_per_patient_histology',cell_floors=[30,100]),paired_minimum=3,assayed_fraction_minimum=.7,tested_pathway_genes_minimum=10,identity_limits='healthy-atlas candidate labels; low-confidence cells remain unassigned; no malignant/KAC/DATP assignment; assay background can affect RNA compatibility',selection_rule='largest QC cell library chosen without subtype or expression results; same chosen library used for every subtype',input_processing_record_sha256=sha256_file(PAPER/'trials/u5_human_full/processing_run_record.json'))
 write_json_atomic(out/'specification.json',spec)
 state=dict(status='running',code=code_identity(ROOT,__file__),started_utc=spec['frozen_utc'],completed=[]);write_json_atomic(record,state);t0=time.monotonic()
 try:
  for unc in [.2,.3]:
   tag=f'unc{int(100*unc):02d}';mats=[];detections=[];metadata=[];genes=None
   for gsm in q.gsm:
    here=PAPER/'cache/u5_human_full'/gsm;m=pd.read_csv(here/(tag+'_units.csv'));z=np.load(here/(tag+'_pseudobulks.npz'))
    if genes is None:genes=z['genes'];assert len(genes)==len(set(genes))
    assert np.array_equal(genes,z['genes']) and z['counts'].shape==(len(m),len(genes))
    assert np.array_equal(z['counts'].sum(axis=1),m.full_library_sum) and (z['detected_cells']<=m.cells.to_numpy()[:,None]).all()
    mats.append(z['counts']);detections.append(z['detected_cells']);metadata.append(m)
   full=np.concatenate(mats);detect=np.concatenate(detections);meta=pd.concat(metadata,ignore_index=True);del mats,detections
   for mode in ['pooled','largest_library']:
    if unc==.3 and mode!='pooled':continue # one-at-a-time sensitivities
    mask=np.ones(len(meta),dtype=bool) if mode=='pooled' else meta.gsm.isin(selected).to_numpy()
    m=meta.loc[mask].reset_index(drop=True);x=full[mask];d=detect[mask];new=[];rows=[];dr=[]
    for key,ix in m.groupby(['patient','histology','label'],sort=True).groups.items():
     g=m.loc[ix];s=x[list(ix)].sum(axis=0,dtype=np.int64);dd=d[list(ix)].sum(axis=0,dtype=np.int64)
     comp=next((k for k,v in TRIAD.items() if key[2] in v),'other')
     new.append(dict(patient=key[0],histology=key[1],label=key[2],comp=comp,view='subtype',cells=int(g.cells.sum()),full_library_sum=int(g.full_library_sum.sum()),libraries=';'.join(sorted(g.gsm.unique())),n_libraries=g.gsm.nunique(),uncertainty=unc,library_mode=mode))
     rows.append(s);dr.append(dd)
    a=pd.DataFrame(new);xx=np.asarray(rows);dd=np.asarray(dr);del rows,dr
    assert np.array_equal(xx.sum(axis=1),a.full_library_sum) and int(xx.sum())==int(x.sum())
    # Broad compartments explicitly exclude pericytes and other myeloid cells.
    broad=[];br=[];bd=[]
    for key,ix in a[a.comp!='other'].groupby(['patient','histology','comp'],sort=True).groups.items():
     g=a.loc[ix];broad.append(dict(patient=key[0],histology=key[1],comp=key[2],label='__broad__',view='broad',cells=int(g.cells.sum()),full_library_sum=int(g.full_library_sum.sum()),libraries=';'.join(sorted({s for v in g.libraries for s in v.split(';')})),n_libraries=len({s for v in g.libraries for s in v.split(';')}),uncertainty=unc,library_mode=mode));br.append(xx[list(ix)].sum(axis=0));bd.append(dd[list(ix)].sum(axis=0))
    a=pd.concat([a,pd.DataFrame(broad)],ignore_index=True);xx=np.concatenate([xx,np.asarray(br)]);dd=np.concatenate([dd,np.asarray(bd)])
    a['unit_id']=[f'{tag}_{mode}_{i:04d}' for i in range(len(a))];assert np.array_equal(xx.sum(axis=1),a.full_library_sum)
    name=tag+'_'+mode;a.to_csv(out/(name+'_units.csv'),index=False);np.savez_compressed(cache/(name+'_all_labels.npz'),genes=genes,counts=xx,detected_cells=dd)
    take=a.comp.ne('other').to_numpy();pd.DataFrame(xx[take].T,index=genes,columns=a.loc[take,'unit_id']).to_csv(cache/(name+'_triad_counts.csv.gz'),compression='gzip',index_label='gene');a.loc[take].to_csv(out/(name+'_triad_units.csv'),index=False)
    state['completed'].append(dict(name=name,all_units=len(a),triad_units=int(take.sum()),exact_count_parity=True));write_json_atomic(record,state);print(name,'units',len(a),flush=True)
   del full,detect,meta
  state['status']='completed'
 except Exception as e:state.update(status='failed',error=repr(e));raise
 finally:state.update(elapsed_seconds=round(time.monotonic()-t0,1),updated_utc=datetime.now(timezone.utc).isoformat());write_json_atomic(record,state)

if __name__=='__main__':main()
