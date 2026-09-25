"""Descriptive early-treatment RNA-compatibility contrasts, fixed directions."""
from pathlib import Path
import sys,json,time
from datetime import datetime,timezone
PAPER=Path(__file__).resolve().parents[1];ROOT=PAPER.parents[1];sys.path.insert(0,str(ROOT))
from analysis.lib.provenance import write_json_atomic,code_identity,sha256_file

def directions(ligand):
    if ligand in ['IL1A','IL1B']:return [('myeloid','fibroblast'),('myeloid','alveolar'),('alveolar','fibroblast')]
    if ligand in ['AREG','HBEGF']:return [('myeloid','fibroblast'),('alveolar','fibroblast')]
    if ligand=='TGFB1':return [('myeloid','fibroblast')]
    if ligand in ['CCL2','CXCL12']:return [('fibroblast','myeloid')]
    if ligand in ['FGF7','FGF10'] or ligand.startswith('WNT'):return [('fibroblast','alveolar')]
    return []

def main():
    import numpy as np,pandas as pd
    out=PAPER/'trials/u4_mouse_niche';cache=PAPER/'cache/u4_mouse_niche'
    source=cache/'normalized_components.csv.gz';d=pd.read_csv(source);res=pd.read_csv(PAPER/'trials/u4_resources/mouse_resource_edges.csv')
    state=dict(status='running',started_utc=datetime.now(timezone.utc).isoformat(),code=code_identity(ROOT,__file__),components_sha256=sha256_file(source),scope='descriptive independent broad-lineage early mouse contrast; no KAC endpoint or causal claim');start=time.monotonic();write_json_atomic(out/'compatibility_run_record.json',state)
    summaries=[];eligibility=[];donors=[]
    try:
        for (floor,prior),df in d.groupby(['cell_floor','prior_count']):
            for view in ['broad','subtype']:
                s=df[df.label.eq('__broad__') if view=='broad' else df.label.ne('__broad__')]
                lookup={k:g.set_index('gene') for k,g in s.groupby(['compartment','label','gsm'])}
                for edge in res.itertuples():
                    for src,tgt in directions(edge.human_ligand):
                        for sl in s.loc[s.compartment==src,'label'].unique():
                            for tl in s.loc[s.compartment==tgt,'label'].unique():
                                keys=dict(resource=edge.resource,family=edge.family,view=view,ligand=edge.ligand,receptor=edge.receptor,human_ligand=edge.human_ligand,human_receptor=edge.human_receptor,source=src,source_label=sl,target=tgt,target_label=tl,cell_floor=floor,prior_count=prior)
                                gsms=sorted(set(s[(s.compartment==src)&(s.label==sl)].gsm)&set(s[(s.compartment==tgt)&(s.label==tl)].gsm));values=[];missing=set()
                                for gsm in gsms:
                                    a=lookup[(src,sl,gsm)];b=lookup[(tgt,tl,gsm)];lg=edge.ligand.split('_');rg=edge.receptor.split('_');missing.update(set(lg)-set(a.index));missing.update(set(rg)-set(b.index))
                                    if not set(lg)<=set(a.index) or not set(rg)<=set(b.index):continue
                                    assert a.treatment.iloc[0]==b.treatment.iloc[0]
                                    x=float(a.loc[lg,'logCPM'].min());y=float(b.loc[rg,'logCPM'].min());values.append(dict(gsm=gsm,treatment=a.treatment.iloc[0],ligand_min=x,receptor_min=y,score=(x+y)/2))
                                n0=sum(v['treatment']=='Control IgG' for v in values);n1=len(values)-n0;status='eligible_descriptive' if min(n0,n1)>=3 else ('missing_components' if missing else 'insufficient_same_animal_coverage')
                                eligibility.append(dict(**keys,n_control=n0,n_antiIL1B=n1,status=status,missing_components=';'.join(sorted(missing))))
                                if status!='eligible_descriptive':continue
                                z=pd.DataFrame(values)
                                def delta(v,col='score'):return float(v.loc[v.treatment!='Control IgG',col].mean()-v.loc[v.treatment=='Control IgG',col].mean())
                                loo=[delta(z.drop(i)) for i in z.index]
                                summaries.append(dict(**keys,n_control=n0,n_antiIL1B=n1,difference_antiIL1B_minus_IgG=delta(z),ligand_component_difference=delta(z,'ligand_min'),receptor_component_difference=delta(z,'receptor_min'),loo_min=min(loo),loo_max=max(loo),inference='descriptive_no_p_value'))
                                if floor==50 and prior==1:donors.extend(dict(**keys,**v) for v in values)
        results=pd.DataFrame(summaries);results.to_csv(out/'compatibility_contrasts.csv',index=False)
        primary=results[(results.cell_floor==50)&(results.prior_count==1)&(results.resource=='consensus')];primary.to_csv(out/'primary_compatibility_contrasts.csv',index=False)
        pd.DataFrame(eligibility).to_csv(cache/'compatibility_eligibility.csv.gz',index=False,compression='gzip');pd.DataFrame(donors).to_csv(out/'primary_compatibility_donor_values.csv',index=False)
        assert np.isfinite(results.difference_antiIL1B_minus_IgG).all()
        state.update(status='completed',primary_contrasts=len(primary),all_contrasts=len(results),method='half-sum of minimum ligand/receptor logCPM; all-gene TMM; 30/50/100 floors and 0.5/1/2 priors; no P values')
    except Exception as exc:state.update(status='failed',error=str(exc));raise
    finally:state.update(elapsed_seconds=round(time.monotonic()-start,1),updated_utc=datetime.now(timezone.utc).isoformat());write_json_atomic(out/'compatibility_run_record.json',state)
    print('Mouse RNA compatibility completed',len(primary),'primary rows')

if __name__=='__main__':main()
