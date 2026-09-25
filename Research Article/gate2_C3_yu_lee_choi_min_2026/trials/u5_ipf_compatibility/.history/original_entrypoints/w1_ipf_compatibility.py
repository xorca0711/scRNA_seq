"""Frozen donor-level RNA compatibility and robustness; no causal/P-value claim."""
from pathlib import Path
import itertools
import json
import subprocess
import sys
import time
from datetime import datetime, timezone
PAPER=Path(__file__).resolve().parents[1]
ROOT=PAPER.parents[1]
sys.path.insert(0,str(ROOT))
from analysis.lib.provenance import sha256_file,write_json_atomic,code_identity


def main():
    import numpy as np
    import pandas as pd
    out=PAPER/'trials/w1_ipf_compatibility';out.mkdir(exist_ok=True)
    if (out/'specification.json').exists():
        raise RuntimeError('Existing frozen run: inspect it instead of silently refreezing')
    sources=[PAPER/'trials/u5_ipf_liana'/('resource_'+r+'.csv') for r in ['consensus','cellchatdb']]
    resources=[]
    for path in sources:
        df=pd.read_csv(path);df['resource']=path.stem.removeprefix('resource_');resources.append(df)
    resources=pd.concat(resources,ignore_index=True)
    resources.to_csv(out/'resource_edges.csv',index=False)
    genes=sorted({g for col in ['ligand','receptor'] for entry in resources[col] for g in entry.split('_')}|{'IL1RAP','IL1RN','IL1R2','SIGIRR'})
    (out/'component_genes.txt').write_text('\n'.join(genes)+'\n')
    inputs=sources+[ROOT/'analysis/corrections/statistics/cache'/(c+s) for c in ['GSE136831','GSE135893'] for s in ['_subtypes_counts.csv.gz','_subtypes_units.csv']]
    spec={'frozen_utc':datetime.now(timezone.utc).isoformat(),'primary':{'resource':'consensus','cell_floor':50,'prior_count':1},
          'sensitivity':{'resource':'cellchatdb','cell_floor':[30,100],'prior_count':[0.5,2]},
          'units':'donor; broad and deposited-subtype views are separate, not independent replications',
          'normalization':'edgeR TMM on all assayed genes within subtype/cohort/floor using all eligible IPF/control donors; log2 CPM with specified scaled prior count',
          'score':'half sum of minimum ligand-subunit sender logCPM and minimum receptor-subunit receiver logCPM',
          'contrast':'mean IPF minus mean control, with individual donor values and leave-one-donor-out range',
          'inference':'descriptive only: disease-label exchangeability and sufficient covariate adjustment not established; no P values or FDR claims',
          'minimum_donors_per_arm':3,'cell_floor_limit':'30-cell sensitivity conditional on original cache broad-compartment >=50 cells',
          'selection':'resource-frozen families; intended source/receiver directions; no effect-based subtype selection',
          'input_hashes':{str(p.relative_to(ROOT)):sha256_file(p) for p in inputs}}
    write_json_atomic(out/'specification.json',spec)
    start=time.monotonic();state={'status':'running','started_utc':spec['frozen_utc'],'code':code_identity(ROOT,__file__),'completed_cohorts':[]}
    write_json_atomic(out/'run_record.json',state)
    try:
        rscript=ROOT/'analysis/corrections/statistics/.tools/R-portable/app/bin/Rscript.exe'
        summaries=[];donor_rows=[];eligible_rows=[]
        for cohort in ['GSE136831','GSE135893']:
            subprocess.run([str(rscript),str(PAPER/'trials/w1_normalize.R'),cohort],cwd=ROOT,check=True)
            df=pd.read_csv(PAPER/'cache/w1_ipf'/(cohort+'_components.csv.gz'))
            labels=df[['comp','label']].drop_duplicates()
            for floor,prior in itertools.product([50,30,100],[1,.5,2]):
                d=df[(df.cell_floor==floor)&(df.prior_count==prior)]
                lookup={key:part.set_index('gene') for key,part in d.groupby(['comp','label','donor'])}
                for view in ['broad','subtype']:
                    labs=labels[labels.label.eq('__broad__') if view=='broad' else labels.label.ne('__broad__')]
                    for _,edge in resources.iterrows():
                        ligand=edge.ligand.split('_');receptor=edge.receptor.split('_')
                        lead=ligand[0]
                        if lead in ['IL1B','IL1A']: directions=[('macrophages','fibroblasts'),('macrophages','AT2'),('AT2','fibroblasts')]
                        elif lead in ['AREG','HBEGF']: directions=[('macrophages','fibroblasts'),('AT2','fibroblasts')]
                        elif lead=='TGFB1': directions=[('macrophages','fibroblasts')]
                        elif lead in ['CCL2','CXCL12']: directions=[('fibroblasts','macrophages')]
                        elif lead in ['FGF7','FGF10']: directions=[('fibroblasts','AT2')]
                        else: continue  # neutrophil/endothelial recipients are outside these source inputs
                        for source,target in directions:
                            for sl,tl in itertools.product(labs[labs.comp==source].label,labs[labs.comp==target].label):
                                keys=dict(cohort=cohort,view=view,resource=edge.resource,family=edge.family,ligand=edge.ligand,receptor=edge.receptor,source=source,source_label=sl,target=target,target_label=tl,cell_floor=floor,prior_count=prior)
                                ds=sorted(set(d[(d.comp==source)&(d.label==sl)].donor)&set(d[(d.comp==target)&(d.label==tl)].donor))
                                vals=[];missing=set()
                                for donor in ds:
                                    a=lookup[(source,sl,donor)];b=lookup[(target,tl,donor)]
                                    missing.update(set(ligand)-set(a.index));missing.update(set(receptor)-set(b.index))
                                    if not set(ligand)<=set(a.index) or not set(receptor)<=set(b.index): continue
                                    disease=a.disease.iloc[0];assert disease==b.disease.iloc[0]
                                    av=float(a.loc[ligand,'logCPM'].min());bv=float(b.loc[receptor,'logCPM'].min())
                                    vals.append(dict(donor=donor,disease=disease,score=(av+bv)/2,ligand_min=av,receptor_min=bv,source_cells=int(a.cells.iloc[0]),target_cells=int(b.cells.iloc[0])))
                                v=pd.DataFrame(vals);n0=sum(x['disease']=='control' for x in vals);n1=sum(x['disease']=='IPF' for x in vals)
                                status='eligible_descriptive' if min(n0,n1)>=3 else ('missing_assay_subunits' if missing else 'insufficient_same_donor_coverage')
                                eligible_rows.append({**keys,'n_control':n0,'n_IPF':n1,'status':status,'missing_genes':';'.join(sorted(missing))})
                                donor_rows.extend({**keys,**x,'eligible_contrast':status=='eligible_descriptive'} for x in vals)
                                if status!='eligible_descriptive': continue
                                def delta(w,col='score'): return float(w.loc[w.disease=='IPF',col].mean()-w.loc[w.disease=='control',col].mean())
                                loo=[delta(v.drop(i)) for i in v.index]
                                summaries.append({**keys,'n_control':n0,'n_IPF':n1,'difference_IPF_minus_control':delta(v),'ligand_component_difference':delta(v,'ligand_min'),'receptor_component_difference':delta(v,'receptor_min'),'loo_min':min(loo),'loo_max':max(loo),'inference':'descriptive_no_p_value'})
            state['completed_cohorts'].append(cohort);write_json_atomic(out/'run_record.json',state)
        summary=pd.DataFrame(summaries);eligibility=pd.DataFrame(eligible_rows)
        summary.to_csv(out/'contrasts.csv',index=False);eligibility.to_csv(out/'eligibility.csv',index=False)
        pd.DataFrame(donor_rows).to_csv(PAPER/'cache/w1_ipf/donor_scores.csv.gz',index=False)
        primary=summary[(summary.resource=='consensus')&(summary.cell_floor==50)&(summary.prior_count==1)]
        primary.to_csv(out/'primary_contrasts.csv',index=False)
        assert np.isfinite(summary[['difference_IPF_minus_control','loo_min','loo_max']]).all().all()
        assert (summary[['n_control','n_IPF']]>=3).all().all()
        state.update(status='completed_descriptive_contrasts',elapsed_seconds=round(time.monotonic()-start,2),primary_contrasts=len(primary),sensitivity_and_primary_contrasts=len(summary),outputs={p.name:sha256_file(p) for p in out.glob('*.csv')})
    except Exception as exc:
        state.update(status='failed',error=str(exc));raise
    finally:
        state['updated_utc']=datetime.now(timezone.utc).isoformat();write_json_atomic(out/'run_record.json',state)
    print(json.dumps({k:v for k,v in state.items() if k not in ['code','outputs']}))


if __name__=='__main__':main()
