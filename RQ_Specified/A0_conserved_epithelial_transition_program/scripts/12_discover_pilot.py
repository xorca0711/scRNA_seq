"""Run the fixed four-contrast discovery; freeze a module or stop without V1 scores."""
from __future__ import annotations
import hashlib,json,subprocess
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd

BASE=Path(__file__).resolve().parents[1]
OUT=BASE/'tables/pilot_v1'
WORK=BASE/'processed/pilot_v1'

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
    return h.hexdigest()

def load_role(role):
    genes=json.loads((WORK/f'{role}_genes.json').read_text())
    cells=pd.read_csv(WORK/f'{role}_cells.tsv',sep='\t')
    x=np.memmap(WORK/f'{role}_counts.bin',mode='r',dtype=np.int32,shape=(len(genes),len(cells)))
    return genes,cells,x

def main():
    cfg_path=BASE/'config/pilot_v1.json';cfg=json.loads(cfg_path.read_text())
    prep=json.loads((OUT/'preparation.json').read_text())
    assert sha(cfg_path)==prep['config_sha256']
    targets=['discovery_gene_effects.tsv','discovery_unit_effects.tsv','frozen_programme.json','discovery_run.json']
    if any((OUT/f).exists() for f in targets):raise SystemExit('Refusing to overwrite discovery')
    for name,record in prep['prepared'].items():assert sha(WORK/name)==record['sha256'],name
    ortho=pd.read_csv(OUT/'ortholog_universe.tsv',sep='\t')
    genes=ortho[['human','mouse']].copy(); gene_stats=[]; unit_rows=[]
    rules=cfg['discovery']; allpass=np.ones(len(genes),dtype=bool); min_effect=np.full(len(genes),np.inf)
    for role in ['D1','D2']:
        native,c,x=load_role(role); index={g:i for i,g in enumerate(native)}
        rows=np.array([index[g] for g in genes['human' if role=='D2' else 'mouse']])
        units=sorted(c.unit.unique());groups=[(u,s) for u in units for s in ['start','intermediate','destination']]
        masks=[np.flatnonzero((c.unit==u)&(c.state==s)) for u,s in groups]
        totals=np.array([c.iloc[ix].library_total.sum() for ix in masks],dtype=float)
        assert np.all(totals>0) and all(len(ix)>=cfg['cell_floor'] for ix in masks)
        cpm=np.empty((len(rows),len(groups)),dtype=float)
        for begin in range(0,len(rows),256):
            block=np.asarray(x[rows[begin:begin+256],:],dtype=np.int64)
            for j,ix in enumerate(masks): cpm[begin:begin+len(block),j]=block[:,ix].sum(axis=1)/totals[j]*1e6
        logs=np.log2(cpm+1).reshape(len(rows),len(units),3)
        det=(cpm.reshape(len(rows),len(units),3)[:,:,1]>=rules['minimum_intermediate_CPM']).mean(axis=1)
        genes[f'{role}_intermediate_detection_fraction']=det
        allpass &= det>=rules['minimum_detected_unit_fraction']
        for end_idx,end in [(0,'start'),(2,'destination')]:
            delta=logs[:,:,1]-logs[:,:,end_idx]
            median=np.median(delta,axis=1); frac=(delta>0).mean(axis=1)
            genes[f'{role}_vs_{end}_median']=median;genes[f'{role}_vs_{end}_positive_fraction']=frac
            passed=(median>=rules['minimum_median_log2CPM_difference_each_contrast'])&(frac>=rules['minimum_positive_unit_fraction'])
            genes[f'{role}_vs_{end}_pass']=passed;allpass &= passed; min_effect=np.minimum(min_effect,median)
            # Freeze all unit differences compactly, one gene per row, for independent verification.
            for j,unit in enumerate(units): genes[f'{role}_{unit}_vs_{end}']=delta[:,j]
        for j,unit in enumerate(units):
            for k,state in enumerate(['start','intermediate','destination']):
                unit_rows.append(dict(role=role,unit=unit,state=state,cells=len(masks[j*3+k]),pseudobulk_total=int(totals[j*3+k])))
        del x,cpm,logs
    exclude=set(cfg['label_exclusions_human'])
    genes['label_excluded']=genes.human.isin(exclude)
    genes['operational_excluded']=genes.human.str.startswith(('MT-','RPL','RPS'))|genes.human.eq('XIST')
    genes['four_contrast_pass']=allpass
    genes['eligible']=allpass&~genes.label_excluded&~genes.operational_excluded
    genes['minimum_median_effect']=min_effect
    candidates=genes[genes.eligible].sort_values(['minimum_median_effect','human'],ascending=[False,True])
    n=len(candidates); passed=n>=rules['minimum_genes']
    selected=candidates.head(rules['maximum_genes']) if passed else candidates.iloc[:0]
    genes['selected']=genes.human.isin(selected.human)
    genes.to_csv(OUT/'discovery_gene_effects.tsv',sep='\t',index=False,float_format='%.12g')
    pd.DataFrame(unit_rows).to_csv(OUT/'discovery_unit_effects.tsv',sep='\t',index=False)
    frozen=dict(status='FROZEN_FOR_TRANSFER' if passed else 'STOP_NO_QUALIFYING_COMMON_MODULE',eligible_genes=n,
                selected_genes=len(selected),human_genes=selected.human.tolist(),mouse_genes=selected.mouse.tolist(),
                unselected_eligible_human_genes=candidates.human.tolist() if not passed else [],
                decision_scope=cfg['adaptive_limits'],config_sha256=sha(cfg_path),ortholog_universe_sha256=sha(OUT/'ortholog_universe.tsv'))
    (OUT/'frozen_programme.json').write_text(json.dumps(frozen,indent=2)+'\n',encoding='utf-8')
    record=dict(completed_utc=datetime.now(timezone.utc).isoformat(),git_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=BASE,text=True).strip(),
                script_sha256=sha(Path(__file__)),config_sha256=sha(cfg_path),preparation_sha256=sha(OUT/'preparation.json'),
                source_programme_effects_previously_inspected=False,V1_expression_effects_read=False,
                status=frozen['status'],eligible_genes=n,selected_genes=len(selected),
                outputs={f:sha(OUT/f) for f in targets if f!='discovery_run.json'})
    (OUT/'discovery_run.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:record[k] for k in ['status','eligible_genes','selected_genes','V1_expression_effects_read']}))

if __name__=='__main__':main()
