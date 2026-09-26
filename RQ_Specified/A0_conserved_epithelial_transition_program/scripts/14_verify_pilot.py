"""Independently verify saved discovery arithmetic and raw-count probe genes."""
from __future__ import annotations
import csv,gzip,hashlib,json,math,statistics
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd

BASE=Path(__file__).resolve().parents[1];OUT=BASE/'tables/pilot_v1';WORK=BASE/'processed/pilot_v1'

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
    return h.hexdigest()

def verify():
    cfg=json.loads((BASE/'config/pilot_v1.json').read_text());rules=cfg['discovery']
    run=json.loads((OUT/'discovery_run.json').read_text());prep=json.loads((OUT/'preparation.json').read_text())
    frozen=json.loads((OUT/'frozen_programme.json').read_text());checks=0
    assert run['config_sha256']==prep['config_sha256']==sha(BASE/'config/pilot_v1.json')
    assert run['script_sha256']==sha(BASE/'scripts/12_discover_pilot.py')
    assert prep['script_sha256']==sha(BASE/'scripts/11_prepare_pilot.py')
    for file,digest in run['outputs'].items():assert sha(OUT/file)==digest;checks+=1
    assert not run['V1_expression_effects_read'];checks+=1
    with gzip.open(OUT/'discovery_gene_effects.tsv.gz','rt') as f: rows=list(csv.DictReader(f,delimiter='\t'))
    keys=list(rows[0]);eligible=[];max_summary_error=0
    for row in rows:
        contrasts=[];passed=True
        for role in ['D1','D2']:
            passed &= float(row[f'{role}_intermediate_detection_fraction'])+1e-12>=rules['minimum_detected_unit_fraction']
            for endpoint in ['start','destination']:
                values=[float(row[k]) for k in keys if k.startswith(role+'_') and k.endswith('_vs_'+endpoint)]
                assert len(values)==prep['roles'][role]['units']
                median=statistics.median(values);fraction=sum(v>0 for v in values)/len(values)
                err=abs(median-float(row[f'{role}_vs_{endpoint}_median']));max_summary_error=max(max_summary_error,err)
                assert err<1e-9 and abs(fraction-float(row[f'{role}_vs_{endpoint}_positive_fraction']))<1e-9
                decision=median>=rules['minimum_median_log2CPM_difference_each_contrast'] and fraction>=rules['minimum_positive_unit_fraction']
                assert decision==(row[f'{role}_vs_{endpoint}_pass']=='True')
                passed &= decision;contrasts.append(median);checks+=3
        assert passed==(row['four_contrast_pass']=='True');checks+=1
        label=row['human'] in cfg['label_exclusions_human']
        operational=row['human'].startswith(('MT-','RPL','RPS')) or row['human']=='XIST'
        assert label==(row['label_excluded']=='True') and operational==(row['operational_excluded']=='True');checks+=2
        good=passed and not label and not operational
        assert good==(row['eligible']=='True') and abs(min(contrasts)-float(row['minimum_median_effect']))<1e-9;checks+=2
        if good:eligible.append(row)
    eligible.sort(key=lambda r:(-float(r['minimum_median_effect']),r['human']))
    assert len(eligible)==frozen['eligible_genes']==run['eligible_genes'];checks+=1
    selected=eligible[:rules['maximum_genes']] if len(eligible)>=rules['minimum_genes'] else []
    assert [r['human'] for r in selected]==frozen['human_genes'];checks+=1
    assert set(frozen['human_genes'])=={r['human'] for r in rows if r['selected']=='True'};checks+=1
    if not selected:
        assert frozen['status']=='STOP_NO_QUALIFYING_COMMON_MODULE'
        assert not (OUT/'transfer_run.json').exists();checks+=2
    # Independent scalar calculation on all chosen genes plus 32 evenly spaced
    # universe genes. No discovery helper or pseudobulk output is reused.
    probes={r['human']:r for r in rows[::max(1,len(rows)//32)]}
    probes.update({r['human']:r for r in selected or eligible})
    max_raw_error=0;raw_checks=0
    for role in ['D1','D2']:
        gene_names=json.loads((WORK/f'{role}_genes.json').read_text())
        index={g:i for i,g in enumerate(gene_names)}
        c=pd.read_csv(WORK/f'{role}_cells.tsv',sep='\t')
        x=np.memmap(WORK/f'{role}_counts.bin',mode='r',dtype=np.int32,shape=(len(gene_names),len(c)))
        for row in probes.values():
            vals=x[index[row['human' if role=='D2' else 'mouse']],:]
            for unit in sorted(c.unit.unique()):
                logs={}
                for state in ['start','intermediate','destination']:
                    ix=np.flatnonzero((c.unit==unit)&(c.state==state))
                    total=sum(int(c.library_total.iloc[i]) for i in ix)
                    count=sum(int(vals[i]) for i in ix)
                    logs[state]=math.log2(1+1e6*count/total)
                for endpoint in ['start','destination']:
                    err=abs((logs['intermediate']-logs[endpoint])-float(row[f'{role}_{unit}_vs_{endpoint}']))
                    max_raw_error=max(max_raw_error,err);assert err<1e-9;raw_checks+=1
        del x
    return dict(status='PASS',completed_utc=datetime.now(timezone.utc).isoformat(),verifier_sha256=sha(Path(__file__)),
                summary_checks=checks,raw_count_probe_genes=len(probes),raw_count_comparisons=raw_checks,
                maximum_summary_error=max_summary_error,maximum_raw_error=max_raw_error,
                source_raw_count_comparisons=prep['raw_count_comparisons'],programme_decision=frozen['status'],
                eligible_genes=frozen['eligible_genes'],selected_genes=frozen['selected_genes'])

if __name__=='__main__':
    dest=OUT/'verification.json'
    if dest.exists():raise SystemExit('Refusing to overwrite verification')
    record=verify();dest.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8');print(json.dumps(record))
