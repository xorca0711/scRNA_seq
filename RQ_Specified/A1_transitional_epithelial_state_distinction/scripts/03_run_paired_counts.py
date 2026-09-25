"""Guarded bulk ATAC/RNA pilot. Default: print design holds, fit nothing.

Execution requires a frozen, approved contrast and a verified, adapted raw-count
TSV (feature ID first, distinct library columns next). Source-file adapters are
intentionally deferred until the deposited formats and identities are audited.
"""
import argparse
import csv
from datetime import datetime,timezone
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import sys

from a1_contract import ROOT,STUDY,digest,readiness
sys.path.insert(0,str(ROOT))
from analysis.lib.provenance import code_identity,write_json_atomic


def read_counts(path,selected):
    opener=gzip.open if path.suffix=='.gz' else open
    with opener(path,'rt',encoding='utf-8-sig',newline='') as handle:
        rows=csv.reader(handle,delimiter='\t')
        header=next(rows)
        if len(header)!=len(set(header)):raise ValueError('Duplicate count columns')
        indexes=[header.index(s['counts_column']) for s in selected]
        seen=set();matrix=[]
        for line in rows:
            if len(line)!=len(header) or not line[0] or line[0] in seen:
                raise ValueError('Invalid/duplicate feature row')
            seen.add(line[0])
            values=[int(line[i]) for i in indexes]
            if any(value<0 for value in values):raise ValueError('Negative count')
            matrix.append([line[0],*values])
    if not matrix:raise ValueError('Empty count matrix')
    if any(sum(row[j+1] for row in matrix)==0 for j in range(len(selected))):
        raise ValueError('Empty sample library')
    return matrix


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--contrast',required=True)
    parser.add_argument('--execute',action='store_true')
    parser.add_argument('--rscript',type=Path)
    args=parser.parse_args()
    contract_path=STUDY/'config/contrasts.json';sample_path=STUDY/'config/samples.json'
    contracts=json.loads(contract_path.read_text())
    contrast=next((c for c in contracts['contrasts'] if c['id']==args.contrast),None)
    if contrast is None:raise ValueError('Unknown contrast')
    reasons,selected=readiness(contrast,json.loads(sample_path.read_text())['samples'])
    if not contracts['execution_authorized']:reasons.append('Execution has not been authorized in the reviewed contract')
    if contrast['status']!='frozen':reasons.append('Contrast is not frozen')
    if not args.execute:
        print(json.dumps(dict(contrast=args.contrast,execution=False,holds=reasons),indent=2));return
    if reasons:raise ValueError('; '.join(reasons))
    if not args.rscript or not args.rscript.is_file():raise ValueError('Provide an existing Rscript executable')
    counts=read_counts(ROOT/contrast['counts_path'],selected)
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    out=STUDY/'reports/paired_counts'/args.contrast/stamp
    out.mkdir(parents=True,exist_ok=False)
    # Count payloads stay in ignored cache; compact results/provenance may be tracked.
    cache=STUDY/'cache/paired_counts'/args.contrast/stamp;cache.mkdir(parents=True,exist_ok=False)
    count_path=cache/'counts.tsv';design_path=out/'samples.tsv'
    with count_path.open('w',newline='',encoding='utf-8') as handle:
        writer=csv.writer(handle,delimiter='\t');writer.writerow(['feature_id',*[s['sample_id'] for s in selected]]);writer.writerows(counts)
    with design_path.open('w',newline='',encoding='utf-8') as handle:
        writer=csv.DictWriter(handle,fieldnames=['sample_id','biological_unit_id','group'],delimiter='\t')
        writer.writeheader();writer.writerows({k:s[k] for k in writer.fieldnames} for s in selected)
    r_code=STUDY/'scripts/04_fit_paired_counts.R'
    spec=dict(contrast=contrast,counts_path=str(count_path.resolve()),samples_path=str(design_path.resolve()),
              counts_sha256=digest(count_path),samples_sha256=digest(design_path),R_sha256=digest(r_code))
    write_json_atomic(out/'input_contract.json',spec)
    # Base R can read this companion without extra JSON/hash dependencies.
    parameters={k:contrast[k] for k in ['status','input_scale','reference','case',
                                      'minimum_independent_pairs','min_count','min_total_count']}
    parameters.update(counts_path=str(count_path.resolve()),samples_path=str(design_path.resolve()),
                      counts_md5=hashlib.md5(count_path.read_bytes()).hexdigest(),
                      samples_md5=hashlib.md5(design_path.read_bytes()).hexdigest())
    with (out/'input_contract.tsv').open('w',newline='',encoding='utf-8') as handle:
        writer=csv.writer(handle,delimiter='\t');writer.writerow(['key','value']);writer.writerows(parameters.items())
    record=dict(status='running',started_utc=stamp,code=code_identity(ROOT,__file__),
                input_contract_sha256=digest(out/'input_contract.json'),R_sha256=digest(r_code),
                R_parameters_sha256=digest(out/'input_contract.tsv'),
                contrasts_sha256=digest(contract_path),sample_manifest_sha256=digest(sample_path),
                original_count_sha256=contrast['counts_sha256'],analysis_scope='paired within-study bulk counts')
    write_json_atomic(out/'run_record.json',record)
    try:
        with (out/'R.log').open('w',encoding='utf-8') as log:
            subprocess.run([str(args.rscript.resolve()),str(r_code),str(out.resolve())],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,check=True)
        record['status']='completed'
    except Exception as error:
        record.update(status='failed',error=str(error));raise
    finally:
        record['finished_utc']=datetime.now(timezone.utc).isoformat()
        record['outputs']=[dict(path=p.relative_to(ROOT).as_posix(),sha256=digest(p)) for p in sorted(out.iterdir()) if p.is_file() and p.name!='run_record.json']
        write_json_atomic(out/'run_record.json',record)
    print(out.relative_to(ROOT).as_posix())


if __name__=='__main__':main()
