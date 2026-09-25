"""Bounded explicit annotation/identity evidence; downloaded code is never run."""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json
import urllib.request
import time

BASE=Path(__file__).resolve().parents[1]


def digest(data):return hashlib.sha256(data).hexdigest()


def main():
    rp=BASE/'reports/closure_source_inventory.json'
    if rp.exists():raise SystemExit('Refusing to overwrite closure inventory')
    tree_path=BASE/'cache/followup_sources/hpcs_code_tree.json'
    tree=json.loads(tree_path.read_text());items=[]
    for name in ['tracing_depletion_analysis/01_Concatenate_cells.ipynb','tracing_depletion_analysis/02_Analyze_control_DT_depletion.ipynb']:
        row=next(r for r in tree['tree'] if r['path']==name)
        items.append(dict(name='hpcs/'+name,url='https://raw.githubusercontent.com/dbetel/HPCS_LUAD/'+tree['sha']+'/'+name,
                          cap=80*1024**2,expected_bytes=row['size'],git_blob_sha=row['sha']))
    with (BASE/'tables/identity_audit/GEO_ENA_sample_audit.tsv').open() as h:
        for row in csv.DictReader(h,delimiter='\t'):
            exp=row['experiment']
            items.append(dict(name='sra/'+exp+'.xml',url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id='+exp+'&rettype=full&retmode=xml',cap=5*1024**2,accession=row['accession'],gsm=row['gsm']))
    records=[]
    for item in items:
        p=BASE/'cache/evidence_closure'/item['name']
        row=dict(item)
        try:
            if p.exists():data=p.read_bytes();row['cache_reused']=True
            else:
                req=urllib.request.Request(item['url'],headers={'User-Agent':'scRNA-seq evidence review'})
                with urllib.request.urlopen(req,timeout=90) as r:data=r.read(item['cap']+1)
                assert len(data)<=item['cap'],'Resource ceiling exceeded'
                if 'expected_bytes' in item:assert len(data)==item['expected_bytes']
                if 'git_blob_sha' in item:assert hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()==item['git_blob_sha']
                p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(data)
            if 'git_blob_sha' in item:assert hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()==item['git_blob_sha']
            row.update(status='retrieved',path=p.relative_to(BASE).as_posix(),bytes=len(data),sha256=digest(data))
            if p.suffix=='.ipynb':
                notebook=json.loads(data)
                source=[dict(cell=i,cell_type=c['cell_type'],source=''.join(c.get('source',[]))) for i,c in enumerate(notebook['cells'])]
                extracted=p.with_suffix('.source.json');extracted.write_text(json.dumps(source,indent=2)+'\n')
                row['extracted_source']=extracted.relative_to(BASE).as_posix();row['extracted_sha256']=digest(extracted.read_bytes())
        except Exception as exc:row.update(status='retrieval_failed',error=str(exc))
        records.append(row);print(item['name'],row['status'],row.get('bytes',row.get('error')),flush=True)
        if item['name'].startswith('sra/'):time.sleep(.4)
    rp.write_text(json.dumps(dict(utc=datetime.now(timezone.utc).isoformat(),code_sha256=digest(Path(__file__).read_bytes()),
        contract_sha256=digest((BASE/'config/evidence_closure_scope.md').read_bytes()),author_commit=tree['sha'],files=records),indent=2)+'\n')


if __name__=='__main__':main()
