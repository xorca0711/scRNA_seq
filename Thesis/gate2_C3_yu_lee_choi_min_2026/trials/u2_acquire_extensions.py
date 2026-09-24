"""Bounded public processed-data pilots; no sequencing or restricted files."""
from pathlib import Path
import hashlib,json,re,sys,urllib.request,time
from datetime import datetime,timezone
PAPER=Path(__file__).resolve().parents[1];ROOT=PAPER.parents[1]
sys.path.insert(0,str(ROOT))
from analysis.lib.provenance import write_json_atomic,code_identity


def main():
    out=PAPER/'trials/u2_extension_pilot';out.mkdir(exist_ok=True)
    cache=PAPER/'cache/extension_pilot';cache.mkdir(parents=True,exist_ok=True)
    d=json.loads((PAPER/'trials/u0_geo_design_audit/GSE308103.json').read_text());by={}
    for s in d['samples']:
        title=s['title'][0];p=re.search(r'patient (\d+)',title)
        hist='normal' if 'normal' in title.lower() else next((h for h in ['AAH','AIS','MIA','LUAD'] if h in title),None)
        if p and hist:by.setdefault(int(p[1]),{}).setdefault(hist,[]).append(s)
    patients=sorted(p for p,x in by.items() if 'normal' in x and 'AAH' in x)[:3]
    assert len(patients)==3
    jobs=[]
    for patient in patients:
        for hist in ['normal','AAH']:
            for s in by[patient][hist]:
                for url in s['supplementary_file']:jobs.append({'series':'GSE308103','gsm':s['accession'],'patient':patient,'histology':hist,'url':url.replace('ftp://','https://'),'role':'paired_nucleus_pilot'})
    for series in ['GSE307534','GSE277777']:
        d=json.loads((PAPER/f'trials/u0_geo_design_audit/{series}.json').read_text())
        for s in d['samples'][:(2 if series=='GSE307534' else 1)]:
            for url in s['supplementary_file']:
                if series=='GSE277777' and not url.endswith('.h5ad'):continue
                jobs.append({'series':series,'gsm':s['accession'],'url':url.replace('ftp://','https://'),'role':'spatial_pair_pilot' if series=='GSE307534' else 'annotated_HPCS_feasibility'})
    manifest={'frozen_utc':datetime.now(timezone.utc).isoformat(),'patients':patients,'per_file_limit':250_000_000,'total_new_limit':800_000_000,'selection':'first three numeric patient labels with normal and AAH; all repeats; first spatial pair; first annotated HPCS object','jobs':jobs}
    write_json_atomic(out/'specification.json',manifest)
    state={'status':'running','code':code_identity(ROOT,__file__),'files':[]};write_json_atomic(out/'run_record.json',state)
    total=0;start=time.monotonic()
    try:
        for job in jobs:
            name=job['url'].rsplit('/',1)[1];dest=cache/name
            if dest.exists():
                state['files'].append({**job,'path':str(dest.relative_to(ROOT)),'bytes':dest.stat().st_size,'sha256':hashlib.sha256(dest.read_bytes()).hexdigest(),'status':'existing'});continue
            with urllib.request.urlopen(job['url'],timeout=90) as response:
                expected=int(response.headers.get('Content-Length',0))
                if expected>manifest['per_file_limit'] or total+expected>manifest['total_new_limit']:raise RuntimeError('Declared download cap exceeded')
                h=hashlib.sha256();n=0;temp=dest.with_suffix(dest.suffix+'.part')
                with temp.open('wb') as handle:
                    while chunk:=response.read(1_048_576):
                        n+=len(chunk)
                        if n>manifest['per_file_limit'] or total+n>manifest['total_new_limit']:raise RuntimeError('Streaming cap exceeded')
                        handle.write(chunk);h.update(chunk)
                if expected and n!=expected:raise RuntimeError('Download byte count mismatch')
                temp.replace(dest);total+=n
            state['files'].append({**job,'path':str(dest.relative_to(ROOT)),'bytes':n,'sha256':h.hexdigest(),'status':'retrieved'})
            write_json_atomic(out/'run_record.json',state);print(job['gsm'],name,n,flush=True)
        state.update(status='completed',new_bytes=total)
    except Exception as exc:state.update(status='failed',error=str(exc));raise
    finally:
        state.update(updated_utc=datetime.now(timezone.utc).isoformat(),elapsed_seconds=round(time.monotonic()-start,1));write_json_atomic(out/'run_record.json',state)


if __name__=='__main__':main()
