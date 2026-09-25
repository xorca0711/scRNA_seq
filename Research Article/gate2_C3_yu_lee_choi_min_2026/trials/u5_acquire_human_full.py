"""Acquire all processed human libraries after pilot format/annotation checks."""
from pathlib import Path
import sys,json,urllib.request,hashlib,time,re
from datetime import datetime,timezone
PAPER=Path(__file__).resolve().parents[1];ROOT=PAPER.parents[1];sys.path.insert(0,str(ROOT))
from analysis.lib.provenance import write_json_atomic,sha256_file,code_identity

def main():
    out=PAPER/'trials/u5_human_full';out.mkdir(exist_ok=True)
    cache=PAPER/'cache/human_full_sources';cache.mkdir(parents=True,exist_ok=True)
    source=json.loads((PAPER/'trials/u0_geo_design_audit/GSE308103.json').read_text())
    jobs=[]
    for s in source['samples']:
        title=s['title'][0];patient=int(re.search(r'patient (\d+)',title)[1]);hist='normal' if 'normal' in title.lower() else next(h for h in ['AAH','AIS','MIA','LUAD'] if h in title)
        for url in s['supplementary_file']:
            assert url.endswith('.raw_counts.mtx.txt.gz')
            jobs.append(dict(gsm=s['accession'],patient=patient,histology=hist,url=url.replace('ftp://','https://')))
    spec=dict(frozen_utc=datetime.now(timezone.utc).isoformat(),series='GSE308103',selection='all 75 deposited processed libraries; no outcome-based selection; repeated tissues remain nested within 23 patients',per_file_cap=500_000_000,total_new_cap=8_000_000_000,jobs=jobs,pilot='eight libraries parsed and reference mapped; all series patient/title/file/demographic crosschecks agree',phenotype_limit='source-deposited histology labels; atlas candidates are not malignant/KAC labels')
    write_json_atomic(out/'acquisition_specification.json',spec)
    state=dict(status='running',started_utc=spec['frozen_utc'],code=code_identity(ROOT,__file__),files=[]);write_json_atomic(out/'acquisition_run_record.json',state);total=0;start=time.monotonic()
    try:
        for job in jobs:
            name=job['url'].rsplit('/',1)[1];pilot=PAPER/'cache/extension_pilot'/name;dest=pilot if pilot.exists() else cache/name
            if not dest.exists():
                temp=dest.with_suffix(dest.suffix+'.part');n=0
                with urllib.request.urlopen(job['url'],timeout=90) as r,temp.open('wb') as h:
                    expected=int(r.headers.get('Content-Length',0))
                    if expected>spec['per_file_cap'] or total+expected>spec['total_new_cap']:raise ValueError('Bounded processed-download cap exceeded')
                    while chunk:=r.read(1048576):
                        n+=len(chunk)
                        if n>spec['per_file_cap'] or total+n>spec['total_new_cap']:raise ValueError('Bounded processed-download stream cap exceeded')
                        h.write(chunk)
                    assert not expected or expected==n
                temp.replace(dest);total+=n
            state['files'].append(dict(**job,path=str(dest.relative_to(ROOT)),bytes=dest.stat().st_size,sha256=sha256_file(dest)))
            state.update(completed_files=len(state['files']),new_bytes=total,elapsed_seconds=round(time.monotonic()-start,1));write_json_atomic(out/'acquisition_run_record.json',state)
            print(job['gsm'],len(state['files']),len(jobs),flush=True)
        state['status']='completed'
    except Exception as exc:state.update(status='failed',error=str(exc));raise
    finally:state.update(updated_utc=datetime.now(timezone.utc).isoformat(),elapsed_seconds=round(time.monotonic()-start,1));write_json_atomic(out/'acquisition_run_record.json',state)

if __name__=='__main__':main()
