"""Bounded processed Visium acquisition; no raw reads or restricted objects."""
from pathlib import Path
import json,sys,re,time,urllib.request,urllib.parse
from datetime import datetime,timezone
PAPER=Path(__file__).resolve().parents[1];ROOT=PAPER.parents[1];sys.path.insert(0,str(ROOT))
from analysis.lib.provenance import write_json_atomic,sha256_file,code_identity

def main():
 out=PAPER/'trials/u5_spatial_context';cache=PAPER/'cache/u5_spatial_sources';out.mkdir(exist_ok=True);cache.mkdir(exist_ok=True)
 jobs=[]
 for series in ['GSE307534','GSE267226','GSE267228']:
  samples=json.loads((PAPER/f'trials/u0_geo_design_audit/{series}.json').read_text())['samples']
  for s in samples:
   title=s['title'][0]
   for url in s['supplementary_file']:
    if series=='GSE307534':
     assert url.endswith('.tar.gz');patient=int(re.search(r'patient (\d+)',title)[1]);hist='normal' if 'normal' in title.lower() else next(h for h in ['AAH','AIS','MIA','LUAD'] if h in title)
     filename=urllib.parse.unquote(url.split('/')[-1]);assert re.search(r'_P'+str(patient)+'_',filename)
     assert hist.lower() in filename.lower();extra=dict(patient=patient,histology=hist)
    elif '.h5' in url.lower():extra=dict(source_title=title)
    else:continue
    jobs.append(dict(series=series,gsm=s['accession'],title=title,url=url.replace('ftp://','https://'),**extra))
 spec=dict(frozen_utc=datetime.now(timezone.utc).isoformat(),jobs=jobs,per_file_cap=600_000_000,total_new_cap=15_000_000_000,selection='all deposited human precursor spatial sections; all post-viral count H5 files; no result-based selection',spatial_limit='independent pathology ROI annotations required for regional inference; source tissue mask allows descriptive whole-section measurements, not cell contacts',analysis_plan=dict(qc='in_tissue and UMI>=500 and detected_genes>=200 and mitochondrial_percentage<=15',whole_section='source-deposited lesion labels, patient-level paired descriptive differences; repeated same-histology sections averaged equally; spots not replicates',main_region_endpoint='gated on independent released spot-level pathology labels; do not derive ROIs from tested genes'))
 write_json_atomic(out/'acquisition_specification.json',spec)
 record=out/'acquisition_run_record.json';assert not record.exists();state=dict(status='running',started_utc=spec['frozen_utc'],code=code_identity(ROOT,__file__),files=[]);write_json_atomic(record,state);t=time.monotonic();total=0
 try:
  for job in jobs:
   name=urllib.parse.unquote(job['url'].split('/')[-1]);pilot=PAPER/'cache/extension_pilot'/name;target=pilot if pilot.exists() else cache/name
   if not target.exists():
    partial=target.with_name(target.name+'.part');n=0
    with urllib.request.urlopen(job['url'],timeout=90) as response,partial.open('wb') as dest:
     while True:
      chunk=response.read(1024*1024)
      if not chunk:break
      n+=len(chunk);total+=len(chunk);assert n<=spec['per_file_cap'] and total<=spec['total_new_cap'];dest.write(chunk)
    partial.replace(target)
   state['files'].append({**job,'path':str(target.relative_to(ROOT)),'bytes':target.stat().st_size,'sha256':sha256_file(target)})
   state.update(completed_files=len(state['files']),new_download_bytes=total,elapsed_seconds=round(time.monotonic()-t,1));write_json_atomic(record,state)
   print(job['gsm'],len(state['files']),len(jobs),flush=True)
  state['status']='completed'
 except Exception as e:state.update(status='failed',error=repr(e));raise
 finally:state.update(updated_utc=datetime.now(timezone.utc).isoformat(),elapsed_seconds=round(time.monotonic()-t,1));write_json_atomic(record,state)

if __name__=='__main__':main()
