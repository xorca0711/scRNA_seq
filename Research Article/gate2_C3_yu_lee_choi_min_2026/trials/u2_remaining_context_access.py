"""Bounded read-only audit of conditional processed-data/coordinate sources."""
from pathlib import Path
import json,sys,re,urllib.request,urllib.parse
from datetime import datetime,timezone
PAPER=Path(__file__).resolve().parents[1];ROOT=PAPER.parents[1];sys.path.insert(0,str(ROOT))
from analysis.lib.provenance import write_json_atomic,sha256_file,code_identity

def main():
 out=PAPER/'trials/u2_remaining_context_access';cache=PAPER/'cache/u2_remaining_context_access';out.mkdir(exist_ok=True);cache.mkdir(exist_ok=True)
 rows=[]
 def fetch(url,name,cap=5000000):
  path=cache/name
  try:
   with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'public-research-data-audit/1.0'}),timeout=45) as r:
    data=r.read(cap+1);assert len(data)<=cap;status=r.status
   path.write_bytes(data);rows.append(dict(url=url,file=str(path.relative_to(ROOT)),status=status,bytes=len(data),sha256=sha256_file(path)))
   return data.decode('utf-8')
  except Exception as e:rows.append(dict(url=url,status='unavailable',error=repr(e)));return None
 repo='MinhoLee-DGU/2023.Han.et.al.Nature';base='https://api.github.com/repos/'+repo
 commit=fetch(base+'/commits/main','han_commit.json');sha=json.loads(commit)['sha'] if commit else None
 tree=fetch(base+'/git/trees/'+sha+'?recursive=1','han_tree.json') if sha else None
 if tree:
  files=[x for x in json.loads(tree)['tree'] if x['type']=='blob'];write_json_atomic(out/'han_repository_file_inventory.json',dict(commit=sha,files=files))
  for name in ['README.md','scRNA/Whole clusters.R','scRNA/Subset epithelial.R','scRNA/AT1 Pseudobulk.R','scRNA/Other models and NDUFS2 cKO Epithelium.R']:
   if any(x['path']==name for x in files):fetch('https://raw.githubusercontent.com/'+repo+'/'+sha+'/'+urllib.parse.quote(name),re.sub('[^A-Za-z0-9_.]','_',name))
 fetch(base+'/issues?state=all&per_page=20','han_issues.json')
 for series in ['GSE267226','GSE267228','GSE222901','GSE300293','GSE307534']:
  text=fetch(f'https://ftp.ncbi.nlm.nih.gov/geo/series/{series[:-3]}nnn/{series}/suppl/',series+'_supplementary_directory.html')
  if text:
   links=re.findall(r'href="([^"#?]+)"',text);write_json_atomic(out/(series+'_series_files.json'),dict(series=series,links=[x for x in links if x not in ['../']]))
  audit=json.loads((PAPER/f'trials/u0_geo_design_audit/{series}.json').read_text());files=[]
  for s in audit['samples']:
   for url in s.get('supplementary_file',[]):
    if url!='NONE':files.append(dict(gsm=s['accession'],title=s['title'],url=url))
  write_json_atomic(out/(series+'_sample_files.json'),files)
 write_json_atomic(out/'run_record.json',dict(status='completed_bounded_access_audit',checked_utc=datetime.now(timezone.utc).isoformat(),code=code_identity(ROOT,__file__),files=rows,scope='public code/directory inventory only; no raw FASTQ acquisition or expression-completion claim'))
 print(json.dumps(dict(requests=len(rows),unavailable=sum(x['status']=='unavailable' for x in rows),han_commit=sha)))

if __name__=='__main__':main()
