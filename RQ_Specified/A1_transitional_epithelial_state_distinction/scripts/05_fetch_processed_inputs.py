"""Fetch selected processed GEO inputs with an explicit byte ceiling and hashes."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import time
import urllib.request
from a1_contract import ROOT,STUDY,digest


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--accessions',nargs='+',default=['GSE273123','GSE154966'])
    parser.add_argument('--max-file-mb',type=int,default=256)
    args=parser.parse_args(); records=[]
    for acc in args.accessions:
        catalog=json.loads((STUDY/'metadata'/f'{acc}.json').read_text(encoding='utf-8'))
        urls=[u for u in catalog['series'].get('supplementary_file',[]) if
              ('count' in u.lower() or u.endswith('.bed.gz')) and not u.endswith('.tar')]
        for source in urls:
            url=source.replace('ftp://','https://',1)
            path=STUDY/'cache/inputs'/acc/url.rsplit('/',1)[1]
            path.parent.mkdir(parents=True,exist_ok=True)
            t0=time.monotonic()
            if not path.exists():
                partial=path.with_suffix(path.suffix+'.partial')
                request=urllib.request.Request(url,headers={'User-Agent':'scRNA_seq-public-reanalysis/1.0'})
                with urllib.request.urlopen(request,timeout=90) as response,partial.open('wb') as handle:
                    length=int(response.headers.get('Content-Length',0))
                    if length>args.max_file_mb*1024**2:raise ValueError('File exceeds explicit ceiling')
                    size=0
                    while chunk:=response.read(1024**2):
                        size+=len(chunk)
                        if size>args.max_file_mb*1024**2:raise ValueError('Download exceeds explicit ceiling')
                        handle.write(chunk)
                partial.replace(path)
            item=dict(accession=acc,url=url,path=path.relative_to(ROOT).as_posix(),bytes=path.stat().st_size,
                      sha256=digest(path),checked_utc=datetime.now(timezone.utc).isoformat(),elapsed_seconds=round(time.monotonic()-t0,2))
            records.append(item);print(acc,path.name,item['bytes'],flush=True)
            out=STUDY/'reports/processed_input_inventory.json'
            old=json.loads(out.read_text())['files'] if out.exists() else []
            combined={r['path']:r for r in old};combined.update({r['path']:r for r in records})
            out.write_text(json.dumps(dict(scope='Downloaded processed assay inputs; sample/assay QC separate',files=list(combined.values())),indent=2)+'\n')


if __name__=='__main__':main()
