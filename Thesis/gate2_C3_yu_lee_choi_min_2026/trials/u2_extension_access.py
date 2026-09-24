"""Bounded public processed-file/header audit; no raw reads or restricted files."""
from pathlib import Path
import gzip
import hashlib
import json
import sys
from datetime import datetime,timezone
from concurrent.futures import ThreadPoolExecutor
import urllib.request
PAPER=Path(__file__).resolve().parents[1]
ROOT=PAPER.parents[1]
sys.path.insert(0,str(ROOT))
from analysis.lib.provenance import write_json_atomic,code_identity


def inspect(url):
    url=url.replace('ftp://','https://')
    row={'url':url}
    try:
        request=urllib.request.Request(url,method='HEAD',headers={'User-Agent':'research-data-audit/1.0'})
        with urllib.request.urlopen(request,timeout=45) as response:
            row.update(http_status=response.status,bytes=int(response.headers.get('Content-Length',0)),content_type=response.headers.get('Content-Type'))
        if '.raw_counts.mtx.txt.gz' in url:
            with urllib.request.urlopen(url,timeout=45) as response:
                compressed=response.read(131072)
            row['compressed_prefix_sha256']=hashlib.sha256(compressed).hexdigest()
            import zlib
            decoder=zlib.decompressobj(16+zlib.MAX_WBITS)
            content=decoder.decompress(compressed,262144).decode('utf-8',errors='replace')
            lines=content.splitlines()
            row['text_prefix']=[line[:700] for line in lines[:4]]
        row['status']='accessible'
    except Exception as exc:row.update(status='failed',error=str(exc))
    return row


def main():
    out=PAPER/'trials/u2_extension_access';out.mkdir(exist_ok=True)
    tasks=[]
    for series in ['GSE308103','GSE307534','GSE267226','GSE267228','GSE277777']:
        d=json.loads((PAPER/f'trials/u0_geo_design_audit/{series}.json').read_text())
        # Inspect at most two libraries per series; preserve the full metadata inventory separately.
        samples=d['samples'][:2]
        for sample in samples:
            for url in sample.get('supplementary_file',[]):
                tasks.append((series,sample['accession'],sample['title'],url))
    with ThreadPoolExecutor(max_workers=3) as pool:
        results=list(pool.map(lambda task:{'series':task[0],'gsm':task[1],'title':task[2],**inspect(task[3])},tasks))
    record={'checked_utc':datetime.now(timezone.utc).isoformat(),'status':'completed_access_audit','scope':'bounded HEAD requests and compressed text prefixes; not expression analysis','code':code_identity(ROOT,__file__),'files':results}
    write_json_atomic(out/'run_record.json',record)
    for row in results:print(json.dumps(row,ensure_ascii=True),flush=True)


if __name__=='__main__':main()
