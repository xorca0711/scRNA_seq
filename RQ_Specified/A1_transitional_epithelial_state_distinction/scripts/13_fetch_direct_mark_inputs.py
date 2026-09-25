"""Bounded public processed-track acquisition; independent of original inventories."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import datetime,timezone
import argparse,hashlib,json,time,urllib.request,subprocess,ssl,urllib.error

BASE=Path(__file__).resolve().parents[1]
def sha(p):
    with p.open('rb') as h:return hashlib.file_digest(h,'sha256').hexdigest()
def candidates():
    rows=[]
    for acc in ['GSE289683','GSE291333','GSE141635','GSE150527']:
        d=json.loads((BASE/'metadata'/f'{acc}.json').read_text(encoding='utf-8'))
        if acc=='GSE141635':
            for u in d['series']['supplementary_file']:
                if any(mark+'.bedgraph.gz' in u for mark in ['H3K4me3','H3K27ac','H3K36me3']):
                    rows.append(dict(accession=acc,gsm=None,title=Path(u).name,url=u.replace('ftp://','https://',1),assembly='mm10'))
        else:
            for rec in d['samples']:
                f=rec.get('fields',rec)
                for k,values in f.items():
                    if not k.startswith('supplementary_file'):continue
                    for u in values:
                        if acc=='GSE150527' and not u.endswith('.tab.gz'):continue
                        rows.append(dict(accession=acc,gsm=f['geo_accession'][0],title=f['title'][0],
                            url=u.replace('ftp://','https://',1),assembly='hg19' if acc=='GSE150527' else 'T2T-CHM13v2.0'))
    references={
      'T2T-CHM13v2.0':'https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_RefSeq_Liftoff_v5.3.gff.gz',
      'mm10':'https://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/refGene.txt.gz',
      'hg19':'https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz',
    }
    rows += [dict(accession='reference',gsm=None,title=Path(u).name,url=u,assembly=a) for a,u in references.items()]
    return rows
def fetch(row):
    p=BASE/'cache/direct_marks'/row['accession']/row['assembly']/row['url'].rsplit('/',1)[1] if row['accession']=='reference' else BASE/'cache/direct_marks'/row['accession']/row['url'].rsplit('/',1)[1]
    p.parent.mkdir(parents=True,exist_ok=True)
    cached=p.exists()
    if not p.exists():
        tmp=p.with_suffix(p.suffix+'.partial')
        for attempt in range(3):
            try:
                req=urllib.request.Request(row['url'],headers={'User-Agent':'scRNA_seq public research/1.0'})
                with urllib.request.urlopen(req,timeout=60) as response,tmp.open('wb') as h:
                    expected=int(response.headers.get('Content-Length',0))
                    if expected>256*1024**2:raise ValueError('File exceeds frozen 256 MiB ceiling')
                    total=0
                    while chunk:=response.read(1024**2):
                        total+=len(chunk)
                        if total>256*1024**2:raise ValueError('File exceeded byte ceiling')
                        h.write(chunk)
                    if expected and total!=expected:raise ValueError('Incomplete response')
                tmp.replace(p)
                break
            except urllib.error.URLError as exc:
                # Windows curl uses the system trust store; certificate validation
                # stays enabled. Fallback is restricted to the official references.
                if row['accession']=='reference' and isinstance(exc.reason,ssl.SSLCertVerificationError):
                    subprocess.run(['curl.exe','--fail','--location','--silent','--show-error',
                        '--max-time','180','--max-filesize',str(256*1024**2),
                        '--output',str(tmp),row['url']],check=True)
                    assert 0<tmp.stat().st_size<=256*1024**2
                    tmp.replace(p)
                    break
                if attempt==2:raise
                time.sleep(2)
            except Exception:
                if attempt==2:raise
                time.sleep(2)
    return {**row,'path':p.relative_to(BASE).as_posix(),'bytes':p.stat().st_size,'sha256':sha(p),
            'cache_reused':cached,'checked_utc':datetime.now(timezone.utc).isoformat()}
def main():
    ap=argparse.ArgumentParser(description=__doc__); ap.add_argument('--execute',action='store_true'); args=ap.parse_args()
    rows=candidates()
    report=BASE/'reports/direct_mark_input_inventory.json'
    # Execute is explicitly bounded to these processed tracks and annotations;
    # no raw reads or whole-series tar archives are downloaded.
    if not args.execute:
        print(json.dumps(rows,indent=2));return
    records={r['url']:r for r in json.loads(report.read_text()).get('files',[])} if report.exists() else {}
    def save():
        report.write_text(json.dumps({'scope':'Processed direct-mark tracks, methylation domains and assembly-matched annotations; no inferred biological replication',
              'max_file_bytes':256*1024**2,'maximum_planned_files':len(rows),'files':list(records.values())},indent=2)+'\n',encoding='utf-8')
    with ThreadPoolExecutor(max_workers=3) as pool:
        tasks={pool.submit(fetch,r):r for r in rows}
        for task in as_completed(tasks):
            row=tasks[task]
            try:
                record=task.result();records[row['url']]=record
                print(row['accession'],Path(record['path']).name,round(record['bytes']/1024**2,1),'MiB',flush=True)
            except Exception as exc:
                records[row['url']]={**row,'status':'unavailable','error':str(exc)}
                print(row['accession'],'unavailable',str(exc),flush=True)
            save()

if __name__=='__main__':main()
