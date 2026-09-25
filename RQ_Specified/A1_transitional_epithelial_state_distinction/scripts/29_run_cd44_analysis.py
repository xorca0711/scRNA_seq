"""Run the newly eligible paired CD44 contrasts after reference/identity checks."""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import json
import subprocess
import time
import urllib.parse
import urllib.request
import numpy as np
import pandas as pd
from a1_robustness import sha

BASE = Path(__file__).resolve().parents[1]
ROOT = BASE.parents[1]


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--rscript',type=Path,required=True)
    args=parser.parse_args()
    out=BASE/'tables/cd44_closure';rp=BASE/'reports/cd44_closure_run.json'
    if out.exists() or rp.exists():raise SystemExit('Refusing to overwrite CD44 run')
    cp=BASE/'config/closure_analysis_contract.json';cfg=json.loads(cp.read_text())['cd44']
    identity=json.loads((BASE/'reports/closure_identity_run.json').read_text())
    for name,digest in identity['output_sha256'].items():assert sha(BASE/name)==digest,name
    count=BASE/'cache/inputs/GSE273123/GSE273123_Count_Matrix.txt.gz'
    assert sha(count)==identity['input_sha256'][count.relative_to(BASE).as_posix()]
    meta_path=BASE/'tables/evidence_closure/cd44_sample_manifest.tsv'
    m=pd.read_csv(meta_path,sep='\t');x=pd.read_csv(count,sep='\t',index_col=0)
    x.index=x.index.astype(str)
    assert list(x.columns)==list(m.sample_id) and x.index.is_unique
    assert np.isfinite(x.to_numpy()).all() and (x.to_numpy()>=0).all() and (x.to_numpy()==np.floor(x.to_numpy())).all()
    assert (x.sum()>=100000).all()
    # Identifier lookup only; no expression-based choice or result inspection.
    cache=BASE/'cache/evidence_closure/gene_lookup';cache.mkdir(exist_ok=True)
    term='('+' OR '.join(s+'[sym]' for s in cfg['focus_symbols'])+') AND 10090[taxid]'
    url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&retmode=json&retmax=100&term='+urllib.parse.quote(term)
    urls=[]
    def get(name,url):
        p=cache/name
        if not p.exists():
            with urllib.request.urlopen(url,timeout=60) as r:data=r.read(20*1024**2+1)
            assert len(data)<=20*1024**2
            p.write_bytes(data)
        urls.append(dict(path=p.relative_to(BASE).as_posix(),url=url,sha256=sha(p)))
        return json.loads(p.read_text())
    ids=get('search.json',url)['esearchresult']['idlist'];time.sleep(.4)
    result=get('summary.json','https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&retmode=json&id='+','.join(ids))['result']
    mapped=[]
    for symbol in cfg['focus_symbols']:
        hits=[(i,result[i]) for i in result['uids'] if result[i]['name']==symbol and result[i]['organism']['taxid']==10090]
        assert len(hits)==1,(symbol,hits)
        ident,record=hits[0];assert ident in x.index
        mapped.append(dict(symbol=symbol,feature_id=ident,NCBI_name=record['name']))
    out.mkdir()
    pd.DataFrame(mapped).to_csv(out/'marker_mapping.tsv',sep='\t',index=False)
    x.to_csv(out/'counts_input.tsv.gz',sep='\t',index_label='feature_id',compression={'method':'gzip','mtime':0})
    # This compressed processed copy is small but ignored by Git; the deposited hash remains authoritative.
    m.to_csv(out/'sample_manifest.tsv',sep='\t',index=False)
    pd.DataFrame([dict(key='counts_md5',value=__import__('hashlib').md5((out/'counts_input.tsv.gz').read_bytes()).hexdigest()),
                  dict(key='samples_md5',value=__import__('hashlib').md5((out/'sample_manifest.tsv').read_bytes()).hexdigest())]).to_csv(out/'input_hashes.tsv',sep='\t',index=False)
    source=BASE/'scripts/30_fit_cd44.R'
    started=datetime.now(timezone.utc).isoformat()
    subprocess.run([str(args.rscript),str(source),str(out)],cwd=ROOT,check=True)
    record=dict(status='completed',started_utc=started,finished_utc=datetime.now(timezone.utc).isoformat(),
                contract_sha256=sha(cp),reference_map_sha256=sha(BASE/'reports/ANALYSIS_REFERENCE_MAP.md'),
                input_sha256={str(p.relative_to(BASE)):sha(p) for p in [count,meta_path]},gene_lookup=urls,
                code_sha256={p.name:sha(p) for p in [Path(__file__),source]},
                output_sha256={p.relative_to(BASE).as_posix():sha(p) for p in out.iterdir()})
    rp.write_text(json.dumps(record,indent=2)+'\n')
    print(pd.read_csv(out/'fit_summary.tsv',sep='\t').to_string(index=False))


if __name__=='__main__':main()
