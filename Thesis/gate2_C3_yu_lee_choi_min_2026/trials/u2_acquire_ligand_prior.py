"""Acquire and hash pinned public NicheNet v2 human prior; never execute downloads."""
from pathlib import Path
import hashlib,json,sys,urllib.request,time
from datetime import datetime,timezone
PAPER=Path(__file__).resolve().parents[1];ROOT=PAPER.parents[1];sys.path.insert(0,str(ROOT))
from analysis.lib.provenance import write_json_atomic,code_identity

def main():
    out=PAPER/'trials/u5_ligand_targets';out.mkdir(exist_ok=True)
    cache=PAPER/'cache/u5_ligand_targets';cache.mkdir(parents=True,exist_ok=True)
    files={'ligand_target_matrix_nsga2r_final.rds':'b09606b04b2d4490418d9028c0e58b9f','lr_network_human_21122021.rds':'2a155f81e9ffffd5d5e709fe66bcc465'}
    state=dict(status='running',started_utc=datetime.now(timezone.utc).isoformat(),code=code_identity(ROOT,__file__),files=[])
    write_json_atomic(out/'prior_acquisition.json',state);start=time.monotonic()
    try:
        for name,expected in files.items():
            url='https://zenodo.org/records/7074291/files/'+name+'?download=1';dest=cache/name
            if not dest.exists():
                temp=dest.with_suffix('.part');n=0
                with urllib.request.urlopen(url,timeout=90) as response,temp.open('wb') as h:
                    while chunk:=response.read(1048576):
                        n+=len(chunk)
                        if n>300_000_000:raise ValueError('300MB file cap exceeded')
                        h.write(chunk)
                temp.replace(dest)
            md5=hashlib.md5();sha=hashlib.sha256()
            with dest.open('rb') as h:
                while chunk:=h.read(1048576):md5.update(chunk);sha.update(chunk)
            assert md5.hexdigest()==expected,(name,'published MD5 mismatch')
            state['files'].append(dict(name=name,url=url,bytes=dest.stat().st_size,md5=md5.hexdigest(),sha256=sha.hexdigest()))
            write_json_atomic(out/'prior_acquisition.json',state);print(name,'verified',flush=True)
        state['status']='completed'
    except Exception as exc:state.update(status='failed',error=str(exc));raise
    finally:state.update(elapsed_seconds=round(time.monotonic()-start,1));write_json_atomic(out/'prior_acquisition.json',state)

if __name__=='__main__':main()
