"""Run the second-batch sensitivities without changing first-batch evidence."""
from pathlib import Path
from datetime import datetime,timezone
import argparse
import hashlib
import json
import subprocess

BASE=Path(__file__).resolve().parents[1]
ROOT=BASE.parents[1]
def sha(p):
    with p.open('rb') as h:return hashlib.file_digest(h,'sha256').hexdigest()

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--rscript',type=Path,required=True)
    args=ap.parse_args()
    out=BASE/'tables/ire1_stability_2026-09-25'
    record_path=BASE/'reports/ire1_stability_run.json'
    if out.exists() or record_path.exists():raise SystemExit('Refusing to overwrite the stability run')
    original=json.loads((BASE/'reports/ire1_run.json').read_text())
    assert original['status']=='completed'
    for name,digest in original['output_sha256'].items():assert sha(BASE/name)==digest,name
    inputs=[BASE/'config/second_batch.md',BASE/'processed/ire1/counts.tsv',BASE/'reports/ire1_run.json']
    inputs += [BASE/name for name in original['output_sha256']]
    script=Path(__file__).with_name('12_fit_ire1_stability.R')
    record={'status':'running','started_utc':datetime.now(timezone.utc).isoformat(),
            'scope':'Post-primary sensitivity analysis; no changed primary result or new confirmation',
            'input_sha256':{p.relative_to(BASE).as_posix():sha(p) for p in inputs},
            'code_sha256':{p.name:sha(p) for p in [Path(__file__),script]},
            'output_directory':out.relative_to(BASE).as_posix()}
    out.mkdir(parents=True)
    record_path.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    log=BASE/'reports/ire1_stability_R.log'
    with log.open('w',encoding='utf-8') as h:
        proc=subprocess.run([str(args.rscript.resolve()),'--vanilla',str(script),str(BASE),str(out)],cwd=ROOT,stdout=h,stderr=subprocess.STDOUT)
    record['exit_code']=proc.returncode
    record['status']='completed' if proc.returncode==0 else 'failed'
    record['finished_utc']=datetime.now(timezone.utc).isoformat()
    record['output_sha256']={p.relative_to(BASE).as_posix():sha(p) for p in sorted(out.iterdir()) if p.is_file()}
    record['log_sha256']=sha(log)
    for name,digest in original['output_sha256'].items():assert sha(BASE/name)==digest,name
    record['first_batch_outputs_unchanged']=True
    record_path.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    print(record['status'],record['finished_utc'],flush=True)
    raise SystemExit(proc.returncode)

if __name__=='__main__':main()
