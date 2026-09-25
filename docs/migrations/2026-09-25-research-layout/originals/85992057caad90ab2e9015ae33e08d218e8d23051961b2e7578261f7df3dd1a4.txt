"""Launch the two declared IPF pathway fits sequentially with provenance."""
from pathlib import Path
import json
import subprocess
import sys
import os
import time
from datetime import datetime,timezone
PAPER=Path(__file__).resolve().parents[1]
ROOT=PAPER.parents[1]
sys.path.insert(0,str(ROOT))
from analysis.lib.provenance import code_identity,sha256_file,write_json_atomic,archive_existing_record


def main():
    contract=json.loads((PAPER/'analysis_contract.json').read_text())
    if not contract['full_run_authorized']: raise RuntimeError('Execution not authorized')
    out=PAPER/'trials/u5_ipf_pathways';out.mkdir(exist_ok=True)
    (out/'.gitignore').write_text('*_DE.csv.gz\n*.log\n')
    spec_path=PAPER/'trials/u5_ipf_spec/specification.json'
    spec=json.loads(spec_path.read_text())
    for item in spec['input_sources']:
        if sha256_file(ROOT/item['path'])!=item['sha256']: raise ValueError('Frozen input changed: '+item['path'])
    record=out/'run_record.json';archive_existing_record(record)
    rscript=ROOT/'analysis/corrections/statistics/.tools/R-portable/app/bin/Rscript.exe'
    analysis=PAPER/'trials/u5_ipf_pathways.R'
    state={'started_utc':datetime.now(timezone.utc).isoformat(),'status':'running','pid':os.getpid(),
           'spec_sha256':sha256_file(spec_path),'R_code_sha256':sha256_file(analysis),
           'code':code_identity(ROOT,__file__),'completed_cohorts':[],
           'interpretation':'previously inspected observational IPF cohorts; not mouse treatment or independent validation'}
    start=time.monotonic();write_json_atomic(record,state)
    try:
        for cohort in spec['cohorts']:
            state['current_cohort']=cohort;write_json_atomic(record,state)
            with (out/(cohort+'.log')).open('w') as log:
                result=subprocess.run([str(rscript),'--vanilla',str(analysis),cohort],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
            if result.returncode: raise RuntimeError(f'{cohort} failed; see its local log')
            state['completed_cohorts'].append(cohort)
            state['elapsed_seconds']=round(time.monotonic()-start,2)
            write_json_atomic(record,state)
            print(cohort+' pathway fits completed',flush=True)
        state['status']='completed_results_require_interpretation'
        state['outputs']=[{'path':x.relative_to(ROOT).as_posix(),'sha256':sha256_file(x)} for x in sorted(out.glob('*')) if x.is_file() and x.name!=record.name and not x.name.endswith('.log')]
    except Exception as exc:
        state.update(status='failed',error=str(exc));raise
    finally:
        state.update(updated_utc=datetime.now(timezone.utc).isoformat(),elapsed_seconds=round(time.monotonic()-start,2))
        write_json_atomic(record,state)


if __name__=='__main__': main()
