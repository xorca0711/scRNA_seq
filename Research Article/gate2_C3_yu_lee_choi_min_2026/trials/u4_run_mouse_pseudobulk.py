from pathlib import Path
import json,sys,subprocess,time,os
from datetime import datetime,timezone
PAPER=Path(__file__).resolve().parents[1];ROOT=PAPER.parents[1];sys.path.insert(0,str(ROOT))
from analysis.lib.provenance import code_identity,sha256_file,write_json_atomic

def main():
    out=PAPER/'trials/u4_mouse_niche';out.mkdir(exist_ok=True)
    record=out/'pseudobulk_run_record.json';assert not record.exists()
    script=PAPER/'trials/u4_mouse_pseudobulk.R'
    state=dict(status='running',pid=os.getpid(),started_utc=datetime.now(timezone.utc).isoformat(),code=code_identity(ROOT,__file__),R_script_sha256=sha256_file(script),resource_spec_sha256=sha256_file(PAPER/'trials/u4_resources/specification.json'));write_json_atomic(record,state);start=time.monotonic()
    try:
        subprocess.run([str(ROOT/'analysis/corrections/statistics/.tools/R-portable/app/bin/Rscript.exe'),str(script)],cwd=ROOT,check=True)
        state['status']='completed_requires_report_and_checks'
    except Exception as exc:state.update(status='failed',error=str(exc));raise
    finally:state.update(elapsed_seconds=round(time.monotonic()-start,1),updated_utc=datetime.now(timezone.utc).isoformat());write_json_atomic(record,state)

if __name__=='__main__':main()
