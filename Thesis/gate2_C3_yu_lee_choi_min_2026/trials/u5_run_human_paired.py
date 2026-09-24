"""Continue the authorized human stages after every library is ready."""
from pathlib import Path
import sys,os,json,time,subprocess
from datetime import datetime,timezone
PAPER=Path(__file__).resolve().parents[1];ROOT=PAPER.parents[1];sys.path.insert(0,str(ROOT))
from analysis.lib.provenance import write_json_atomic,sha256_file,code_identity

def main():
 out=PAPER/'trials/u5_human_niche';out.mkdir(exist_ok=True);cache=PAPER/'cache/u5_human_niche';cache.mkdir(exist_ok=True)
 record=out/'paired_pipeline_run_record.json';assert not record.exists();scripts=['u5_human_aggregate.py','u5_human_paired_pathways.R','u5_human_compatibility.py','u5_human_ligand_refits.R']
 spec=dict(frozen_utc=datetime.now(timezone.utc).isoformat(),scripts={s:sha256_file(PAPER/'trials'/s) for s in scripts},stages=['all-library aggregation','paired-patient TMM/voom/CAMERA with global BH across eligible histology/view/set combinations per sensitivity family','paired descriptive RNA compatibility','conditional paired target refits'],primary=dict(uncertainty=.2,library_mode='pooled',cell_floor=50,correlation='estimated'),sensitivity='one-at-a-time confidence .3, largest-library choice, cell floor 30/100, fixed correlation .01; compatibility prior .5/2',no_results_based_changes=True)
 write_json_atomic(out/'paired_pipeline_specification.json',spec);state=dict(status='running',pid=os.getpid(),started_utc=spec['frozen_utc'],stage='await_all_library_inputs',code=code_identity(ROOT,__file__),completed=[]);write_json_atomic(record,state);t=time.monotonic()
 try:
  while True:
   upstream=json.loads((PAPER/'trials/u5_human_full/processing_run_record.json').read_text())
   if upstream['status']=='completed_inputs_ready':break
   if upstream['status']=='failed':raise RuntimeError('Human full input processing failed')
   time.sleep(5)
  for name in scripts:
   assert sha256_file(PAPER/'trials'/name)==spec['scripts'][name],name+' changed after launch'
   state['stage']=name;write_json_atomic(record,state)
   if name.endswith('.R'):cmd=[str(ROOT/'analysis/corrections/statistics/.tools/R-portable/app/bin/Rscript.exe'),str(PAPER/'trials'/name)]
   else:cmd=[sys.executable,str(ROOT/'analysis/scripts/run_with_environment.py'),'--site-packages',str(ROOT/'.venv-x64/Lib/site-packages'),str(PAPER/'trials'/name)]
   with (cache/(name+'.log')).open('w',encoding='utf-8') as log:subprocess.run(cmd,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,check=True)
   state['completed'].append(name);state['elapsed_seconds']=round(time.monotonic()-t,1);write_json_atomic(record,state);print('Completed human stage',name,flush=True)
  state.update(status='completed_computations_require_reports',stage='source_filter_and_reports')
 except Exception as e:state.update(status='failed',error=repr(e));raise
 finally:state.update(updated_utc=datetime.now(timezone.utc).isoformat(),elapsed_seconds=round(time.monotonic()-t,1));write_json_atomic(record,state)

if __name__=='__main__':main()
