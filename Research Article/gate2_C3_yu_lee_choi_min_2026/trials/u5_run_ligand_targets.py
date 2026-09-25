"""Freeze and launch donor-level target refits; expression filters follow separately."""
from pathlib import Path
import json,sys,os,subprocess,time
from datetime import datetime,timezone
PAPER=Path(__file__).resolve().parents[1];ROOT=PAPER.parents[1];sys.path.insert(0,str(ROOT))
from analysis.lib.provenance import write_json_atomic,code_identity,sha256_file

def main():
    out=PAPER/'trials/u5_ligand_targets';cohort=sys.argv[1]
    here=out/cohort;here.mkdir(exist_ok=True)
    assert json.loads((out/'prior_acquisition.json').read_text())['status']=='completed'
    script=PAPER/'trials/u5_ligand_target_refits.R'
    spec=dict(frozen_utc=datetime.now(timezone.utc).isoformat(),cohort=cohort,scope='extension after observational DE inspection; descriptive ligand prioritization, no independent validation',metric='NicheNet-prior Pearson statistic, base-R implementation; not full nichenetr package',prior_record_sha256=sha256_file(out/'prior_acquisition.json'),receiver_targets='donor-level edgeR TMM / limma voom DE q<.05; up and down separate; minimum ten mapped genes',background='all tested receiver genes intersect prior rows; mapping coverage reported',loo='repeat expression filtering, TMM, voom, DE, BH and target selection after excluding each donor; minimum three units per arm retained',cell_floor=50,duplicate_AT2_broad_view='omitted because identical to sole subtype',candidate_ranking='all prior ligands scored here; source-filtered planned ligand family rankings produced separately; unmeasured source/receptor genes remain unevaluable',refit_script_sha256=sha256_file(script))
    write_json_atomic(here/'specification.json',spec)
    state=dict(status='running',pid=os.getpid(),started_utc=spec['frozen_utc'],code=code_identity(ROOT,__file__));start=time.monotonic();write_json_atomic(here/'run_record.json',state)
    try:
        subprocess.run([str(ROOT/'analysis/corrections/statistics/.tools/R-portable/app/bin/Rscript.exe'),str(script),cohort],cwd=ROOT,check=True)
        state['status']='completed_refits_candidate_filtering_pending'
    except Exception as exc:state.update(status='failed',error=str(exc));raise
    finally:state.update(updated_utc=datetime.now(timezone.utc).isoformat(),elapsed_seconds=round(time.monotonic()-start,1));write_json_atomic(here/'run_record.json',state)

if __name__=='__main__':main()
