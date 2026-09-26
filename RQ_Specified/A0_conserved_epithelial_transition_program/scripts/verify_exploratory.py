"""Verify executed estimates, frozen definitions and provenance; close out A0."""
import ast
import json
import re
import sys
from datetime import datetime, timezone
from urllib.parse import unquote
import numpy as np
import pandas as pd
import scipy
from scipy import sparse
from exploratory_common import BASE, ROOT, PROC, TABLES, sha, save_json, rank_scores, paired_effects, learn_program, excluded, resolve_source_path


def assert_effects(observed,expected):
    keys=['unit','endpoint']
    left=observed.sort_values(keys).reset_index(drop=True)
    right=expected.sort_values(keys).reset_index(drop=True)
    assert left[keys].equals(right[keys])
    for c in ['n_start','n_intermediate','n_end','all_states_at_least_30']:
        assert left[c].equals(right[c]),c
    for c in ['difference','intermediate_mean','endpoint_mean','standardized_difference']:
        np.testing.assert_allclose(left[c],right[c],rtol=1e-10,atol=1e-12,equal_nan=True)


def main():
    now=datetime.now(timezone.utc).isoformat()
    save_json(BASE/'exploratory_validation.json',{'status':'IN_PROGRESS','started_at_utc':now})
    save_json(BASE/'exploratory_execution_record.json',{'status':'IN_PROGRESS','started_at_utc':now})
    config=json.loads((BASE/'exploratory_config.json').read_text())
    program=json.loads((BASE/'frozen_repair_program.json').read_text())
    variant=json.loads((BASE/'frozen_control_variant.json').read_text())
    manifest=json.loads((BASE/'developmental_expression_manifest.json').read_text())
    assert program['config_sha256']==sha(BASE/'exploratory_config.json')
    assert manifest['frozen_program_sha256']==sha(BASE/'frozen_repair_program.json')
    assert len(program['genes'])==50 and len(set(program['genes']))==50
    assert len(config['reference_genes'])==len(set(config['reference_genes']))==1000
    assert not any(excluded(g,config) for g in program['genes'])
    generic=set().union(*(set(config['control_modules'][n]) for n in config['generic_control_names']))
    assert variant['genes']==[g for g in program['genes'] if g not in generic] and len(variant['genes'])==31
    for name,digest in config['source_hashes'].items():assert sha(resolve_source_path(name))==digest,name
    records=[]
    for provenance in sorted((BASE/'cache/expression').glob('*.provenance.json')):
        record=json.loads(provenance.read_text());path=provenance.with_name(provenance.name.removesuffix('.provenance.json'))
        assert path.stat().st_size==record['bytes'] and sha(path)==record['sha256']
        records.append({'file':path.relative_to(BASE).as_posix(),**record})
    for record in manifest['chunks']:
        path=BASE/record['file'];assert path.stat().st_size==record['bytes'] and sha(path)==record['sha256']
        records.append(record)
    assert len(records)==74
    save_json(BASE/'expression_source_manifest.json',{'verified_at_utc':now,'files':records,
        'expression_files_and_api_chunks':len(records),'developmental_input':'SCT normalized data-slot export',
        'additional_source_hashes':config['source_hashes'],'metadata_source_manifest':'source_manifest.json'})

    # Fixtures use explicit expected values, including all-zero ties and overlap.
    x=np.array([[0.,0,0,0],[3,1,0,2],[0,4,1,3]])
    genes=['q1','q2','b1','b2'];q=['q1','q2'];bg=['b1','b2']
    score,_=rank_scores(x,genes,q,bg)
    np.testing.assert_allclose(score,[.5,.75,.5])
    np.testing.assert_allclose(rank_scores(x*np.array([2,7,.3])[:,None],genes,q,bg)[0],score)
    np.testing.assert_allclose(rank_scores(sparse.csr_matrix(x),genes,q,['q1']+bg)[0],score)
    _,coverage=rank_scores(x,genes,q+['absent'],bg)
    assert coverage['present_genes']==2 and coverage['feature_coverage']==2/3

    facts=pd.read_csv(TABLES/'repair_pseudobulk_inventory.csv')
    pb=np.load(PROC/'repair_pseudobulk.npz')
    repair_genes=json.loads((PROC/'repair_genes.json').read_text())
    selected,_=learn_program(facts,pb['logcpm'],pb['detection'],repair_genes,config)
    assert selected==program['genes']
    assert len(facts)==27 and (facts.cells>=30).all() and facts.unit.nunique()==9
    heldout=config['repair_units'][0]
    selected,_=learn_program(facts,pb['logcpm'],pb['detection'],repair_genes,config,omit=heldout)
    folds=pd.read_csv(TABLES/'repair_fold_programs.csv')
    assert selected==folds.loc[folds.heldout_unit==heldout,'gene'].tolist()
    perturbed=pb['logcpm'].copy();perturbed[facts.unit.eq(heldout)]=10000
    changed_detection=pb['detection'].copy();changed_detection[facts.unit.eq(heldout)]=0
    other,_=learn_program(facts,perturbed,changed_detection,repair_genes,config,omit=heldout)
    assert other==selected,'Held-out data leaked into discovery'
    cv=pd.read_csv(TABLES/'repair_heldout_effects.csv')
    e1=json.loads((BASE/'stage_E1_result.json').read_text())
    for endpoint,d in cv.groupby('endpoint'):
        assert len(d)==9 and (d.difference>0).sum()==e1['heldout_results'][endpoint]['positive_mice']
        np.testing.assert_allclose(d.difference.median(),e1['heldout_results'][endpoint]['median_difference'])
    repair_obs=pd.read_csv(PROC/'repair_scored_obs.csv')
    repair_scores=np.load(PROC/'repair_cell_scores.npz')
    assert_effects(paired_effects(repair_obs,repair_scores['candidate']),pd.read_csv(TABLES/'repair_apparent_effects.csv'))
    depth=pd.read_csv(PROC/'depth_matched_cells.csv')
    assert not depth.cell_barcode.duplicated().any()
    counts=depth.groupby(['unit','state']).size()
    assert (counts>=30).all() and (counts.groupby('unit').nunique()==1).all()
    matched_scores=repair_scores['candidate'][repair_obs.set_index('cell_barcode').index.get_indexer(depth.cell_barcode)]
    d=depth.merge(repair_obs[['cell_barcode','time']],on='cell_barcode',validate='one_to_one')
    assert_effects(paired_effects(d,matched_scores),pd.read_csv(TABLES/'repair_depth_matched_effects.csv'))

    transfer_evidence=[]
    for dataset,states,metadata_file,index_name,qualified in [
        ('intestine',config['transfer']['intestinal_states'],'v1_haber_verified_mouse_coverage.csv','mouse_id',2),
        ('development',config['transfer']['developmental_states'],'d2_negretti_encoded_group_coverage.csv','encoded_library_group',1)]:
        obs=pd.read_csv(PROC/f'{dataset}_scored_obs.csv')
        assert obs.cell_id.is_unique
        scores=np.load(PROC/f'{dataset}_transfer_scores.npz')
        effects=pd.read_csv(TABLES/f'{dataset}_transfer_effects.csv')
        coverage=pd.read_csv(TABLES/f'{dataset}_transfer_coverage.csv')
        prior=pd.read_csv(BASE/'tables'/metadata_file).set_index(index_name)
        assert (coverage[coverage.module=='candidate'].feature_coverage>=.8).all()
        assert (coverage.reference_coverage==1).all()
        for module in effects.module.unique():
            estimate=paired_effects(obs,scores[module],start=states[0],intermediate=states[1],end=states[2],dataset=dataset)
            assert_effects(estimate,effects[effects.module==module])
        for unit,frame in obs[obs.state.isin(states)].groupby('unit'):
            for state in states:assert (frame.state==state).sum()==prior.loc[unit,state]
        primary=effects[effects.module=='candidate']
        assert primary[primary.all_states_at_least_30].unit.nunique()==qualified
        # Independently compute direct pairwise ranks for representative stored cells.
        matrix=np.load(PROC/f'{dataset}_scoring_matrix.npz')['matrix']
        g=json.loads((PROC/f'{dataset}_scoring_genes.json').read_text());lookup={v:i for i,v in enumerate(g)}
        q=[lookup[v] for v in program['genes'] if v in lookup]
        bg=[lookup[v] for v in config['reference_genes'] if v in lookup and v not in program['genes']]
        indexes=np.linspace(0,len(obs)-1,12,dtype=int)
        for i in indexes:
            a=matrix[i,q][:,None];b=matrix[i,bg][None,:]
            expected=((a>b)+.5*(a==b)).mean()
            np.testing.assert_allclose(expected,scores['candidate'][i],rtol=1e-12,atol=1e-12)
        transfer_evidence.append({'dataset':dataset,'qualified_units':qualified,'independent_score_reference_cells':len(indexes),
            'effect_rows_reconciled':len(effects),'cell_counts_match_pre_expression_inventory':True})

    draws=pd.read_csv(PROC/'intestine_depth_matching_draws.csv')
    summary=pd.read_csv(TABLES/'intestine_depth_matched_sensitivity.csv')
    for _,r in summary.iterrows():
        d=draws[(draws.unit==r.unit)&(draws.module==r.module)&(draws.endpoint==r.endpoint)]
        assert len(d)==200 and d.iteration.nunique()==200
        np.testing.assert_allclose([d.difference.median(),d.difference.quantile(.05),d.difference.quantile(.95)],
            [r.median_difference,r.draw_q05,r.draw_q95],atol=1e-12)
    figures=json.loads((BASE/'exploratory_figure_review.json').read_text())
    assert len(figures['figures'])==3
    for r in figures['figures']:assert r['status']=='visually_checked' and sha(BASE/r['path'])==r['sha256']
    files=[p for p in BASE.rglob('*') if p.is_file() and not any(x in p.relative_to(BASE).parts for x in ['cache','processed','__pycache__'])]
    syntax=0;links=0
    for path in files:
        if path.suffix=='.py':ast.parse(path.read_text(encoding='utf-8'));syntax+=1
        if path.suffix=='.json':json.loads(path.read_text(encoding='utf-8'))
        if path.suffix=='.md':
            for target in re.findall(r'\[[^\]]*\]\(([^)]+)\)',path.read_text(encoding='utf-8')):
                if target.startswith(('http:','https:','mailto:','#')):continue
                assert (path.parent/unquote(target.split('#',1)[0].strip('<>'))).exists(),(path.name,target)
                links+=1
    validation={'status':'PASS','completed_at_utc':datetime.now(timezone.utc).isoformat(),
        'scope':'amended exploratory pilot; separate from pilot_v1; not independent replication','frozen_config_and_program_hashes_match':True,
        'expression_sources_verified':len(records),'additional_config_source_hashes_verified':len(config['source_hashes']),
        'rank_tie_scaling_sparse_overlap_and_missing_feature_fixtures':'PASS','full_gene_selection_reproduced':True,
        'representative_heldout_fold_reproduced_and_perturbation_invariant':True,'saved_repair_and_transfer_effects_reconciled':True,
        'transfer_checks':transfer_evidence,'matching_draw_summaries_reconciled':True,'figures_visually_checked':3,
        'python_scripts_parsed':syntax,'local_markdown_links_checked':links,
        'ingestion_checks_reused':'Original ingestion checked matrix/annotation alignment, all intestinal raw integer counts, all repair nonzeros and every developmental chunk feature/row mapping; unchanged input hashes verified here.'}
    save_json(BASE/'exploratory_validation.json',validation)
    decisions=json.loads((BASE/'stage_decisions.json').read_text())
    decisions['decisions']=[d for d in decisions['decisions'] if d['stage']!='E4']
    decisions['decisions'].append({'stage':'E4','reviewed_at_utc':validation['completed_at_utc'],'decision':'bounded_exploratory_pilot_complete',
        'reason':'All warranted stages executed; numerical, provenance, freeze, link and figure checks pass. Broad expansion stopped at E3a; original replication gaps remain.',
        'evidence':'exploratory_validation.json'})
    save_json(BASE/'stage_decisions.json',decisions)
    current={'current_status':'EXPLORATORY_PILOT_COMPLETE_NARROW','biological_decision':'NARROW_REPAIR_ASSOCIATED_PROGRAM; BROAD_TRANSFER_UNSUPPORTED',
        'exploratory_pilot_complete':True,'original_three_context_pilot_complete':False,'full_pilot_complete':False,
        'primary_analysis_ready':False,'legacy_field_scope':'full_pilot_complete and primary_analysis_ready refer to the original PLAN.md design',
        'expression_program_learned':True,'transfer_test_performed':True,'completed':['E1 repair discovery and mouse-held-out checks','E2 source specificity and technical challenges',
            'E3 frozen descriptive developmental and intestinal transfer','E3a targeted post-transfer intestinal depth sensitivity','E4 report, figures and validation'],
        'not_established':['universal process','causal transition regulator','robust cross-tissue conservation'],
        'remaining_original_design_requirements':'feasibility_readiness.json','current_report':'reports/EXPLORATORY_PILOT_REPORT.md',
        'required_work_remaining_in_authorized_exploratory_scope':[],'future_investment':'Do not broaden this signature analysis now. A separate bounded independent normal-transition validation could revisit weaker shared components.'}
    save_json(BASE/'exploratory_readiness.json',current)
    historic=json.loads((BASE/'exploratory_decisions.json').read_text());historic['current_stage']='E4_exploratory_complete_narrow'
    historic['exploratory_pilot_complete']=True;save_json(BASE/'exploratory_decisions.json',historic)
    files=[p for p in BASE.rglob('*') if p.is_file() and not any(x in p.relative_to(BASE).parts for x in ['cache','processed','__pycache__'])]
    artifacts={p.relative_to(BASE).as_posix():{'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(files) if p.name!='exploratory_execution_record.json'}
    execution={'completed_at_utc':validation['completed_at_utc'],'scope':'amended exploratory E1-E4; snapshot includes historical context artifacts',
        'status':current['current_status'],'validation':'PASS','python_executable':sys.executable,'python_version':sys.version.split()[0],
        'numpy_version':np.__version__,'pandas_version':pd.__version__,'scipy_version':scipy.__version__,
        'numerical_packages_source':'.venv-x64/Lib/site-packages','frozen_config_sha256':sha(BASE/'exploratory_config.json'),
        'frozen_program_sha256':sha(BASE/'frozen_repair_program.json'),'artifacts':artifacts}
    save_json(BASE/'exploratory_execution_record.json',execution)
    print(json.dumps({'status':'PASS','sources':len(records),'artifacts':len(artifacts),'figures':3,'decision':'NARROW','required_work_remaining':0}))


if __name__=='__main__':main()
