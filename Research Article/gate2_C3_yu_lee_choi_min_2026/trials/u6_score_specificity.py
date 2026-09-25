"""New HPCS/repair specificity extension; existing trials remain unchanged."""
from pathlib import Path
import sys,os,json,time,importlib.util,itertools
from datetime import datetime,timezone
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ.setdefault(k,'2')
PAPER=Path(__file__).resolve().parents[1];ROOT=PAPER.parents[1];sys.path.insert(0,str(ROOT))
from analysis.lib.provenance import write_json_atomic,code_identity,sha256_file

def main():
    import numpy as np,pandas as pd,anndata as ad
    out=PAPER/'trials/u6_specificity';cache=PAPER/'cache/u6_specificity';cache.mkdir(parents=True,exist_ok=True)
    record=out/'run_record.json';assert not record.exists()
    config=json.loads((out/'module_specification.json').read_text());modules=config['modules']
    helper_path=ROOT/'Research Article/epithelial_state_specificity/run_es1.py';loader=importlib.util.spec_from_file_location('existing_es1_readonly',helper_path);helper=importlib.util.module_from_spec(loader);loader.loader.exec_module(helper)
    genes=sorted({g for m in modules for g in m['genes']}|set(helper.AIRWAY+helper.ALVEOLAR+helper.LABELS));pos={g:i for i,g in enumerate(genes)}
    state=dict(status='running',pid=os.getpid(),started_utc=datetime.now(timezone.utc).isoformat(),code=code_identity(ROOT,__file__),specification_sha256=sha256_file(out/'module_specification.json'),readonly_helper_sha256=sha256_file(helper_path),inputs=[]);write_json_atomic(record,state);start=time.monotonic()
    rows=[];gene_rows=[];retention=[]
    def score_named(counts,totals,present,group_labels,meta):
        keep=totals>=2000;retention.append(dict(**meta,n_input=len(totals),n_depth=int(keep.sum())))
        for seed in [17,29]:
            sampled=helper.joint_sample(counts[:,keep],totals[keep],2000,seed);groups=np.asarray(group_labels)[keep]
            for group in sorted(set(groups)):
                selected=groups==group;n=int(selected.sum());detection=np.mean(sampled[:,selected]>0,axis=1)
                for m in modules:
                    gs=[g for g in m['genes'] if g in present];fraction=len(gs)/len(m['genes']);value=float(np.mean(detection[[pos[g] for g in gs]])) if gs else np.nan
                    rows.append(dict(**meta,seed=seed,group=group,n_cells=n,module=m['name'],kind=m['kind'],n_genes=len(gs),total_genes=len(m['genes']),coverage=fraction,score=value,interpretable=bool(n>=30 and fraction>=.8),label_overlap=';'.join(m.get('label_overlap',[]))))
    try:
        overlap=[]
        for a,b in itertools.combinations(modules,2):
            x,y=set(a['genes']),set(b['genes']);overlap.append(dict(module_a=a['name'],module_b=b['name'],intersection=len(x&y),jaccard=len(x&y)/len(x|y),genes=';'.join(sorted(x&y))))
        pd.DataFrame(overlap).to_csv(out/'module_overlap.csv',index=False)
        descriptors=[('GSE247130','GSE247130_Aggregate.P9_filtered_feature_bc_matrix.h5','P9',['P9_control','P9_Cebpa_mutant'],['control','Cebpa_mutant'],'uninjured'),('GSE247130','GSE247130_Aggregate7-wk_filtered_feature_bc_matrix.h5','7wk',['7wk_control','7wk_Cebpa_mutant'],['control','Cebpa_mutant'],'uninjured'),('GSE247130','GSE247130_Aggregate_SeV_filtered_feature_bc_matrix.h5','adult_14dpi',['SeV_control','SeV_Cebpa_mutant'],['control','Cebpa_mutant'],'SeV'),('GSE310539','GSE310539_totalaggr_filtered_feature_bc_matrix.h5','adult',['wildtype_PBS','wildtype_SeV','AP1mut_PBS','AP1mut_SeV'],['wildtype','wildtype','AP1mut','AP1mut'],None)]
        for cohort,name,stage,units,genotypes,injury in descriptors:
            path=ROOT/'raw_data'/cohort/name;state['inputs'].append(dict(path=str(path.relative_to(ROOT)),sha256=sha256_file(path)))
            counts,totals,suffix,present=helper.read_10x(path,genes)
            for i,(unit,genotype) in enumerate(zip(units,genotypes),1):
                selected=suffix==str(i);meta=dict(cohort=cohort,unit=unit,unit_type='pooled_library',genotype=genotype,stage=stage,injury=injury or ('SeV' if unit.endswith('_SeV') else 'uninjured'),independence='pooled well, not replicated animals; new HPCS scoring, original ES1 files unchanged')
                helper.score_unit(counts[:,selected],totals[selected],genes,present,meta,modules,rows,gene_rows,retention)
            del counts;print(cohort,stage,'new specificity scores done',flush=True)
        # Deposited HPCS counts and author labels; aggregate sorting populations
        # into their shared identifier, never treat six sort labels as six mice.
        path=PAPER/'cache/extension_pilot/GSM8529524_IGO14143-slc4a11-traced.h5ad';state['inputs'].append(dict(path=str(path.relative_to(ROOT)),sha256=sha256_file(path)))
        a=ad.read_h5ad(path);x=a.layers['counts'].tocsr();assert np.all(x.data>=0) and np.all(x.data==np.round(x.data));assert a.var_names.is_unique
        counts=np.zeros((len(genes),a.n_obs),dtype=np.int32);present=set(a.var_names)
        for j,g in enumerate(genes):
            if g in present:counts[j]=x[:,a.var_names.get_loc(g)].toarray().ravel().astype(np.int32)
        totals=np.asarray(x.sum(axis=1)).ravel().astype(np.int64);identifiers=a.obs.Classification.astype(str).str.extract(r'^(AV\d+)',expand=False);assert identifiers.notna().all()
        a.obs.assign(source_identifier=identifiers).groupby(['source_identifier','Classification','cell type'],observed=True).size().rename('cells').to_csv(out/'HPCS_source_identifier_coverage.csv')
        for identifier in sorted(identifiers.unique()):
            take=identifiers.eq(identifier).to_numpy();meta=dict(cohort='GSE277777',unit=identifier,unit_type='source_identifier_at_most_two_animals',genotype='KP_Slc4a11_MCD_Hipp11_reporter',stage='16wk_2wk_trace',injury='neoplastic_sorted',independence='multiple sorting gates share source identifier; author state/signature not independent validation')
            score_named(counts[:,take],totals[take],present,a.obs.loc[take,'cell type'].astype(str).to_numpy(),meta)
        del a,x,counts
        # Independently annotated early NNK alveolar compartments.
        ann=pd.read_csv(PAPER/'cache/u3_lineage_annotation/cell_annotations.csv.gz',index_col=0);ann=ann[ann.compartment=='alveolar']
        for gsm,g in ann.groupby('gsm'):
            path=PAPER/'cache/u3_qc'/(gsm+'_qc.h5ad');state['inputs'].append(dict(path=str(path.relative_to(ROOT)),sha256=sha256_file(path)));a=ad.read_h5ad(path);a=a[g.index].copy();symbols=a.var.gene_symbol.astype(str).to_numpy();present=set(symbols)
            counts=np.zeros((len(genes),a.n_obs),dtype=np.int32)
            for j,gene in enumerate(genes):
                ix=np.flatnonzero(symbols==gene)
                if len(ix):counts[j]=np.asarray(a.X[:,ix].sum(axis=1)).ravel().astype(np.int32)
            totals=np.asarray(a.X.sum(axis=1)).ravel().astype(np.int64);meta=dict(cohort='GSE300288',unit=gsm,unit_type='retained_mouse_library',genotype='NNK_model',stage='3mo_post_NNK',injury=g.treatment.iloc[0],independence='source-verified mouse library; independent broad annotation, KAC not assigned')
            score_named(counts,totals,present,g.loc[a.obs_names,'label'].to_numpy(),meta);del a,counts
        result=pd.DataFrame(rows);result.to_csv(out/'module_scores.csv',index=False);pd.DataFrame(retention).to_csv(out/'retention.csv',index=False)
        pd.DataFrame(gene_rows).to_csv(cache/'pooled_well_gene_detection.csv.gz',index=False,compression='gzip')
        assert result.score.dropna().between(0,1).all() and result.coverage.between(0,1).all()
        assert not result.duplicated(['cohort','unit','seed','group','module']).any()
        state.update(status='completed_scores_require_report',rows=len(result),source_identifiers_HPCS=2,technical_seeds=[17,29],independence='no cell-level P values; no pseudo-replication of sorting gates or pooled wells')
    except Exception as exc:state.update(status='failed',error=type(exc).__name__+': '+str(exc));raise
    finally:state.update(updated_utc=datetime.now(timezone.utc).isoformat(),elapsed_seconds=round(time.monotonic()-start,1));write_json_atomic(record,state)

if __name__=='__main__':main()
