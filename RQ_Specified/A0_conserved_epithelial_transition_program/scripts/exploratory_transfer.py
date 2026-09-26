"""E3: apply frozen scores to source-defined developmental and intestinal states."""
import csv
import gzip
import json
import struct
import numpy as np
import pandas as pd
from exploratory_common import BASE, PROC, TABLES, save_json, sha, rank_scores, paired_effects
from audit_extended_coverage import decode_annotations


def decode_numeric_packet(payload):
    u=lambda p:struct.unpack_from('<I',payload,p)[0]
    h=lambda p:struct.unpack_from('<H',payload,p)[0]
    def field(t,i):
        v=t-struct.unpack_from('<i',payload,t)[0]
        offset=h(v+4+2*i) if 4+2*i<h(v) else 0
        return t+offset if offset else None
    ptr=lambda p:p+u(p)
    def array(t,kind):
        v=ptr(field(t,0));n=u(v)
        if kind==5:return json.loads(payload[v+4:v+4+n].decode('utf-8'))
        dtype={1:'<f4',2:'<i4',3:'<u4',4:'<f8'}[kind]
        return np.frombuffer(payload,dtype=dtype,count=n,offset=v+4)
    root=u(0);nrows=u(field(root,0));ncols=u(field(root,1))
    vector=ptr(field(root,2));assert u(vector)==ncols
    columns=[]
    for i in range(ncols):
        col=ptr(vector+4+4*i);columns.append(array(ptr(field(col,1)),payload[field(col,0)]))
    assert all(len(x)==nrows for x in columns)
    ci=array(ptr(field(root,4)),payload[field(root,3)]) if field(root,3) else None
    ri=array(ptr(field(root,6)),payload[field(root,5)]) if field(root,5) else np.arange(nrows)
    return np.column_stack(columns).astype(np.float32),ci,ri


def get_modules(config):
    primary=json.loads((BASE/'frozen_repair_program.json').read_text())['genes']
    variant=json.loads((BASE/'frozen_control_variant.json').read_text())['genes']
    return {'candidate':primary,'candidate_without_generic_controls':variant,
        **{n:config['control_modules'][n] for n in config['generic_control_names']},
        'ADI_published_holdout':config['control_modules']['ADI_published_holdout']}


def load_intestine(config,modules):
    needed=set(config['reference_genes']).union(*(set(v) for v in modules.values()))
    mapping=config['transfer']['intestinal_batches_to_mouse']
    with gzip.open(BASE/'cache/expression/GSE92332_atlas_UMIcounts.txt.gz','rt') as f:
        header=next(csv.reader([next(f)],delimiter='\t'))
        ids=[h for h in header if h]
        assert len(ids)==7216
        parts=[h.split('_',2) for h in ids]
        obs=pd.DataFrame(parts,columns=['batch','barcode','state'])
        obs['cell_id']=ids;obs['unit']=obs.batch.map(mapping);obs['time']='homeostasis'
        indexes=np.flatnonzero(obs.unit.notna().to_numpy())
        obs=obs.iloc[indexes].reset_index(drop=True)
        totals=np.zeros(len(indexes),np.float64);detected=np.zeros(len(indexes),int)
        genes=[];values=[];all_features=[]
        for line in f:
            gene,raw=line.rstrip('\n').split('\t',1);gene=gene.strip('"')
            vector=np.fromstring(raw,sep='\t',dtype=np.float32)
            assert len(vector)==len(ids)
            assert np.isfinite(vector).all() and (vector>=0).all() and np.equal(vector,np.floor(vector)).all()
            v=vector[indexes];totals+=v;detected+=v>0;all_features.append(gene)
            if gene in needed:genes.append(gene);values.append(v)
    assert len(all_features)==len(set(all_features)) and len(obs)==3240
    obs['total_umi']=totals;obs['detected_genes']=detected
    assert (totals>0).all()
    x=np.column_stack(values)
    return x,obs,genes,{'source_features':len(all_features),'scoring_features':len(genes),
        'mapped_mouse_cells':len(obs),'mapped_mice':obs.unit.nunique(),'expression_type':'raw integer UMI counts',
        'feature_alignment':'source cell header and row symbols; all input counts checked'}


def load_development(config):
    manifest=json.loads((BASE/'developmental_expression_manifest.json').read_text())
    assert manifest['frozen_program_sha256']==sha(BASE/'frozen_repair_program.json')
    feature=pd.read_csv(TABLES/'developmental_feature_index.csv').set_index('feature_index')
    obs=decode_annotations((BASE/'cache/sources/Negretti_obs.bin').read_bytes())
    matrices=[];genes=[]
    for record in manifest['chunks']:
        path=BASE/record['file'];assert sha(path)==record['sha256']
        x,cols,rows=decode_numeric_packet(path.read_bytes())
        assert x.shape[0]==len(obs) and np.array_equal(rows,np.arange(len(obs)))
        assert cols is not None
        names=feature.loc[np.asarray(cols,dtype=int),'gene'].tolist()
        assert set(names)==set(record['requested_genes'])
        assert np.isfinite(x).all() and (x>=0).all()
        matrices.append(x);genes.extend(names)
    assert len(genes)==len(set(genes))==manifest['n_requested_features']
    x=np.column_stack(matrices)
    obs=obs.rename(columns={'name_0':'cell_id','celltype':'state','timepoint':'time'})
    obs['unit']=obs.time+':'+obs.cell_id.str.split('_',n=1).str[1]
    keep=obs.time.isin(config['transfer']['developmental_ages']).to_numpy()
    return x[keep],obs.loc[keep].reset_index(drop=True),genes,{
        'source_features':len(feature),'scoring_features':len(genes),'postnatal_cells':int(keep.sum()),
        'encoded_groups':int(obs.loc[keep,'unit'].nunique()),'independent_animals':'not_recovered',
        'expression_type':'author SCT normalized data-slot export; not raw counts',
        'feature_alignment':'every numeric API column mapped through source var index; every row index reconciled to annotation order'}


def analyze(dataset,x,obs,genes,states,modules,config):
    score_arrays={};effects=[];coverage=[];detect=[]
    lookup={g:i for i,g in enumerate(genes)}
    for name,module in modules.items():
        if name=='candidate_without_generic_controls' and len(module)<10:continue
        scores,c=rank_scores(x,genes,module,config['reference_genes'])
        coverage.append({'dataset':dataset,'module':name,**c})
        if name=='candidate':
            assert c['feature_coverage']>=.8 and c['reference_coverage']>=.8
        score_arrays[name]=scores
        result=paired_effects(obs,scores,start=states[0],intermediate=states[1],end=states[2],dataset=dataset)
        result['module']=name;effects.append(result)
        present=[lookup[g] for g in module if g in lookup]
        fraction=(x[:,present]>0).mean(axis=1)
        for (unit,state),f in obs[obs.state.isin(states)].groupby(['unit','state']):
            detect.append({'dataset':dataset,'module':name,'unit':unit,'state':state,'cells':len(f),
                'mean_fraction_of_present_module_genes_detected':float(fraction[f.index].mean())})
    effect=pd.concat(effects,ignore_index=True)
    effect.to_csv(TABLES/f'{dataset}_transfer_effects.csv',index=False)
    pd.DataFrame(coverage).to_csv(TABLES/f'{dataset}_transfer_coverage.csv',index=False)
    pd.DataFrame(detect).to_csv(TABLES/f'{dataset}_module_detection.csv',index=False)
    inventory=obs.groupby(['unit','time','state']).size().rename('cells').reset_index()
    inventory.to_csv(TABLES/f'{dataset}_state_inventory.csv',index=False)
    np.savez_compressed(PROC/f'{dataset}_transfer_scores.npz',**score_arrays)
    np.savez_compressed(PROC/f'{dataset}_scoring_matrix.npz',matrix=x)
    save_json(PROC/f'{dataset}_scoring_genes.json',genes)
    obs.to_csv(PROC/f'{dataset}_scored_obs.csv',index=False)
    primary=effect[effect.module=='candidate']
    summary={}
    for endpoint,f in primary.groupby('endpoint'):
        q=f[f.all_states_at_least_30]
        summary[endpoint]={'all_groups':len(f),'positive_all_groups':int((f.difference>0).sum()),
            'qualified_groups':len(q),'positive_qualified_groups':int((q.difference>0).sum()),
            'qualified_effects':q[['unit','difference','standardized_difference']].to_dict('records')}
    return {'coverage':next(x for x in coverage if x['module']=='candidate'),'effects':summary,
        'unit_interpretation':'verified mice; two pass floor' if dataset=='intestine' else 'encoded capture groups, not verified animals; one passes floor'}


def main():
    assert json.loads((BASE/'stage_E2_result.json').read_text())['checkpoint_passed']
    config=json.loads((BASE/'exploratory_config.json').read_text());modules=get_modules(config)
    result={'stage':'E3','program_retuned':False,'results':{},'input_audits':{}}
    x,obs,genes,audit=load_intestine(config,modules)
    result['results']['intestine']=analyze('intestine',x,obs,genes,config['transfer']['intestinal_states'],modules,config)
    result['input_audits']['intestine']=audit
    print('Intestinal transfer complete; no gene or state definitions changed.',flush=True)
    x,obs,genes,audit=load_development(config)
    result['results']['development']=analyze('development',x,obs,genes,config['transfer']['developmental_states'],modules,config)
    result['input_audits']['development']=audit
    result['interpretation_status']='requires checkpoint review; descriptive transfer does not establish broad conservation'
    save_json(BASE/'stage_E3_result.json',result)
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':
    main()
