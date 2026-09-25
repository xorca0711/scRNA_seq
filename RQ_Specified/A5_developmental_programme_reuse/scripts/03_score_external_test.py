"""Apply fixed raw-depth gates and compute expected depth-matched detection scores."""
import argparse
import gzip
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]


def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1 << 20),b''):
            h.update(b)
    return h.hexdigest()


def detection_probability(N, k, budget):
    """Exact expectation of without-replacement downsampling, evaluated in log space."""
    import numpy as np
    from scipy.special import gammaln
    N, k = np.broadcast_arrays(np.asarray(N,dtype=float),np.asarray(k,dtype=float))
    assert np.all(N >= budget) and np.all(k >= 0) and np.all(k <= N)
    p = np.ones_like(N)
    ok = N-k >= budget
    nn, kk = N[ok], k[ok]
    log_absence = (gammaln(nn-kk+1)-gammaln(nn-kk-budget+1)
                   -gammaln(nn+1)+gammaln(nn-budget+1))
    p[ok] = -np.expm1(np.minimum(log_absence,0))
    p[k==0] = 0
    assert np.isfinite(p).all() and ((p>=0)&(p<=1)).all()
    return p


def main():
    import numpy as np
    import pandas as pd
    import scipy
    from scipy.io import mmread
    from scipy.sparse import coo_matrix
    ap = argparse.ArgumentParser()
    ap.add_argument('--source-root',type=Path,required=True)
    args = ap.parse_args()
    out = HERE/'tables/test_v1'
    out.mkdir(exist_ok=True)
    if any(out.iterdir()):
        raise SystemExit('Refusing to overwrite A5 results')
    config_path = HERE/'config/strunz_test_contract.json'
    cfg = json.loads(config_path.read_text())
    frozen_path = HERE/'tables/external_test_modules.json'
    modules = json.loads(frozen_path.read_text())['modules']
    frozen = json.loads((ROOT/'RQ_Specified/A5_A11_shared_component_contract/tables/frozen_modules.json').read_text())
    for k in ['development_specific','shared_remodelling','lesion_specific']:
        modules['contract_'+k] = frozen['modules'][k]['genes']
    member = pd.read_csv(ROOT/'RQ_Specified/A5_A11_shared_component_contract/tables/module_membership.tsv',sep='\t')
    d = member[member.module=='development_specific']
    modules['Strunz_filtered_development_51'] = d.loc[~d.in_AT1_published_400 & ~d.in_AT2_published_400,'gene'].tolist()
    modules['contract_shared_D_I'] = member.loc[member.in_D & member.in_I,'gene'].tolist()
    modules['contract_shared_I_L'] = member.loc[member.in_I & member.in_L,'gene'].tolist()
    src = args.source_root/'RQ_Specified/A5_developmental_programme_reuse/cache'
    metadata = src/'GSE141259_HighResolution_cellinfo.csv.gz'
    gene_path = src/'GSE141259_HighResolution_genes.txt.gz'
    audit = json.loads((HERE/'tables/audit_run.json').read_text())
    for p in [metadata,gene_path]:
        assert sha(p)==audit['sources'][p.name]['sha256']
    cache=HERE/'cache'
    matrix_path=cache/'GSE141259_HighResolution_rawcounts.mtx.gz'
    barcode_path=cache/'GSE141259_HighResolution_barcodes.txt.gz'
    download=json.loads((HERE/'tables/count_retrieval.json').read_text())
    for p in [matrix_path,barcode_path]:
        assert sha(p)==download['files'][p.name]['sha256']
    genes=[s.strip() for s in gzip.open(gene_path,'rt')]
    barcodes=[s.strip() for s in gzip.open(barcode_path,'rt')]
    assert len(set(genes))==len(genes) and len(set(barcodes))==len(barcodes)
    c=pd.read_csv(metadata,sep='\t').set_index('cell_barcode')
    assert c.index.is_unique and set(c.index)==set(barcodes)
    c=c.loc[barcodes].copy()
    assert c.groupby('identifier').sample_id.nunique().max()==1
    assert c.groupby('sample_id').identifier.nunique().max()==1
    c['day']=c.time_point.str.extract(r'(\d+)').astype(int)
    print('Reading raw Strunz matrix',flush=True)
    with gzip.open(matrix_path,'rb') as f:
        x=mmread(f).tocsc()
    x.sum_duplicates(); x.eliminate_zeros()
    assert x.shape==(len(genes),len(barcodes))
    assert np.isfinite(x.data).all() and (x.data>=0).all() and (x.data==np.floor(x.data)).all()
    library=np.asarray(x.sum(axis=0)).ravel()
    c['raw_umi']=library
    c['depth_pass']=library>=cfg['depth_umi']
    in_window=c.day.between(cfg['day_min'],cfg['day_max']) & ~c.sample_id.str.startswith(cfg['exclude_sample_prefix'])
    labels=[cfg['case_label'],cfg['reference_label'],cfg['secondary_reference_label']]
    groups=[]
    for (sample,day,label),q in c[in_window & c.cell_type.isin(labels)].groupby(['sample_id','day','cell_type']):
        groups.append({'sample_id':sample,'day':int(day),'label':label,'cells_before_depth':len(q),
                       'cells_after_depth':int(q.depth_pass.sum()),'median_raw_umi':float(q.raw_umi.median())})
    units=pd.DataFrame(groups)
    units.to_csv(out/'cell_depth_gates.tsv',sep='\t',index=False)
    elig={}
    for ref in labels[1:]:
        a=set(units.loc[(units.label==cfg['case_label']) & (units.cells_after_depth>=cfg['cell_floor']),'sample_id'])
        b=set(units.loc[(units.label==ref) & (units.cells_after_depth>=cfg['cell_floor']),'sample_id'])
        elig[ref]=sorted(a&b)
    coverage=[]; universe=set(genes)
    for k,gs in modules.items():
        present=sorted(set(gs)&universe)
        coverage.append({'module':k,'source_genes':len(gs),'assayed_genes':len(present),
                         'fraction':len(present)/len(gs),'eligible':len(present)/len(gs)>=cfg['coverage_floor']})
        modules[k]=present
    cov=pd.DataFrame(coverage);cov.to_csv(out/'coverage.tsv',sep='\t',index=False)
    gate={'eligible_mice':elig,'primary_n':len(elig[cfg['reference_label']]),
          'primary_coverage':bool(cov.set_index('module').loc[cfg['primary_module'],'eligible']),
          'count_shape':list(x.shape),'raw_count_sum':int(library.sum()),'depth_umi':cfg['depth_umi']}
    (out/'gates.json').write_text(json.dumps(gate,indent=2)+'\n')
    if gate['primary_n']<cfg['unit_floor'] or not gate['primary_coverage']:
        print('Primary ineligible; no scores',gate,flush=True)
        return
    take=np.flatnonzero(in_window.to_numpy() & c.cell_type.isin(labels).to_numpy() & c.depth_pass.to_numpy())
    union=sorted(set().union(*map(set,modules.values())))
    gene_index={g:i for i,g in enumerate(genes)}
    small=x[[gene_index[g] for g in union],:][:,take].tocoo()
    values=detection_probability(library[take][small.col],small.data,cfg['depth_umi'])
    prob=coo_matrix((values,(small.row,small.col)),shape=small.shape).tocsr()
    del x, small
    selected=c.iloc[take]; uidx={g:i for i,g in enumerate(union)}
    score_rows=[]
    for k,gs in modules.items():
        if not len(gs): continue
        cell_scores=np.asarray(prob[[uidx[g] for g in gs],:].sum(axis=0)).ravel()/len(gs)*100
        frame=selected[['sample_id','day','cell_type']].copy();frame['score']=cell_scores
        for (sample,day,label),q in frame.groupby(['sample_id','day','cell_type']):
            score_rows.append({'sample_id':sample,'day':int(day),'label':label,'module':k,'cells':len(q),'score':float(q.score.mean())})
    scores=pd.DataFrame(score_rows); scores.to_csv(out/'unit_scores.tsv',sep='\t',index=False)
    pairs=[]
    for ref,samples in elig.items():
        for sample in samples:
            for k in modules:
                a=scores[(scores.sample_id==sample)&(scores.module==k)&(scores.label==cfg['case_label'])]
                b=scores[(scores.sample_id==sample)&(scores.module==k)&(scores.label==ref)]
                assert len(a)==1 and len(b)==1
                pairs.append({'sample_id':sample,'day':int(a.iloc[0].day),'reference':ref,'module':k,
                              'difference_pp':float(a.iloc[0].score-b.iloc[0].score)})
    pd.DataFrame(pairs).to_csv(out/'paired_differences.tsv',sep='\t',index=False)
    record={'completed_utc':datetime.now(timezone.utc).isoformat(),'preregistration_commit':'eb5317e',
            'code_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
            'script_sha256':sha(Path(__file__)),'config_sha256':sha(config_path),
            'inputs':{str(p):sha(p) for p in [metadata,gene_path,matrix_path,barcode_path,frozen_path]},
            'outputs':{p.name:sha(p) for p in out.iterdir()},'software':{'numpy':np.__version__,'scipy':scipy.__version__,'pandas':pd.__version__}}
    (out/'score_run.json').write_text(json.dumps(record,indent=2)+'\n')
    print('A5 scored; complete mice by reference:',{k:len(v) for k,v in elig.items()},flush=True)


if __name__=='__main__':
    main()
