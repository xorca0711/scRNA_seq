"""ES1: frozen-panel specificity at a joint RNA depth budget; descriptive only."""
from __future__ import annotations
import argparse
import hashlib
import itertools
import json
import platform
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
OUT=HERE/'results'
BUDGET=2000
SEEDS=(17,29)
MIN_GROUP=30
AIRWAY=['Scgb1a1','Scgb3a2','Foxj1','Krt5']
ALVEOLAR=['Sftpc','Sftpb','Lamp3']
LABELS=['Cldn4','Krt8','Sftpc']

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as fh:
        for b in iter(lambda:fh.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()

def pathmeta(path):
    return dict(path=path.relative_to(ROOT).as_posix(),bytes=path.stat().st_size,
                sha256=sha(path))

def decode(x):
    return x.decode() if isinstance(x,bytes) else str(x)

def read_10x(path,genes):
    """One bounded sparse column pass; retain selected RNA rows, not peaks."""
    with h5py.File(path,'r') as h:
        m=h['matrix']; f=m['features']
        names=[decode(x) for x in f['name'][:]]
        types=np.array([decode(x) for x in f['feature_type'][:]])
        rna=types=='Gene Expression'
        slot=np.full(len(names),-1,dtype=int)
        lookup={g:i for i,g in enumerate(genes)}
        for i,g in enumerate(names):
            if rna[i] and g in lookup:slot[i]=lookup[g]
        bc=[decode(x) for x in m['barcodes'][:]]
        ptr=m['indptr'][:]; n=len(bc)
        selected=np.zeros((len(genes),n),dtype=np.int32); totals=np.zeros(n,dtype=np.int64)
        for start in range(0,n,2000):
            end=min(n,start+2000);lo,hi=ptr[start],ptr[end]
            idx=m['indices'][lo:hi];d=m['data'][lo:hi]
            assert np.all(d>=0) and np.all(d==np.floor(d))
            col=np.repeat(np.arange(start,end),np.diff(ptr[start:end+1]))
            keep=rna[idx]
            totals[start:end]=np.bincount(col[keep]-start,weights=d[keep],minlength=end-start).astype(np.int64)
            target=slot[idx];hit=target>=0
            # Gene symbols with multiple feature IDs are one measured gene;
            # sum their raw counts instead of choosing one feature arbitrarily.
            np.add.at(selected,(target[hit],col[hit]),d[hit])
        return selected,totals,np.array([x.rsplit('-',1)[1] for x in bc]),set(g for g,t in zip(names,rna) if t)

def joint_sample(counts,totals,budget,seed):
    """Selected categories plus residual follow one multivariate hypergeometric."""
    assert np.all(totals>=budget) and np.all(counts.sum(axis=0)<=totals)
    rng=np.random.default_rng(seed)
    remaining=totals.copy();draws=np.full(len(totals),budget,dtype=np.int64)
    sampled=np.zeros_like(counts)
    for j,row in enumerate(counts):
        mask=(row>0)&(draws>0)
        sampled[j,mask]=rng.hypergeometric(row[mask],remaining[mask]-row[mask],draws[mask])
        remaining-=row
        draws-=sampled[j]
    assert np.all(sampled<=counts) and np.all(sampled.sum(axis=0)<=budget)
    assert np.all(draws<=remaining)
    return sampled

def score_unit(counts,totals,genes,present,meta,modules,rows,gene_rows,retention):
    pos={g:i for i,g in enumerate(genes)}
    airway=counts[[pos[g] for g in AIRWAY]].sum(axis=0)
    alv=counts[[pos[g] for g in ALVEOLAR]].sum(axis=0)
    keep=(totals>=BUDGET)&(airway<=alv)
    retention.append(dict(meta,n_input=len(totals),n_depth=int((totals>=BUDGET).sum()),
        n_retained=int(keep.sum()),retained_fraction=float(keep.mean()) if len(keep) else np.nan,
        median_UMI=float(np.median(totals)) if len(totals) else np.nan))
    for seed in SEEDS:
        ds=joint_sample(counts[:,keep],totals[keep],BUDGET,seed)
        trans=(ds[pos['Cldn4']]>0)&(ds[pos['Krt8']]>0)
        masks={'two_marker':trans,'reference':(ds[pos['Sftpc']]>0)&~trans}
        for group,mask in masks.items():
            n=int(mask.sum())
            det=(ds[:,mask]>0).mean(axis=1) if n else np.full(len(genes),np.nan)
            for module in modules:
                gs=[g for g in module['genes'] if g in present]
                coverage=len(gs)/len(module['genes'])
                rows.append(dict(meta,seed=seed,group=group,n_cells=n,module=module['name'],
                    kind=module['kind'],n_genes=len(gs),total_genes=len(module['genes']),coverage=coverage,
                    score=float(np.mean(det[[pos[g] for g in gs]])) if gs else np.nan,
                    interpretable=bool(n>=MIN_GROUP and coverage>=.8),
                    label_overlap=';'.join(module.get('label_overlap',[]))))
            # Gene-level evidence is retained for the primary technical draw.
            if seed==SEEDS[0]:
                for g,d in zip(genes,det):
                    gene_rows.append(dict(meta,group=group,n_cells=n,gene=g,present=g in present,detection=float(d)))

def external(genes,modules,rows,gene_rows,retention,inputs):
    """Read raw post-QC counts, not the log-normalized epithelial cache."""
    import anndata as ad
    path=ROOT/'Thesis/gate1_01_niethamer_2025/GSE262927/processed/postQC.h5ad'
    units=ROOT/'Thesis/gate1_01_niethamer_2025/trials/w1_amac_pseudobulk_de/w1_units.csv'
    if not path.exists():return {'status':'unavailable','reason':'postQC raw counts missing'}
    inputs.append(pathmeta(path));inputs.append(pathmeta(units))
    um=pd.read_csv(units).set_index('animal')
    with h5py.File(path,'r') as h:
        obs=ad.io.read_elem(h['obs']);var=ad.io.read_elem(h['var'])
        names=var.index.astype(str).tolist()
        labels=obs['author_celltype'].astype(str)
        # Fixed named alveolar cell types, not a score-fitted or trajectory-derived selection.
        allowed={'AT1','AT2','AT1_AT2','Alveolar_transitional'}
        mask=obs['has_author_metadata'].astype(bool).to_numpy() & labels.isin(allowed).to_numpy()
        sel=np.flatnonzero(mask); sub=obs.iloc[sel].copy()
        m=h['X'];ptr=m['indptr'][:]
        lookup={g:i for i,g in enumerate(genes)};slot=np.array([lookup.get(g,-1) for g in names])
        counts=np.zeros((len(genes),len(sel)),dtype=np.int32);totals=np.zeros(len(sel),dtype=np.int64)
        # Read only contiguous blocks containing selected cells, then slice rows in memory.
        for start in range(0,len(obs),2000):
            end=min(len(obs),start+2000);where=np.flatnonzero((sel>=start)&(sel<end))
            if not len(where):continue
            lo,hi=ptr[start],ptr[end]
            d=m['data'][lo:hi];idx=m['indices'][lo:hi]
            assert np.all(d>=0) and np.all(d==np.floor(d)),'external object is not raw counts'
            block=sp.csr_matrix((d,idx,ptr[start:end+1]-lo),shape=(end-start,len(names)))[sel[where]-start]
            totals[where]=np.asarray(block.sum(axis=1)).ravel().astype(np.int64)
            valid=np.flatnonzero(slot>=0)
            counts[slot[valid,None],where]=block[:,valid].T.toarray().astype(np.int32)
    inventory=sub.groupby(['sample_id','sacrifice_day','sex','author_celltype'],observed=True).size().rename('n_cells').reset_index()
    inventory.to_csv(OUT/'external_population_inventory.csv',index=False)
    for animal,idx in sub.groupby('sample_id',observed=True).indices.items():
        idx=np.asarray(idx);o=sub.iloc[idx[0]]
        genotype=str(um.loc[animal,'genotype']) if animal in um.index else 'metadata unavailable'
        meta=dict(cohort='GSE262927',unit=str(animal),unit_type='animal',genotype=genotype,
                  stage=str(o['sacrifice_day']),injury=str(o['condition']),
                  independence='external-study, previously analyzed; no neonatal arm')
        score_unit(counts[:,idx],totals[idx],genes,set(names),meta,modules,rows,gene_rows,retention)
    return {'status':'scored','animal_count':int(sub.sample_id.nunique()),'input_cells':len(sel),
            'allowed_author_labels':sorted(allowed),'available_author_labels':sorted(labels.unique()),
            'caution':'external to multiome studies, not untouched held-out data; different injury; no neonatal arm; cell floor decides evaluability'}

def write_summary(scores,effects,ext):
    main=effects[(effects.seed==17)&(effects.module=='ADI_published_holdout')]
    second=effects[(effects.seed==29)&(effects.module=='ADI_published_holdout')].set_index('unit')
    lines=['# ES1 results: epithelial-state specificity','',
     'Descriptive RNA measurements at 2,000 UMI per retained cell. No fate, reversibility, chromatin or causal result is inferred.',
     '',f'Frozen definitions: {len(scores.module.unique())} panels/modules, including the complete 400-gene published ADI, AT2 and AT1 marker lists and label-free variants. Complete DATP/PATS and developmental maturation signatures remain unavailable; short panels, published marker lists and full Hallmark modules are distinguished.',
     '', '## Complete published ADI list, excluding labeling genes','',
     'The original list contains 400 genes. Removing Cldn4 and Krt8 leaves 398; 390 are present in each multiome feature universe. All reported effects below are labelled-minus-reference detection points, not percent change in expression.',
     '', '| Unit | Labelled cells, seeds 17 / 29 | Difference, seeds 17 / 29 | Both seeds evaluable |',
     '|---|---:|---:|---|']
    for r in main[main.cohort!='GSE262927'].itertuples():
        s=second.loc[r.unit]
        lines.append(f'| {r.unit} | {r.n_two_marker} / {s.n_two_marker} | {100*r.difference:+.2f} / {100*s.difference:+.2f} | {bool(r.interpretable and s.interpretable)} |')
    p=main.set_index('unit')
    lines+=['',f'In control-genotype wells the ADI contrast is {100*p.loc["P9_control","difference"]:.2f} / {100*second.loc["P9_control","difference"]:.2f} points at P9 and {100*p.loc["SeV_control","difference"]:.2f} / {100*second.loc["SeV_control","difference"]:.2f} in injured adults. This is descriptive evidence that the same two-transcript call has different panel specificity across contexts; it does not identify an injury-specific program or a developmental mechanism.',
      '', 'The seven-week control crosses the cell floor between seeds (30 versus 26 labelled cells). Its contrast is not robustly evaluable. Wildtype PBS and AP1-mutant PBS do not reach the floor in either seed.']
    em=main[main.cohort=='GSE262927']
    lines+=['',f'External sample check: {ext.get("animal_count",0)} animals scored; {int(em.interpretable.sum())} meet both 30-cell group floors for this panel.',
     'This is an external study consistency check on previously analyzed data, not independent confirmation of developmental-program reuse. No neonatal samples are present.']
    for r in em[em.interpretable].itertuples():
        s=second.loc[r.unit]
        lines.append(f'{r.unit} at {r.stage} dpi has {r.n_two_marker}/{r.n_reference} labelled/reference cells in seed 17 and {s.n_two_marker}/{s.n_reference} in seed 29, with ADI differences {100*r.difference:+.2f}/{100*s.difference:+.2f} points. A single eligible animal does not provide replicated confirmation.')
    lines += [
     '', '## Interpretation limits','',
     '- Full DATP/PATS and developmental maturation signatures and a development-by-injury factorial cohort are still missing. The ADI source caps its published list at 400 markers; no marker list is an exhaustive program.',
     '- Each multiome condition is one pooled library. Technical seeds are sensitivity checks, never replicate animals.',
     '- The full five-gene DATP marker panel overlaps the label and is only a circular control. Primary specificity readings use label-free panels.',
     '- The reference requires Sftpc detection and thus operationally enriches AT2 cells; these are conditional group contrasts, not unbiased state prevalence.',
     '- A score within one panel is comparable across wells at this budget; absolute scores across different panels are not comparable measures of program strength.',
     '- Low cell coverage means insufficient measurement, not absence of a biological program.',
     '', '## Files','',
     '`module_scores.csv`, `within_unit_effects.csv`, `between_well_contrasts.csv`, `retention.csv`, `module_overlap.csv`, `seed_sensitivity.csv`, `external_population_inventory.csv`, `verification.json`, `run_record.json` and `es1_specificity.png` are compact tracked outputs. The larger `gene_detection.csv` is a regenerable local artifact; source XLSX and plot PDF are also untracked. Source factual gene lists and hashes are preserved in `modules.json`.']
    (OUT/'SUMMARY.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')

def plot(effects):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    choose=['DATP_label_free_panel','DATP_PATS_holdout_panel','ADI_published_holdout','AT2_published_holdout',
            'AT1_published_holdout','AT1_late_panel','HALLMARK_P53_PATHWAY','HALLMARK_HYPOXIA','HALLMARK_INFLAMMATORY_RESPONSE']
    e=effects[(effects.seed==17)&(effects.cohort!='GSE262927')]
    order=['P9_control','P9_Cebpa_mutant','7wk_control','7wk_Cebpa_mutant','SeV_control','SeV_Cebpa_mutant',
           'wildtype_PBS','wildtype_SeV','AP1mut_PBS','AP1mut_SeV']
    p=e.pivot(index='unit',columns='module',values='difference').reindex(index=order,columns=choose)*100
    ok=e.pivot(index='unit',columns='module',values='interpretable').reindex(index=order,columns=choose)
    data=p.to_numpy().copy();data[~ok.fillna(False).to_numpy(dtype=bool)]=np.nan
    fig,(ax,bx)=plt.subplots(1,2,figsize=(15,7),gridspec_kw={'width_ratios':[2.6,1]},layout='constrained')
    cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#e4e4e4')
    limit=5*np.ceil(np.nanmax(np.abs(data))/5)
    im=ax.imshow(data,cmap=cmap,vmin=-limit,vmax=limit,aspect='auto')
    labels=['P9 control','P9 Cebpa mutant','7wk control *','7wk Cebpa mutant','Injured control','Injured Cebpa mutant',
            'WT PBS','WT injured','AP1 mutant PBS','AP1 mutant injured']
    ax.set_yticks(range(len(order)),labels)
    names=['DATP\n3 holdout genes','DATP/PATS\n6 holdout genes','ADI\npublished list','AT2\npublished list','AT1\npublished list',
           'AT1 late\n4 genes','P53\nfull Hallmark','Hypoxia\nfull Hallmark','Inflammation\nfull Hallmark']
    ax.set_xticks(range(len(choose)),names,rotation=40,ha='right')
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            ax.text(j,i,'<30 cells' if np.isnan(data[i,j]) else f'{data[i,j]:+.1f}',ha='center',va='center',fontsize=7,
                    color='white' if not np.isnan(data[i,j]) and abs(data[i,j])>10 else '#222222')
    ax.set_title('A  Within-well RNA differences\n* 7wk control falls below the cell floor in seed 29',loc='left')
    fig.colorbar(im,ax=ax,shrink=.7,label='Labelled minus reference (detection points)')
    ext=effects[(effects.cohort=='GSE262927')&(effects.seed==17)&(effects.module=='ADI_published_holdout')]
    ext=ext.copy();ext['day']=pd.to_numeric(ext.stage);ext=ext.sort_values(['day','unit'])
    good=ext[ext.interpretable]
    other=effects[(effects.cohort=='GSE262927')&(effects.seed==29)&(effects.module=='ADI_published_holdout')].set_index('unit')
    assert (ext.n_reference>=MIN_GROUP).all() and (other.n_reference>=MIN_GROUP).all()
    x=np.arange(len(ext))
    bx.barh(x-.16,ext.n_two_marker,height=.3,color='#2166ac',label='Seed 17')
    bx.barh(x+.16,other.loc[ext.unit,'n_two_marker'],height=.3,color='#92c5de',label='Seed 29')
    bx.set_yticks(x,[f'{r.unit.replace("EEM-scRNA-","")} / d{r.day:g}' for r in ext.itertuples()],fontsize=7)
    bx.axvline(MIN_GROUP,color='#b2182b',ls='--',lw=1,label='30-cell floor')
    bx.invert_yaxis();bx.legend(fontsize=7,loc='lower right')
    bx.set_title(f'B  External-study coverage\n{len(good)} / {len(ext)} animals evaluable',loc='left')
    bx.set_xlabel('Two-marker positive cells per animal\nAll reference groups pass the 30-cell floor')
    fig.suptitle('Epithelial-state specificity at 2,000 RNA UMI\nMarker panels do not identify developmental program identity or future fate',fontsize=14)
    fig.savefig(OUT/'es1_specificity.png',dpi=180)
    fig.savefig(OUT/'es1_specificity.pdf')
    plt.close(fig)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--replot',action='store_true',help='Regenerate the figure from frozen score tables only')
    args=parser.parse_args()
    if args.replot:
        path=OUT/'run_record.json';record=json.loads(path.read_text(encoding='utf-8'))
        effects=pd.read_csv(OUT/'within_unit_effects.csv')
        plot(effects)
        write_summary(pd.read_csv(OUT/'module_scores.csv'),effects,record['external'])
        record['replot']={'utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
            'script_sha256':sha(Path(__file__)),
            'scope':'Figure and table-derived summary only; numerical tables unchanged. Includes published ADI/AT2/AT1 list panels.'}
        record['outputs']=[pathmeta(p) for p in sorted(OUT.iterdir()) if p.is_file() and p.name!='run_record.json']
        path.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
        print('Replotted frozen results; numeric tables unchanged')
        return
    OUT.mkdir(exist_ok=True)
    config=json.loads((HERE/'modules.json').read_text(encoding='utf-8'));modules=config['modules']
    genes=sorted(set(g for m in modules for g in m['genes'])|set(AIRWAY+ALVEOLAR+LABELS))
    overlaps=[]
    for a,b in itertools.combinations(modules,2):
        x,y=set(a['genes']),set(b['genes'])
        overlaps.append(dict(module_a=a['name'],module_b=b['name'],intersection=len(x&y),jaccard=len(x&y)/len(x|y),genes=';'.join(sorted(x&y))))
    pd.DataFrame(overlaps).to_csv(OUT/'module_overlap.csv',index=False)
    inputs=[pathmeta(HERE/'PLAN.md'),pathmeta(HERE/'modules.json')];rows=[];gene_rows=[];retention=[]
    descriptors=[('GSE247130','GSE247130_Aggregate.P9_filtered_feature_bc_matrix.h5','P9',
                  ['P9_control','P9_Cebpa_mutant'],['control','Cebpa_mutant'],'uninjured'),
                 ('GSE247130','GSE247130_Aggregate7-wk_filtered_feature_bc_matrix.h5','7wk',
                  ['7wk_control','7wk_Cebpa_mutant'],['control','Cebpa_mutant'],'uninjured'),
                 ('GSE247130','GSE247130_Aggregate_SeV_filtered_feature_bc_matrix.h5','adult_14dpi',
                  ['SeV_control','SeV_Cebpa_mutant'],['control','Cebpa_mutant'],'SeV'),
                 ('GSE310539','GSE310539_totalaggr_filtered_feature_bc_matrix.h5','adult',
                  ['wildtype_PBS','wildtype_SeV','AP1mut_PBS','AP1mut_SeV'],
                  ['wildtype','wildtype','AP1mut','AP1mut'],None)]
    for dep,file,stage,units,genotypes,inj in descriptors:
        path=ROOT/'raw_data'/dep/file;print('Reading',file,flush=True)
        inputs.append(pathmeta(path))
        counts,totals,suffix,present=read_10x(path,genes)
        for i,(unit,genotype) in enumerate(zip(units,genotypes),1):
            sel=suffix==str(i)
            meta=dict(cohort=dep,unit=unit,unit_type='pooled_library',genotype=genotype,stage=stage,
                      injury=inj or ('SeV' if unit.endswith('_SeV') else 'uninjured'),
                      independence='same-laboratory studies; not independent replication')
            score_unit(counts[:,sel],totals[sel],genes,present,meta,modules,rows,gene_rows,retention)
        del counts
    print('Checking external raw counts and animal groups',flush=True)
    ext=external(genes,modules,rows,gene_rows,retention,inputs)
    scores=pd.DataFrame(rows);scores.to_csv(OUT/'module_scores.csv',index=False)
    pd.DataFrame(gene_rows).to_csv(OUT/'gene_detection.csv',index=False)
    pd.DataFrame(retention).to_csv(OUT/'retention.csv',index=False)
    keys=['cohort','unit','unit_type','genotype','stage','injury','seed','module','kind']
    effects=scores[scores.group=='two_marker'].merge(scores[scores.group=='reference'],on=keys,suffixes=('_two_marker','_reference'))
    effects['difference']=effects.score_two_marker-effects.score_reference
    effects['interpretable']=effects.interpretable_two_marker & effects.interpretable_reference
    effects=effects.rename(columns={'n_cells_two_marker':'n_two_marker','n_cells_reference':'n_reference'})
    effects.to_csv(OUT/'within_unit_effects.csv',index=False)
    contrasts=[]
    pairs=[('P9_control','7wk_control','development_within_control'),('P9_Cebpa_mutant','7wk_Cebpa_mutant','development_within_mutant'),
       ('SeV_control','7wk_control','adult_injury_within_control'),('SeV_Cebpa_mutant','7wk_Cebpa_mutant','adult_injury_within_mutant'),
       ('P9_Cebpa_mutant','P9_control','genotype_P9'),('7wk_Cebpa_mutant','7wk_control','genotype_7wk'),
       ('SeV_Cebpa_mutant','SeV_control','genotype_injured'),('wildtype_SeV','wildtype_PBS','injury_wildtype'),('AP1mut_SeV','AP1mut_PBS','injury_AP1mut')]
    for a,b,name in pairs:
        x=effects[effects.unit==a].merge(effects[effects.unit==b],on=['seed','module'],suffixes=('_a','_b'))
        for r in x.itertuples():
            contrasts.append(dict(contrast=name,unit_a=a,unit_b=b,module=r.module,seed=r.seed,
                delta_within_well_effect=r.difference_a-r.difference_b,
                both_evaluable=r.interpretable_a and r.interpretable_b,
                inference='descriptive only; one pooled library per condition'))
    pd.DataFrame(contrasts).to_csv(OUT/'between_well_contrasts.csv',index=False)
    write_summary(scores,effects,ext);plot(effects)
    versions={m.__name__:m.__version__ for m in [np,pd,h5py]}
    record={'completed_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'python':sys.version,'platform':platform.platform(),
      'versions':versions,'script_sha256':sha(Path(__file__)),'budget':BUDGET,'seeds':SEEDS,'minimum_group':MIN_GROUP,
      'inputs':inputs,'external':ext,'module_provenance':config,'outputs':[]}
    for path in sorted(OUT.iterdir()):
        if path.is_file() and path.name!='run_record.json':record['outputs'].append(pathmeta(path))
    (OUT/'run_record.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    print('Completed',len(scores),'module/group rows;',ext,flush=True)

if __name__=='__main__':main()
