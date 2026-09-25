"""Descriptive reconstruction from measured tracing labels and author RNA states."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE=Path(__file__).resolve().parents[1]


def sha(p):
    with p.open('rb') as h:return hashlib.file_digest(h,'sha256').hexdigest()


def main():
    out=BASE/'tables/hpcs_source_composition';figdir=BASE/'figures/hpcs_source_composition'
    rp=BASE/'reports/hpcs_source_composition_run.json'
    if any(p.exists() for p in [out,figdir,rp]):raise SystemExit('Refusing to overwrite source reconstruction')
    cp=BASE/'config/hpcs_descendant_reconstruction.json';cfg=json.loads(cp.read_text())
    recovery=json.loads((BASE/'reports/hpcs_metadata_recovery.json').read_text())
    p=BASE/cfg['input'];assert sha(p)==recovery['output_sha256']
    obs=pd.read_csv(p,index_col=0,low_memory=False)
    assert obs.index.is_unique and len(obs)==28402
    data=obs[obs.batch.isin(cfg['select_batches'])].copy()
    assert len(data)==5333 and set(data.Group)==set(cfg['groups'])
    assert set(data.Sorting_Groups)=={cfg['require_sorting_group']}
    assert data.Reporter.isna().all() and data['cell type'].notna().all()
    expected={'6wk_3d':586,'8wk':742,'12wk_3d':877,'14wk':991,'Hopx_12wk3d':1226,'Hopx_12wk14d':911}
    assert data.groupby('Group').size().to_dict()==expected
    source_hashes={cfg['input']:sha(p),'config/hpcs_descendant_reconstruction.json':sha(cp)}
    geo=json.loads((BASE/'metadata/GSE277777.json').read_text())
    source_rows=[];rows=[]
    states=['HPCS','AT1-like','AT2-like','Lung endoderm-like','Hybrid lung/gastric-like','Highly proliferative','EMT','Ribosome']
    assert set(data['cell type'])==set(states)
    for group,info in cfg['groups'].items():
        notebooks={}
        for name in info['notebooks']:
            path=next((BASE/'cache/followup_sources/hpcs/tracing'/name).glob('*.ipynb'))
            n=json.loads(path.read_text());source_hashes[path.relative_to(BASE).as_posix()]=sha(path)
            # Only the early stored hash-feature outputs establish assignment.
            output='\n'.join(''.join(o.get('text',[]))+' '.join(o.get('data',{}).get('text/plain',[]))
                for cell in n['cells'][:36] for o in cell.get('outputs',[]))
            notebooks[name]=(path,output)
        for label,f in data[data.Group==group].groupby('Classification',sort=True):
            matches=[name for name,(path,text) in notebooks.items() if label in text]
            assert len(matches)==1,(label,matches)
            name=matches[0];path=notebooks[name][0]
            samples=[s for s in geo['samples'] if any(name in v and 'feature_bc_matrix' in v
                for k,values in s.items() if k.startswith('supplementary_file') for v in values)]
            assert len(samples)==1,(name,len(samples))
            sample=samples[0]
            source_rows.append(dict(source_label=label,group=group,driver=info['driver'],chase_days=info['chase_days'],
                source_age_label=info['source_age_label'],source_library=name,gsm=sample['accession'],
                GEO_title='; '.join(sample['title']),retained_cells=len(f),
                permanent_trace='author trace-sorted assignment',current_mScarlet='unavailable for traced rows',
                unit_status='source alias; independent mouse/pool membership not verified',
                source_notebook=path.relative_to(BASE).as_posix()))
            counts=f['cell type'].value_counts()
            for state in states:
                count=int(counts.get(state,0))
                rows.append(dict(source_label=label,group=group,driver=info['driver'],chase_days=info['chase_days'],
                    state=state,numerator=count,denominator=len(f),fraction=count/len(f)))
    table=pd.DataFrame(rows);manifest=pd.DataFrame(source_rows)
    assert len(manifest)==22 and len(table)==176
    assert manifest.source_label.is_unique
    assert table.groupby('source_label').numerator.sum().equals(manifest.set_index('source_label').retained_cells.sort_index())
    np.testing.assert_allclose(table.groupby('source_label').fraction.sum(),1)
    hpcs=table[table.state=='HPCS'].set_index('source_label')
    author_early={'BO1535_B0302':.894737,'BO1540_B0303':.938356,'BR1128_B0307':.964286,
        'BR1615_B0304':.924051,'BS1265_B0305':.932584,'BL1121_B0304':.736842,'BL1124_B0305':.761194,
        'BL1233_B0306':.465517,'BL1237_B0302':.782427,'BL1606_B0303':.570621,'BL1646_B0301':.754902}
    for label,value in author_early.items():np.testing.assert_allclose(hpcs.loc[label,'fraction'],value,atol=5.1e-7,rtol=0)
    for group,state,mean in [('Hopx_12wk3d','HPCS',.04767397154940407),
        ('Hopx_12wk14d','HPCS',.0432409972299169),('Hopx_12wk3d','AT1-like',.47159479260113635),
        ('Hopx_12wk14d','AT1-like',.6506824477461597)]:
        np.testing.assert_allclose(table[(table.group==group)&(table.state==state)].fraction.mean(),mean,atol=1e-12)
    summary=table.groupby(['group','driver','chase_days','state'],sort=False).agg(
        source_labels=('source_label','nunique'),total_cells=('denominator','sum'),state_cells=('numerator','sum'),
        source_mean_fraction=('fraction','mean'),source_min_fraction=('fraction','min'),source_max_fraction=('fraction','max')).reset_index()
    summary['cell_pooled_fraction']=summary.state_cells/summary.total_cells
    out.mkdir(parents=True);figdir.mkdir(parents=True)
    table.to_csv(out/'source_state_counts.tsv',sep='\t',index=False)
    manifest.to_csv(out/'source_manifest.tsv',sep='\t',index=False)
    summary.to_csv(out/'group_descriptive_summary.tsv',sep='\t',index=False)
    plt.rcParams.update({'font.size':9,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
    colors=['#a64e80','#4b9c85','#447da8','#caa044','#8165a3','#db7859','#78816d','#b6b8bb']
    fig,axes=plt.subplots(3,2,figsize=(13,11),sharey=True)
    for ax,(group,info) in zip(axes.flat,cfg['groups'].items()):
        f=table[table.group==group]
        pivot=f.pivot(index='source_label',columns='state',values='fraction')[states]
        bottom=np.zeros(len(pivot))
        for state,color in zip(states,colors):
            ax.bar(range(len(pivot)),pivot[state],bottom=bottom,color=color,label=state,width=.72)
            bottom+=pivot[state].to_numpy()
        labels=[label.split('_')[0]+'\nN='+str(manifest.set_index('source_label').loc[label,'retained_cells']) for label in pivot.index]
        ax.set_xticks(range(len(pivot)),labels);ax.set_ylim(0,1)
        ax.set_title(f'{info["driver"]} · {info["chase_days"]}-day chase · source group {group}',fontsize=10)
        ax.set_ylabel('Fraction of retained source cells')
    handles,labels=axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.94),ncol=4,frameon=False)
    fig.suptitle('HPCS tracing: source-level descendant RNA-state composition',fontsize=16,y=.98)
    fig.text(.5,.025,'GSE277777 · author trace assignments and cell-type labels · N = retained cells, not biological replicates\nSource aliases are not independently verified mice/pools. Current mScarlet status is unavailable; no p-values or transition rates.',ha='center')
    fig.tight_layout(rect=(0,.07,1,.875))
    for ext in ['png','svg']:fig.savefig(figdir/f'a1_hpcs_source_composition.{ext}',dpi=200,bbox_inches='tight')
    plt.close(fig)
    record=dict(status='completed_descriptive_source_reproduction',utc=datetime.now(timezone.utc).isoformat(),
        code_sha256=sha(Path(__file__)),input_sha256=source_hashes,selected_cells=len(data),source_labels=22,
        independent_biological_units_verified=False,current_mScarlet_available=False,
        saved_author_group_counts_verified=True,author_early_HPCS_fractions_verified=11,author_Hopx_means_verified=4,
        output_sha256={p.relative_to(BASE).as_posix():sha(p) for folder in [out,figdir] for p in sorted(folder.iterdir())})
    rp.write_text(json.dumps(record,indent=2)+'\n')
    print(summary[summary.state.isin(['HPCS','AT1-like'])].to_string(index=False))


if __name__=='__main__':main()
