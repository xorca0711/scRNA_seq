"""Independently check donor-summary arithmetic and render fixed core examples."""
from pathlib import Path
import json,sys
from datetime import datetime,timezone
PAPER=Path(__file__).resolve().parents[1];ROOT=PAPER.parents[1]


def main():
    import pandas as pd,numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    out=PAPER/'trials/u5_ipf_compatibility'
    d=pd.read_csv(out/'primary_contrasts.csv')
    units=pd.read_csv(PAPER/'cache/u5_ipf_compatibility/donor_scores.csv.gz')
    units=units[(units.resource=='consensus')&(units.cell_floor==50)&(units.prior_count==1)&units.eligible_contrast]
    assert np.allclose(units.score,(units.ligand_min+units.receptor_min)/2)
    keys=['cohort','view','resource','family','ligand','receptor','source','source_label','target','target_label','cell_floor','prior_count']
    assert not d.duplicated(keys).any()
    grouped={key:part for key,part in units.groupby(keys,dropna=False)}
    for row in d.itertuples(index=False):
        part=grouped[tuple(getattr(row,k) for k in keys)]
        assert not part.donor.duplicated().any()
        x=part.loc[part.disease=='IPF','score'].to_numpy();y=part.loc[part.disease=='control','score'].to_numpy()
        assert (len(x),len(y))==(row.n_IPF,row.n_control)
        assert np.isclose(x.mean()-y.mean(),row.difference_IPF_minus_control)
        loo=np.r_[[(x.sum()-z)/(len(x)-1)-y.mean() for z in x],[x.mean()-(y.sum()-z)/(len(y)-1) for z in y]]
        assert np.isclose(loo.min(),row.loo_min) and np.isclose(loo.max(),row.loo_max)
    allfits=pd.read_csv(out/'contrasts.csv')
    base=allfits[(allfits.resource=='consensus')&(allfits.cell_floor==50)&(allfits.prior_count==1)]
    cross=base[base.view=='broad'].pivot(index=['family','ligand','receptor','source','target'],columns='cohort',values='difference_IPF_minus_control').dropna()
    cross['same_direction']=np.sign(cross.GSE136831)==np.sign(cross.GSE135893)
    cross.to_csv(out/'broad_cross_cohort_directions.csv')
    robustness=[]
    for key,part in allfits[allfits.resource=='consensus'].groupby(keys[:10]):
        primary=part[(part.cell_floor==50)&(part.prior_count==1)]
        if len(primary)!=1:continue
        effect=float(primary.difference_IPF_minus_control.iloc[0])
        robustness.append(dict(zip(keys[:10],key))|{'primary_effect':effect,'eligible_configurations':len(part),'sensitivity_min':part.difference_IPF_minus_control.min(),'sensitivity_max':part.difference_IPF_minus_control.max(),'all_available_same_sign':bool((np.sign(part.difference_IPF_minus_control)==np.sign(effect)).all())})
    pd.DataFrame(robustness).to_csv(out/'robustness.csv',index=False)
    palette=json.loads((ROOT/'analysis/config/palette.json').read_text())
    planned=[('IL1B','IL1R1_IL1RAP','macrophages','fibroblasts'),('IL1A','IL1R1_IL1RAP','macrophages','fibroblasts'),('TGFB1','TGFBR1_TGFBR2','macrophages','fibroblasts'),('AREG','EGFR','macrophages','fibroblasts'),('HBEGF','EGFR','macrophages','fibroblasts'),('IL1B','IL1R1_IL1RAP','macrophages','AT2'),('IL1A','IL1R1_IL1RAP','macrophages','AT2')]
    fig,axes=plt.subplots(1,2,figsize=(12,6.5),sharey=True,facecolor=palette['surface'])
    plotted=[]
    for ax,cohort in zip(axes,['GSE136831','GSE135893']):
        ax.set_facecolor(palette['surface']);ax.axvline(0,color=palette['axis'],lw=1)
        for i,(lig,rec,src,tgt) in enumerate(planned):
            sel=d[(d.cohort==cohort)&(d.view=='broad')&(d.ligand==lig)&(d.receptor==rec)&(d.source==src)&(d.target==tgt)]
            if sel.empty:
                ax.text(.02,i,'No eligible primary resource row',fontsize=8,color=palette['muted'],va='center',transform=ax.get_yaxis_transform());continue
            r=sel.iloc[0];color=palette['categorical']['1']
            ax.plot([r.loo_min,r.loo_max],[i,i],color=color,lw=2)
            ax.scatter(r.difference_IPF_minus_control,i,c=color,s=28,zorder=3)
            ax.text(.98,i+.18,f'n IPF/control {r.n_IPF}/{r.n_control}',transform=ax.get_yaxis_transform(),fontsize=7,ha='right',color=palette['ink_2'])
            plotted.append(r.to_dict())
        ax.set_title(cohort,fontsize=12,color=palette['ink']);ax.set_xlabel('IPF minus control RNA-compatibility score\n(log2 CPM-based units)',fontsize=9)
        ax.set_xlim(-2.3,2.3);ax.set_ylim(len(planned)-.5,-.5)
        ax.grid(axis='y',color=palette['grid'],lw=.5);ax.set_axisbelow(True)
        for sp in ax.spines.values():sp.set_visible(False)
        ax.tick_params(colors=palette['ink'],length=0,labelsize=8)
    axes[0].set_yticks(range(len(planned)),[f'{lig} / {rec}\n{src} to {tgt}' for lig,rec,src,tgt in planned],fontsize=8)
    fig.suptitle('IL-1 RNA compatibility differs between IPF cohorts',fontsize=15,x=.02,ha='left',color=palette['ink'])
    fig.text(.02,.04,'Points: donor-mean contrast. Lines: leave-one-donor-out range, NOT confidence intervals. Primary floor: 50 cells; prior.count = 1.\nSelected canonical examples were specified by ligand family; all eligible resource rows are in the tables. Descriptive, unadjusted observational contrasts.\nRNA compatibility is not measured signalling. Broad-compartment mixture and sparse controls remain limitations.',fontsize=8,color=palette['ink_2'])
    fig.subplots_adjust(left=.29,right=.98,top=.89,bottom=.19,wspace=.14)
    f=PAPER/'figures/ipu5_rna_compatibility'
    for ext in ['png','svg']:fig.savefig(f.with_suffix('.'+ext),dpi=220,facecolor=fig.get_facecolor())
    plt.close(fig);pd.DataFrame(plotted).to_csv(out/'figure_core_examples.csv',index=False)
    summary={'checked_utc':datetime.now(timezone.utc).isoformat(),'status':'passed','primary_rows_recomputed':len(d),'common_broad_rows':len(cross),'common_broad_same_direction':int(cross.same_direction.sum()),'scope':'independent donor-summary arithmetic and sensitivity table checks; not causal or statistical validation'}
    (out/'validation.json').write_text(json.dumps(summary,indent=2)+'\n')
    sizes=d.groupby(['cohort','view','family']).size()
    report=['# IPF donor-level RNA compatibility','', 'Completed in two previously inspected observational cohorts. These are descriptive contrasts; no disease-label permutation or FDR claim is made because an adequate exchangeability/covariate model was not established.','',
            '## Main result','', 'The canonical macrophage-to-fibroblast IL1B/IL1R1-IL1RAP contrast is near zero in GSE136831 and positive in GSE135893. Macrophage-to-AT2 results also differ. The data therefore do not establish a uniform increase in the proposed IL-1 circuit across IPF cohorts. Low fibroblast control counts, compartment mixture and unmeasured clinical covariates remain alternative explanations.','',
            '| Cohort | View | Family | Eligible primary contrasts |','|---|---|---|---:|']
    report += [f'| {c} | {v} | {fam} | {n} |' for (c,v,fam),n in sizes.items()]
    report += ['', '## Robustness and interpretation','',
               'The run evaluated both resources, 30/50/100-cell floors and prior.count 0.5/1/2. `robustness.csv` reports the range and sign consistency over available consensus configurations. These ranges are not confidence intervals. The 30-cell analysis is conditional on the older cache retaining donors with at least 50 broad-compartment cells. Broad and subtype views are not independent replications.', '',
               'All primary contrasts were independently reconstructed from individual donor scores, including the leave-one-out ranges. Ligand and receptor-component differences are retained so a combined score cannot conceal opposing changes. Missing assays or compartments were not zero-imputed; unsupported neutrophil/endothelial directions were not evaluated with the triad inputs.', '',
               '![Core RNA compatibility examples](../../figures/ipu5_rna_compatibility.png)','',
               'Sources: [frozen specification](specification.json), [all contrasts](contrasts.csv), [eligibility](eligibility.csv), [robustness](robustness.csv), [table checks](validation.json), [run record](run_record.json).', '',
               'The initial work-package label W1 was changed to the established U5 naming after execution. The original run record and source snapshots are preserved in `.history/original_entrypoints`; [relocation](path_relocation.json) records the change.']
    (out/'REPORT.md').write_text('\n'.join(report)+'\n')
    print(json.dumps(summary))


if __name__=='__main__':main()
