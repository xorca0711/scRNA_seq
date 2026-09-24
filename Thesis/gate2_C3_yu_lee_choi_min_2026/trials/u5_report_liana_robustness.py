"""Check donor/pair membership and report full-cell LR sensitivity and sources."""
from pathlib import Path
import json,sys
from datetime import datetime,timezone
PAPER=Path(__file__).resolve().parents[1];ROOT=PAPER.parents[1]

def main():
    import numpy as np,pandas as pd,matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from scipy.stats import spearmanr
    out=PAPER/'trials/u5_liana_robustness'
    palette=json.loads((ROOT/'analysis/config/palette.json').read_text())
    fig,axes=plt.subplots(2,2,figsize=(12,9),gridspec_kw={'width_ratios':[1,1.3]});fig.set_facecolor(palette['surface'])
    sensitivity=[];source_summary=[];fig_values=[];cohort_sizes=[];triads=[];checks=[]
    for ci,cohort in enumerate(['GSE136831','GSE135893']):
        here=out/cohort;state=json.loads((here/'run_record.json').read_text());assert state['status']=='completed'
        d=pd.read_csv(here/'donor_lr_scores.csv.gz');pairs=pd.read_csv(here/'fixed_pairs_by_config.csv');calls=pd.read_csv(here/'call_coverage.csv')
        coverage=pd.read_csv(here/'donor_subtype_coverage.csv')
        assert len(d)==state['rows'] and len(calls)==state['completed_calls']
        assert not d.duplicated(['config','resource','donor','source','target','ligand_complex','receptor_complex']).any()
        assert d.merge(pairs,on=['config','donor','disease','source','target'],how='left',indicator=True,validate='many_to_one')._merge.eq('both').all()
        spec=json.loads((here/'specification.json').read_text())
        for conf in spec['configs']:
            p=pairs[pairs.config==conf['name']]
            count={(r.donor,r.celltype):r.cells for r in coverage.itertuples()}
            assert all(count[(r.donor,r.source)]>=conf['floor'] and count[(r.donor,r.target)]>=conf['floor'] for r in p.itertuples())
        keys=['source','target','ligand_complex','receptor_complex']
        for (resource,donor),group in d.groupby(['resource','donor']):
            base=group[group.config=='all_cells_primary'].set_index(keys).lr_probs
            primary_called=((calls.donor==donor)&(calls.resource==resource)&(calls.config=='all_cells_primary')).any()
            for conf in spec['configs'][1:]:
                alternative=group[group.config==conf['name']].set_index(keys).lr_probs
                if not ((calls.donor==donor)&(calls.resource==resource)&(calls.config==conf['name'])).any():continue
                common=base.index.intersection(alternative.index);union=base.index.union(alternative.index)
                r=spearmanr(base.loc[common],alternative.loc[common]).statistic if len(common)>2 and base.loc[common].nunique()>1 and alternative.loc[common].nunique()>1 else np.nan
                sensitivity.append(dict(cohort=cohort,resource=resource,donor=donor,config=conf['name'],primary_pair_eligible=bool(primary_called),primary_edges=len(base) if primary_called else np.nan,alternative_edges=len(alternative),common_edges=len(common) if primary_called else np.nan,union_edges=len(union) if primary_called else np.nan,availability_jaccard=len(common)/len(union) if primary_called and len(union) else np.nan,common_edge_rank_spearman=r if primary_called else np.nan))
        # Joint triads are explicitly distinct from isolated pair availability.
        for floor in [30,50,100]:
            eligible=coverage[coverage.cells>=floor]
            for donor,g in eligible.groupby('donor'):
                triads.append(dict(cohort=cohort,donor=donor,disease=g.disease.iloc[0],floor=floor,complete_triad={'myeloid','fibroblast','epithelial'}<=set(g.compartment)))
        profiles=pd.read_csv(here/'all_label_source_recipient_profiles.csv.gz')
        assert profiles.detected_cells.le(profiles.cells).all() and profiles.detection_fraction.between(0,1).all()
        source=profiles[profiles.cells>=50]
        summary=source.groupby(['disease','label','gene']).agg(donors=('donor','nunique'),mean_of_donor_expression=('mean_normalized_10000','mean'),mean_of_donor_detection=('detection_fraction','mean')).reset_index();summary['cohort']=cohort;source_summary.append(summary)
        compact=source[source.gene.isin(['IL1B','IL1A','IL1R1','IL1RAP','IL1RN','IL1R2','SIGIRR','NLRP3','PYCARD','CASP1','GSDMD'])].copy();compact['cohort']=cohort
        compact.to_csv(here/'il1_context_donor_values.csv',index=False)
        primary_calls=calls[calls.config=='all_cells_primary']
        cohort_sizes.append(dict(cohort=cohort,source_annotated_donors=profiles.donor.nunique(),source_annotated_cells=int(profiles.drop_duplicates(['donor','label']).cells.sum()),primary_LR_donors=primary_calls.donor.nunique(),native_calls=len(calls),native_rows=len(d),empty_calls=int(calls.retained_edges.eq(0).sum())))
        ax=axes[ci,0];rows=pd.DataFrame(sensitivity);rows=rows[(rows.cohort==cohort)&(rows.resource=='consensus')]
        configs=[x['name'] for x in spec['configs'][1:]]
        for i,c in enumerate(configs):
            z=rows[rows.config==c];jitter=np.linspace(-.12,.12,len(z))
            ax.scatter(z.availability_jaccard,i+jitter,s=15,color=palette['categorical']['1'],alpha=.6)
            ax.plot([z.availability_jaccard.median()]*2,[i-.22,i+.22],color=palette['ink'],lw=2)
        ax.set_yticks(range(len(configs)),['500 cells / seed 1','500 cells / seed 2','30-cell floor','100-cell floor','5% detection','20% detection']);ax.set_xlim(-.03,1.03);ax.invert_yaxis()
        ax.set_xlabel('Edge-set overlap with full-cell primary (Jaccard)');ax.set_title(cohort+' | LR sensitivity',loc='left',fontsize=11)
        ax=axes[ci,1];s=source[(source.gene=='IL1B')&~source.label.str.contains('multiplet|doublet|unknown|unassigned',case=False,regex=True)];ranks=s.groupby('label').agg(mean=('mean_normalized_10000','mean'),n=('donor','nunique'));labels=ranks[ranks.n>=3].sort_values('mean',ascending=False).head(8).index.tolist()
        for i,label in enumerate(labels):
            for disease,offset,slot in [('control',-.14,'1'),('IPF',.14,'2')]:
                z=s[(s.label==label)&(s.disease==disease)].copy();z['cohort']=cohort;fig_values.append(z)
                ax.scatter(z.mean_normalized_10000,i+offset+np.linspace(-.045,.045,len(z)),s=17,color=palette['categorical'][slot],alpha=.75,label=disease if i==0 else None)
        ax.set_yticks(range(len(labels)),labels);ax.invert_yaxis();ax.set_xlabel('Mean IL1B counts per cell, normalized to 10,000');ax.set_title(cohort+' | observed IL1B sources',loc='left',fontsize=11);ax.legend(frameon=False,fontsize=8)
        checks.append(dict(cohort=cohort,rows=len(d),pair_membership='passed',cell_floors='passed',profiles='passed'))
    for ax in axes.ravel():
        ax.set_facecolor(palette['surface']);ax.grid(axis='x',color=palette['grid'],lw=.5);ax.set_axisbelow(True);ax.tick_params(labelsize=8,length=0,colors=palette['ink'])
        for spine in ax.spines.values():spine.set_visible(False)
    fig.suptitle('Full-cell LR robustness and the cellular sources of IL1B RNA',x=.02,ha='left',fontsize=15,color=palette['ink'])
    fig.text(.02,.015,'Each dot represents one donor. Left: consensus resource; black ticks show donor medians. Comparisons require eligible pairs in both settings.\nRight: eight labels with the highest pooled donor-mean IL1B expression among labels represented in ≥3 donors; each donor-label has ≥50 cells.\nMultiplet/unknown labels are excluded from source plots. Plots are descriptive, not confirmatory tests. Cell composition and inhibitory context have separate tables.',fontsize=8,color=palette['ink_2'])
    fig.subplots_adjust(left=.16,right=.98,top=.92,bottom=.13,wspace=.8,hspace=.35)
    for ext in ['png','svg']:fig.savefig(PAPER/'figures'/('ipf_full_cell_robustness_and_sources.'+ext),dpi=220,facecolor=fig.get_facecolor())
    pd.DataFrame(sensitivity).to_csv(out/'donor_sensitivity.csv',index=False)
    pd.concat(source_summary).to_csv(out/'all_label_context_summary.csv',index=False)
    pd.concat(fig_values).to_csv(out/'figure_IL1B_source_donor_values.csv',index=False)
    pd.DataFrame(triads).to_csv(out/'complete_triad_coverage.csv',index=False)
    pd.DataFrame(cohort_sizes).to_csv(out/'cohort_completion.csv',index=False)
    s=pd.DataFrame(sensitivity).groupby(['cohort','resource','config']).agg(donors=('donor','nunique'),median_edge_overlap=('availability_jaccard','median'),median_common_edge_rank_r=('common_edge_rank_spearman','median')).reset_index();s.to_csv(out/'sensitivity_summary.csv',index=False)
    tri=pd.DataFrame(triads);ntri=tri[(tri.floor==50)&tri.complete_triad].groupby(['cohort','disease']).size()
    qa=dict(status='passed',checked_utc=datetime.now(timezone.utc).isoformat(),checks=checks);(out/'validation.json').write_text(json.dumps(qa,indent=2)+'\n')
    lines=['# Full-cell IPF LR robustness and source context','', 'Both cohorts were evaluated with native LIANA CellChat-like magnitude. All-cell inputs retain each cell’s original full-library total and full-assay maximum, preserving normalization and native magnitude scaling. Two-donor/two-resource comparisons against full-gene inputs had zero score differences. The computational maximum column is not a gene and never appears in an interaction.','', '| Cohort | Annotated IPF/control donors | Annotated cells | Primary LR donors | Native calls | Empty calls |','|---|---:|---:|---:|---:|---:|']
    lines += [f"| {r['cohort']} | {r['source_annotated_donors']} | {r['source_annotated_cells']} | {r['primary_LR_donors']} | {r['native_calls']} | {r['empty_calls']} |" for r in cohort_sizes]
    lines += ['', '## Interpretation','', 'The full-cell, two-seed 500-cell cap, 30/50/100-cell floor, 5/10/20% detection and two-resource outputs are all retained. Sensitivities vary one choice at a time. Availability overlap and rank correlation are different checks; correlation is calculated only on common edges and cannot conceal absent edges. No missing edge is imputed as zero. The cap seeds use stable donor/subtype seeds and are not an exact replay of the original pilot sampling.','', 'All annotated compartments were scanned for source, receptor, inhibitor and processing-gene RNA. The tables distinguish within-label expression/detection from proportions of recovered cells. Neither quantity measures mature cytokine secretion, receptor activation or absolute tissue abundance. WNT-family resource edges are explicitly exploratory.','', 'Native LIANA 1.10.0 raises an indexing error when no edge passes expression eligibility. In the second cohort these cases were independently checked across every resource subunit and fixed pair, then recorded as empty calls. Completed cached calls were retained; the original failed run and code are preserved under `.history`.','', 'Complete ≥50-cell triads, distinct from pair coverage: '+', '.join(f'{c} {d}: {n}' for (c,d),n in ntri.items())+'. Joint associations require at least ten complete patients, so this coverage must be checked for each proposed association.','', '![Full-cell robustness and sources](../../figures/ipf_full_cell_robustness_and_sources.png)','', 'Tables: [sensitivity summary](sensitivity_summary.csv), [donor sensitivity](donor_sensitivity.csv), [source/context summary](all_label_context_summary.csv), [triad coverage](complete_triad_coverage.csv), [independent checks](validation.json). Per-cohort folders contain donor-level IL-1 context, cell composition, native score caches and exact resource definitions. This package does not complete the human-lesion, mouse-perturbation or ligand-target arms.']
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n');print(json.dumps(qa))

if __name__=='__main__':main()
