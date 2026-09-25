"""Validate and report the evaluable early mouse niche analyses and gates."""
from pathlib import Path
import json
from datetime import datetime,timezone
PAPER=Path(__file__).resolve().parents[1];ROOT=PAPER.parents[1]

def main():
    import numpy as np,pandas as pd,matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    from scipy.stats import false_discovery_control
    out=PAPER/'trials/u4_mouse_niche'
    for name in ['liana_run_record.json','compatibility_run_record.json']:
        assert json.loads((out/name).read_text())['status']=='completed'
    pw=pd.read_csv(out/'camera.csv');de=pd.read_csv(out/'DE_target_eligibility.csv')
    for _,g in pw.groupby(['cell_floor','correlation_setting']):
        q=false_discovery_control(g.PValue.to_numpy());assert np.allclose(q,g.FDR_global_early,rtol=1e-10,atol=1e-12)
        assert np.allclose(np.minimum(1,2*q),g.FDR_reserving_late_family,rtol=1e-10,atol=1e-12)
    c=pd.read_csv(out/'primary_compatibility_contrasts.csv');v=pd.read_csv(out/'primary_compatibility_donor_values.csv')
    keys=['resource','view','ligand','receptor','source','source_label','target','target_label']
    lookup={k:g for k,g in v.groupby(keys)}
    for row in c.itertuples(index=False):
        g=lookup[tuple(getattr(row,k) for k in keys)];ctl=g.treatment.eq('Control IgG')
        assert g.gsm.is_unique and min(ctl.sum(),(~ctl).sum())>=3
        delta=lambda z,col='score':z.loc[z.treatment!='Control IgG',col].mean()-z.loc[z.treatment=='Control IgG',col].mean()
        assert np.isclose(delta(g),row.difference_antiIL1B_minus_IgG,atol=1e-10)
        assert np.isclose(delta(g,'ligand_min'),row.ligand_component_difference,atol=1e-10)
        assert np.isclose(delta(g,'receptor_min'),row.receptor_component_difference,atol=1e-10)
        loo=[delta(g.drop(i)) for i in g.index];assert np.allclose([min(loo),max(loo)],[row.loo_min,row.loo_max],atol=1e-10)
    cov=pd.read_csv(PAPER/'trials/u3_lineage_annotation/sample_subtype_coverage.csv')
    broad=cov.groupby(['gsm','treatment','compartment']).cells.sum().unstack(fill_value=0)
    alveolar=[]
    for floor in [50,100,200]:
        for treatment,g in broad.groupby(level='treatment'):
            alveolar.append(dict(cell_floor=floor,treatment=treatment,eligible_animals=int((g.alveolar>=floor).sum()),total_animals=len(g)))
    pd.DataFrame(alveolar).to_csv(out/'alveolar_fraction_eligibility.csv',index=False)
    # Protect the narrative from reversing the treatment-arm counts.
    e100=pd.DataFrame(alveolar).query('cell_floor == 100').set_index('treatment')
    assert e100.loc['Control IgG','eligible_animals']==3
    assert e100.loc[e100.index!='Control IgG','eligible_animals'].item()==2
    palette=json.loads((ROOT/'analysis/config/palette.json').read_text());fig,axes=plt.subplots(1,2,figsize=(12,6),gridspec_kw={'width_ratios':[1,1.3]});fig.set_facecolor(palette['surface'])
    b=broad.reset_index().sort_values(['treatment','gsm']);numbers=b[['alveolar','fibroblast','myeloid']].to_numpy()
    cmap=LinearSegmentedColormap.from_list('coverage',[palette['surface'],palette['sequential_ramp'][2]])
    ax=axes[0];ax.imshow(np.log1p(numbers),cmap=cmap,vmin=0,vmax=np.log1p(numbers.max()),aspect='auto')
    for i in range(len(b)):
        for j in range(3):ax.text(j,i,str(numbers[i,j])+(' *' if j==0 and numbers[i,j]<100 else ''),ha='center',va='center',fontsize=10,color=palette['ink'])
    ax.set_yticks(range(len(b)),[r.gsm+' | '+('IgG' if r.treatment=='Control IgG' else 'anti-IL1β') for r in b.itertuples()]);ax.set_xticks(range(3),['Alveolar','Fibroblast','CSF1R+\nmyeloid']);ax.set_title('Observed cells per animal',loc='left',fontsize=11)
    ax=axes[1];planned=[('IL1B','IL1R1_IL1RAP','myeloid','fibroblast'),('IL1A','IL1R1_IL1RAP','myeloid','fibroblast'),('TGFB1','TGFBR1_TGFBR2','myeloid','fibroblast'),('IL1B','IL1R1_IL1RAP','myeloid','alveolar')];plotted=[]
    for i,(lig,rec,src,tgt) in enumerate(planned):
        z=c[(c.view=='broad')&(c.human_ligand==lig)&(c.human_receptor==rec)&(c.source==src)&(c.target==tgt)]
        if z.empty:ax.text(.02,i,'Ineligible/resource unavailable',transform=ax.get_yaxis_transform(),fontsize=8);continue
        assert len(z)==1;r=z.iloc[0];ax.plot([r.loo_min,r.loo_max],[i,i],color=palette['categorical']['1'],lw=2);ax.scatter(r.difference_antiIL1B_minus_IgG,i,s=30,color=palette['categorical']['1']);plotted.append(r.to_dict())
    ax.axvline(0,color=palette['axis'],lw=.8);ax.set_yticks(range(len(planned)),[f'{l}/{r}\n{s} → {t}' for l,r,s,t in planned]);ax.invert_yaxis();ax.set_ylim(len(planned)-.5,-.5);ax.set_xlabel('anti-IL1β minus IgG RNA compatibility\n(log2 CPM-based units)',fontsize=9);ax.set_title('Descriptive niche contrasts',loc='left',fontsize=11);ax.grid(axis='x',color=palette['grid'],lw=.5)
    for ax in axes:
        ax.tick_params(length=0,labelsize=8,colors=palette['ink']);ax.set_facecolor(palette['surface'])
        for sp in ax.spines.values():sp.set_visible(False)
    fig.suptitle('Early mouse niche analysis: eligible contrasts and primary phenotype limits',x=.02,ha='left',fontsize=14,color=palette['ink'])
    fig.text(.02,.02,'* Below the primary 100-alveolar-cell floor: only 3 IgG and 2 anti-IL1β animals remain. Author KAC annotations are also unavailable.\nCSF1R+ myeloid is not a pure macrophage population. Alveolar includes mixed AT2/AT1 and AT1-like cells, without a KAC/DATP assignment.\nRight: 50-cell niche floor; point = animal-mean difference; line = leave-one-animal-out range, not a confidence interval.\nRNA does not measure neutralization efficacy. Seven early libraries represent retained animals; later treatment history remains unresolved.',fontsize=8,color=palette['ink_2'])
    fig.subplots_adjust(left=.17,right=.98,top=.88,bottom=.25,wspace=1.05)
    for ext in ['png','svg']:fig.savefig(PAPER/'figures'/('mouse_early_niche_and_phenotype_gates.'+ext),dpi=220,facecolor=fig.get_facecolor())
    pd.DataFrame(plotted).to_csv(out/'figure_compatibility_values.csv',index=False)
    primary=pw[(pw.cell_floor==50)&(pw.correlation_setting=='estimated')];sens=pw[(pw.cell_floor==50)&(pw.correlation_setting=='fixed001')]
    lt_ineligible=bool((de[['q05_up','q05_down']]<10).all().all())
    qa=dict(status='passed',checked_utc=datetime.now(timezone.utc).isoformat(),primary_compatibility_rows_recomputed=len(c),pathway_BH_families_checked=6,primary_pathway_tests=len(primary),primary_pathway_q05=int((primary.FDR_global_early<.05).sum()),fixed001_sensitivity_q05=int((sens.FDR_global_early<.05).sum()),conditional_LT_ineligible=lt_ineligible,primary_KAC_fraction='ineligible: fewer than three treated animals at >=100 alveolar cells, and source KAC classifier unavailable')
    (out/'validation.json').write_text(json.dumps(qa,indent=2)+'\n')
    lines=['# Early mouse niche results and unreleased primary phenotype','', 'The evaluable early-endpoint niche analyses are complete using independent, treatment-blind broad-lineage review. They do not reconstruct the author KAC annotation and do not replace the planned KAC-fraction endpoint.','', '## Eligibility changes the question','', 'At the prespecified 100-alveolar-cell floor, only three IgG controls and two treated animals remain, below the three-per-arm requirement. The source-frozen KAC classifier is also unavailable. The 50-cell denominator sensitivity retains all seven animals, but it cannot repair the missing KAC definition. The later endpoint remains unresolved because treatment start/history is not mapped to libraries.','', 'Fibroblast clusters have consistent stromal markers. The broad CSF1R-positive myeloid group is mixed monocyte/inflammatory myeloid, not a macrophage-only population. The more specific macrophage-candidate cluster has only 19 cells overall and is ineligible for subtype inference. Smooth-muscle, ambiguous stromal and mixed-lineage clusters were not forced into fibroblast/macrophage labels.','', '## Completed niche outputs','', f'Native LIANA evaluated 96 sample/resource/configuration calls with independent subunit-expression checks. The early RNA-compatibility analysis produced {len(c)} eligible primary descriptive rows, with resource, cell-floor and prior-count sensitivities and animal-level component values. These are associations in independently reviewed compartments, not a completed two-endpoint causal analysis.','', f'Estimated-correlation CAMERA evaluated {len(primary)} primary early-endpoint set/view combinations; {qa["primary_pathway_q05"]} pass global early-family q < 0.05. The fixed-0.01 sensitivity has {qa["fixed001_sensitivity_q05"]} q < 0.05 and does not replace the primary result. A separate conservative adjustment reserves an equally large unobserved late-endpoint family; it is not a claim that the late analysis ran. Broad and subtype views overlap.', '', 'All eligible primary receiver DE fits have fewer than ten q < 0.05 targets in each direction, so the conditional mouse ligand-target arm is ineligible without lowering the threshold. No lack-of-effect or successful-neutralization claim follows from these RNA results.','', '![Early mouse niche and phenotype gates](../../figures/mouse_early_niche_and_phenotype_gates.png)','', 'Tables: [primary compatibility](primary_compatibility_contrasts.csv), [animal values](primary_compatibility_donor_values.csv), [all pathway tests](camera.csv), [niche eligibility](eligibility.csv), [phenotype denominator gate](alveolar_fraction_eligibility.csv), [DE target eligibility](DE_target_eligibility.csv), [checks](validation.json). Raw counts, annotations and resource hashes remain linked through the run records and U3/U4 specifications.']
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n');print(json.dumps(qa))

if __name__=='__main__':main()
