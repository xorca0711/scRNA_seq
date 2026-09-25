"""Within-IPF-donor state specificity; deposited states remain distinct."""
from pathlib import Path
import json
PAPER=Path(__file__).resolve().parents[1];ROOT=PAPER.parents[1]

def main():
 import numpy as np,pandas as pd,matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 out=PAPER/'trials/u6_ipf_specificity';s=pd.read_csv(out/'donor_state_scores.csv');assert not s.duplicated(['cohort','cell_floor','prior_count','donor','disease','label','module']).any()
 planned=[('GSE136831','Aberrant_Basaloid','ATII','primary'),('GSE135893','KRT5-/KRT17+','AT2','primary'),('GSE135893','Transitional AT2','AT2','secondary')];all_pairs=[];summaries=[]
 for cohort,case,reference,scope in planned:
  a=s[(s.cohort==cohort)&s.disease.eq('IPF')&s.eligible];z=a[a.label==case].merge(a[a.label==reference],on=['cohort','cell_floor','prior_count','donor','module'],suffixes=('_case','_reference'),validate='one_to_one');z['difference']=z.score_case-z.score_reference;z['case']=case;z['reference']=reference;z['scope']=scope;all_pairs.append(z)
  for key,g in z.groupby(['cell_floor','prior_count','module']):
   d=g.difference.to_numpy();loo=(d.sum()-d)/(len(d)-1) if len(d)>1 else np.array([np.nan]);summaries.append(dict(cohort=cohort,case=case,reference=reference,scope=scope,cell_floor=key[0],prior_count=key[1],module=key[2],n_paired_donors=len(g),eligible=len(g)>=3,mean_difference=d.mean(),positive_donors=int((d>0).sum()),loo_min=np.nanmin(loo) if len(d)>1 else np.nan,loo_max=np.nanmax(loo) if len(d)>1 else np.nan))
 pairs=pd.concat(all_pairs,ignore_index=True);r=pd.DataFrame(summaries)
 # Retain explicit empty/undercovered planned cells rather than dropping them.
 extra=[]
 for cohort,case,reference,scope in planned:
  for (floor,prior,module),g in s[s.cohort==cohort].groupby(['cell_floor','prior_count','module']):
   present=r[(r.cohort==cohort)&(r.case==case)&(r.cell_floor==floor)&(r.prior_count==prior)&(r.module==module)]
   if len(present):continue
   a=g[g.disease=='IPF'];n=len(set(a[a.label==case].donor)&set(a[a.label==reference].donor));coverage=bool(g.eligible.all())
   extra.append(dict(cohort=cohort,case=case,reference=reference,scope=scope,cell_floor=floor,prior_count=prior,module=module,n_paired_donors=n,eligible=False,mean_difference=np.nan,positive_donors=np.nan,loo_min=np.nan,loo_max=np.nan,reason='no_complete_same_donor_state_pairs' if coverage else 'insufficient_assayed_source_gene_fraction'))
 r=pd.concat([r,pd.DataFrame(extra)],ignore_index=True);r['reason']=r.reason.fillna('measured_with_pair_count_gate')
 pairs.to_csv(out/'within_donor_state_values.csv',index=False);r.to_csv(out/'state_contrast_summary.csv',index=False)
 # Also retain the narrower annotated-AT2 IPF/control comparison separately.
 disease=[]
 for cohort,label in [('GSE136831','ATII'),('GSE135893','AT2')]:
  a=s[(s.cohort==cohort)&s.label.eq(label)&s.eligible]
  for key,g in a.groupby(['cell_floor','prior_count','module']):
   x=g[g.disease=='IPF'];y=g[g.disease=='control']
   disease.append(dict(cohort=cohort,label=label,cell_floor=key[0],prior_count=key[1],module=key[2],n_IPF=x.donor.nunique(),n_control=y.donor.nunique(),eligible=min(x.donor.nunique(),y.donor.nunique())>=3,difference_IPF_minus_control=x.score.mean()-y.score.mean()))
 pd.DataFrame(disease).to_csv(out/'AT2_disease_descriptive_contrasts.csv',index=False)
 primary=r[(r.cell_floor==50)&(r.prior_count==1)];v=pairs[(pairs.cell_floor==50)&(pairs.prior_count==1)]
 for row in primary[primary.mean_difference.notna()].itertuples():
  a=v[(v.cohort==row.cohort)&(v.case==row.case)&(v.module==row.module)];assert a.donor.is_unique and len(a)==row.n_paired_donors;assert np.isclose(a.difference.mean(),row.mean_difference)
 modules=['HPCS_author_top100','HPCS_without_ADI_or_operational_markers','ADI_published_holdout','Han_ISR_source_symbols','Han_ISR_without_HPCS_ADI_or_operational_markers'];nick=['HPCS source','HPCS minus ADI/labels','ADI holdout','ISR source','ISR minus HPCS/ADI/labels'];pal=json.loads((ROOT/'analysis/config/palette.json').read_text());fig,axes=plt.subplots(1,3,figsize=(15,6));fig.set_facecolor(pal['surface'])
 for ax,(cohort,case,reference,scope) in zip(axes,planned):
  for i,module in enumerate(modules):
   a=v[(v.cohort==cohort)&(v.case==case)&v.module.eq(module)]
   ax.scatter(a.difference,i+np.linspace(-.13,.13,len(a)),color=pal['categorical']['1'],s=22)
   if len(a)>=3:ax.scatter(a.difference.mean(),i,s=48,marker='D',facecolor=pal['surface'],edgecolor=pal['ink'],lw=.8)
  n=v[(v.cohort==cohort)&(v.case==case)&v.module.eq(modules[0])].donor.nunique();ax.set_yticks(range(len(modules)),nick);ax.invert_yaxis();ax.axvline(0,color=pal['axis'],lw=.8);ax.set_xlabel('Within-donor mean log2 CPM difference');ax.set_title(f'{cohort} | {case}\nminus {reference} | n={n} | {scope}',loc='left',fontsize=10);ax.grid(axis='x',color=pal['grid'],lw=.5);ax.set_facecolor(pal['surface']);ax.tick_params(labelsize=8,length=0)
  if n==0:ax.text(.5,.5,'No complete donor pairs\nat the 50-cell floor',transform=ax.transAxes,ha='center',va='center',fontsize=10,bbox=dict(facecolor=pal['surface'],edgecolor='none'));ax.set_xticks([])
  for spine in ax.spines.values():spine.set_visible(False)
 fig.suptitle('Fibrosis specificity: compare deposited epithelial states within each donor',x=.02,ha='left',fontsize=14)
 fig.text(.02,.02,'Dots = donors with both source-labelled states; diamond = mean where at least three paired donors remain. Primary floor: 50 cells per state.\nStrict one-to-one orthologs and original source-list coverage are retained. Full-assay TMM normalization uses retained epithelial pseudobulks.\nAberrant basaloid, KRT5-/KRT17+ and transitional AT2 are separate deposited labels; none is relabelled KAC, HPCS or malignant.',fontsize=9,color=pal['ink_2']);fig.tight_layout(rect=[0,.16,1,.91])
 for ext in ['png','svg']:fig.savefig(PAPER/'figures'/('ipf_epithelial_state_specificity.'+ext),dpi=220,facecolor=pal['surface'])
 plt.close(fig)
 report=['# Fibrosis epithelial state specificity','','Both full IPF matrices were streamed into donor-by-deposited-state raw-count pseudobulks. Every epithelial unit total agrees with the earlier independent full-library stream. Source-defined HPCS, ADI, alveolar, stress and ISR programs use the frozen strict human/mouse ortholog map. Exact source-list denominators remain visible; missing genes are not interpreted as unexpressed.','','The primary state comparisons are Aberrant_Basaloid versus ATII in GSE136831 and KRT5-/KRT17+ versus AT2 in GSE135893. Transitional AT2 versus AT2 is a separate secondary comparison. These source labels describe different populations; similar scores are not evidence that they represent one shared state.','','Primary normalization uses all assayed genes in retained epithelial donor/state pseudobulks, with edgeR TMM and log2 CPM prior 1. The 30/100-cell and 0.5/2-prior analyses remain sensitivities. Comparisons retain paired IPF donors; no pooled-cell P values or cross-cohort equivalence test is used.','','## Reduced HPCS results','']
 for row in primary[primary.module=='HPCS_without_ADI_or_operational_markers'].itertuples():
  if pd.isna(row.mean_difference):report.append(f'- {row.cohort}, {row.case} versus {row.reference}: unevaluable at the primary setting ({row.reason}); no zero effect is imputed.')
  else:report.append(f'- {row.cohort}, {row.case} versus {row.reference}: {row.n_paired_donors} paired donors, mean {row.mean_difference:+.3f} log2 CPM units; {int(row.positive_donors)}/{row.n_paired_donors} individual differences positive. Declared pair-count gate: {"eligible" if row.eligible else "not met; individual values remain descriptive"}.')
 report += ['', 'These within-state RNA-program contrasts do not establish IL-1 dependence, fibrosis-to-cancer progression or malignant transformation. The annotated AT2 IPF/control contrast is provided separately and cannot stand in for an aberrant-state contrast.', '', '![Fibrosis state specificity](../../figures/ipf_epithelial_state_specificity.png)', '', '[Donor/state scores](donor_state_scores.csv), [paired values](within_donor_state_values.csv), [state summaries and sensitivities](state_contrast_summary.csv), [AT2 disease contrasts](AT2_disease_descriptive_contrasts.csv), [frozen labels and questions](specification.json).']
 (out/'REPORT.md').write_text('\n'.join(report)+'\n',encoding='utf-8');qa=dict(status='passed',primary_state_module_contrasts=len(primary),eligible_primary_state_module_contrasts=int(primary.eligible.sum()),all_paired_effects_reconstructed=True,no_cell_level_inference=True);(out/'validation.json').write_text(json.dumps(qa,indent=2)+'\n');print(json.dumps(qa))

if __name__=='__main__':main()
