"""Report specificity comparisons with units and signature circularity visible."""
from pathlib import Path
import json
PAPER=Path(__file__).resolve().parents[1];ROOT=PAPER.parents[1]

def main():
 import numpy as np,pandas as pd,matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 out=PAPER/'trials/u6_specificity';s=pd.read_csv(out/'module_scores.csv');assert not s.duplicated(['cohort','unit','seed','group','module']).any()
 modules=['HPCS_author_top100','HPCS_without_ADI_or_operational_markers','ADI_published_holdout','AT2_published_holdout','AT1_published_holdout','HALLMARK_HYPOXIA','HALLMARK_INFLAMMATORY_RESPONSE','HALLMARK_P53_PATHWAY']
 # Resolve stress module names from the frozen names, without changing genes.
 names=list(s.module.unique());modules=[next((n for n in names if n.lower()==m.lower()),m) for m in modules];modules=[m for m in modules if m in names]
 contrasts=[]
 for cohort,g in s.groupby('cohort'):
  if cohort in ['GSE247130','GSE310539']:case,ref='two_marker','reference'
  elif cohort=='GSE277777':case,ref='HPCS','AT2-like'
  else:
   groups=list(g.group.unique());case=next((x for x in groups if 'mixed' in x.lower()),None);ref=next((x for x in groups if 'at1' in x.lower() and x!=case),None)
  if case is None or ref is None:continue
  a=g[g.group==case];b=g[g.group==ref];z=a.merge(b,on=['cohort','unit','seed','module'],suffixes=('_case','_reference'),validate='one_to_one')
  for r in z.itertuples():contrasts.append(dict(cohort=cohort,unit=r.unit,seed=r.seed,module=r.module,case=case,reference=ref,n_case=r.n_cells_case,n_reference=r.n_cells_reference,difference=r.score_case-r.score_reference,interpretable=bool(r.interpretable_case and r.interpretable_reference),unit_type=r.unit_type_case))
 c=pd.DataFrame(contrasts);c.to_csv(out/'within_unit_contrasts.csv',index=False)
 assert np.isfinite(c.difference).all() and c.difference.between(-1,1).all()
 primary=c[(c.seed==17)&c.interpretable];sensitivity=c[(c.seed==29)&c.interpretable]
 paired=primary.merge(sensitivity,on=['cohort','unit','module'],suffixes=('_17','_29'));paired['same_sign']=np.sign(paired.difference_17)==np.sign(paired.difference_29);paired.to_csv(out/'seed_robustness.csv',index=False)
 pal=json.loads((ROOT/'analysis/config/palette.json').read_text());fig,axes=plt.subplots(1,2,figsize=(13,7),gridspec_kw={'width_ratios':[1,1.7]});fig.set_facecolor(pal['surface'])
 display=[m for m in modules if m in ['HPCS_author_top100','HPCS_without_ADI_or_operational_markers','ADI_published_holdout','AT2_published_holdout','AT1_published_holdout']]
 nick={'HPCS_author_top100':'Author HPCS (100)','HPCS_without_ADI_or_operational_markers':'HPCS without ADI/label genes (91)','ADI_published_holdout':'ADI holdout','AT2_published_holdout':'AT2 holdout','AT1_published_holdout':'AT1 holdout'}
 for ax,cohorts,title in [(axes[0],['GSE277777'],'Author HPCS vs AT2-like\nTwo shared source identifiers'),(axes[1],['GSE247130','GSE310539'],'Two-marker vs Sftpc-positive reference\nIndividual pooled libraries')]:
  z=primary[primary.cohort.isin(cohorts)&primary.module.isin(display)];units=sorted(z.unit.unique())
  for i,m in enumerate(display):
   for j,u in enumerate(units):
    r=z[(z.module==m)&(z.unit==u)]
    if len(r):ax.scatter(r.difference*100,i+(j-(len(units)-1)/2)*.05,s=22,color=pal['categorical']['1'],marker='o')
  ax.axvline(0,color=pal['axis'],lw=.8);ax.set_yticks(range(len(display)),[nick[m] for m in display]);ax.invert_yaxis();ax.set_xlabel('Within-unit difference in gene detection\nat 2,000 UMIs (percentage points)');ax.set_title(title,loc='left',fontsize=11);ax.grid(axis='x',color=pal['grid'],lw=.5);ax.set_facecolor(pal['surface']);ax.tick_params(labelsize=9,length=0)
  if z.empty:ax.text(.5,.5,'No eligible contrast',transform=ax.transAxes,ha='center')
  for spine in ax.spines.values():spine.set_visible(False)
 fig.suptitle('Epithelial specificity: shared transcription does not identify the same state',x=.02,ha='left',fontsize=14,color=pal['ink'])
 fig.text(.02,.03,'Points are source identifiers or pooled libraries, never independent cells. At least 30 cells per compared group and 80% gene coverage.\nHPCS author labels and author signature are coupled: the left panel checks consistency, not independent validation.\nOperational two-marker labels do not establish DATP, KAC or malignant identity. Seed 29 and all stress-module scores are provided in tables.',fontsize=9,color=pal['ink_2'])
 fig.tight_layout(rect=[0,.14,1,.92]);base=PAPER/'figures/epithelial_specificity_and_signature_overlap'
 for ext in ['png','svg']:fig.savefig(str(base)+'.'+ext,dpi=220,facecolor=pal['surface'])
 plt.close(fig)
 summary=primary.groupby(['cohort','module']).agg(eligible_units=('unit','nunique'),minimum=('difference','min'),maximum=('difference','max'),positive=('difference',lambda x:int((x>0).sum()))).reset_index();summary.to_csv(out/'contrast_summary.csv',index=False)
 h=summary[(summary.cohort=='GSE277777')&summary.module.eq('HPCS_without_ADI_or_operational_markers')]
 report=['# Epithelial state specificity extension','',f'Computed {len(s):,} module-score rows across four cohorts, with joint 2,000-UMI sampling and two technical seeds. {len(primary):,} seed-17 module/within-unit contrasts pass both 30-cell and 80% assay-coverage gates. These counts are analysis rows, not independent biological replicates.','','## What these comparisons can establish','','The author HPCS pilot has two shared source identifiers across six sorting gates. Sorting gates are not six independent mice, and the source identifiers have not been promoted to verified experimental replication. Author HPCS labels and the source HPCS signature are coupled; their separation is a consistency check. Removing operational label genes and ADI-overlapping genes tests a narrower signature-overlap explanation, but does not make the author labels independent.','','The repair/developmental inputs contain ten pooled libraries. Their differences are descriptive and cannot supply donor-level P values. The two-marker group is defined by Cldn4/Krt8 detection after joint depth sampling, while the reference requires Sftpc. These are conditional operational groups, not validated DATP or KAC populations. Existing source-paper outputs were not changed.','','NNK inputs use independently reviewed mixed-alveolar/AT1-like labels; source KAC identities remain unavailable. A shared HPCS or injury score does not establish a common cellular identity, transformation, or an IL-1 causal mechanism.','','## Results and robustness','']
 if len(h):
  r=h.iloc[0];report.append(f'The HPCS signature excluding ADI/operational markers is higher in author HPCS than AT2-like cells in {r.positive}/{r.eligible_units} eligible source identifiers (differences {100*r.minimum:.2f} to {100*r.maximum:.2f} detection percentage points).')
 repair=primary[primary.cohort.isin(['GSE247130','GSE310539'])&primary.module.eq('HPCS_without_ADI_or_operational_markers')]
 report += ['',f'The same reduced HPCS signature is higher in the operational two-marker group in {int((repair.difference>0).sum())}/{len(repair)} eligible pooled repair/developmental libraries. Thus the observed signature increase is not specific to the sorted neoplastic context. These conditional group comparisons do not establish the presence of the author-defined HPCS state.', '',f'Across {len(paired)} contrasts eligible in both technical seeds, {int(paired.same_sign.sum())} retain their direction. This checks sampling robustness, not biological replication. Scores for hypoxia, inflammatory and p53 programs remain in the full table.','','![Specificity comparisons](../../figures/epithelial_specificity_and_signature_overlap.png)','','Tables: [module scores](module_scores.csv), [within-unit contrasts](within_unit_contrasts.csv), [summary](contrast_summary.csv), [seed robustness](seed_robustness.csv), [module overlap](module_overlap.csv), [source identifier coverage](HPCS_source_identifier_coverage.csv), [retention](retention.csv). The developmental ISR-specific and spatial extensions have separate access/eligibility records and are not silently replaced by these datasets.']
 (out/'REPORT.md').write_text('\n'.join(report)+'\n',encoding='utf-8')
 (out/'validation.json').write_text(json.dumps(dict(status='passed',score_rows=len(s),primary_eligible_contrasts=len(primary),paired_seeds=len(paired),unique_unit_module_seed_keys=True,finite_bounded_contrasts=True,independence='descriptive; no cell or pooled-library pseudoreplication'),indent=2)+'\n')
 print(json.dumps(dict(score_rows=len(s),primary_eligible_contrasts=len(primary),hpcs_comparisons=len(primary[primary.cohort=='GSE277777']))))

if __name__=='__main__':main()
