"""Report the separately frozen ISR extension without inventing Han inputs."""
from pathlib import Path
import json
PAPER=Path(__file__).resolve().parents[1];ROOT=PAPER.parents[1]

def main():
 import numpy as np,pandas as pd,matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 out=PAPER/'trials/u6_isr_extension';s=pd.read_csv(out/'module_scores.csv');contrasts=[]
 for cohort,g in s.groupby('cohort'):
  if cohort in ['GSE247130','GSE310539']:case,ref='two_marker','reference'
  elif cohort=='GSE277777':case,ref='HPCS','AT2-like'
  else:
   groups=list(g.group.unique());case=next((x for x in groups if 'mixed' in x.lower()),None);ref=next((x for x in groups if 'at1' in x.lower() and x!=case),None)
  if case is None or ref is None:continue
  z=g[g.group==case].merge(g[g.group==ref],on=['cohort','unit','seed','module'],suffixes=('_case','_reference'),validate='one_to_one')
  for r in z.itertuples():contrasts.append(dict(cohort=cohort,unit=r.unit,seed=r.seed,module=r.module,case=case,reference=ref,n_case=r.n_cells_case,n_reference=r.n_cells_reference,difference=r.score_case-r.score_reference,eligible=bool(r.interpretable_case and r.interpretable_reference)))
 c=pd.DataFrame(contrasts);c.to_csv(out/'within_unit_contrasts.csv',index=False);primary=c[(c.seed==17)&c.eligible];assert np.isfinite(c.difference).all()
 paired=primary.merge(c[(c.seed==29)&c.eligible],on=['cohort','unit','module'],suffixes=('_17','_29'));paired['same_direction']=np.sign(paired.difference_17)==np.sign(paired.difference_29);paired.to_csv(out/'seed_robustness.csv',index=False)
 summary=primary.groupby(['cohort','module']).agg(eligible_units=('unit','nunique'),positive=('difference',lambda v:int((v>0).sum())),minimum=('difference','min'),maximum=('difference','max')).reset_index();summary.to_csv(out/'contrast_summary.csv',index=False)
 pal=json.loads((ROOT/'analysis/config/palette.json').read_text());fig,ax=plt.subplots(figsize=(12,8));fig.set_facecolor(pal['surface']);ax.set_facecolor(pal['surface']);units=primary[['cohort','unit']].drop_duplicates().sort_values(['cohort','unit']);modules=list(s.module.unique())
 for j,m in enumerate(modules):
  for i,r in enumerate(units.itertuples(index=False)):
   z=primary[(primary.cohort==r.cohort)&(primary.unit==r.unit)&(primary.module==m)]
   if len(z):ax.scatter(z.difference*100,i+(j-.5)*.18,color=pal['categorical'][str(j+1)],s=32,marker=['o','s'][j],label=['Source ISR (129 genes)','ISR excluding HPCS/ADI/label genes (105)'][j] if i==0 else None)
 ax.axvline(0,color=pal['axis'],lw=.8);ax.set_yticks(range(len(units)),[f'{r.cohort} | {r.unit}' for r in units.itertuples(index=False)]);ax.invert_yaxis();ax.set_xlabel('Within-unit difference in gene detection at 2,000 UMIs (percentage points)');ax.set_title('ISR specificity: direction depends on the compared context and signature',loc='left',fontsize=13);ax.grid(axis='x',color=pal['grid'],lw=.5);ax.legend(loc='lower right',fontsize=8,frameon=False);ax.tick_params(labelsize=8,length=0)
 for sp in ax.spines.values():sp.set_visible(False)
 fig.text(.02,.02,'Repair/development: two-marker minus Sftpc reference; sorted cancer: author HPCS minus AT2-like; NNK: mixed alveolar minus AT1-like.\nEach row is a source identifier, pooled library or retained mouse library. This RNA signature does not measure ISR biochemical activity.\nOriginal Han developmental matrices were not found in the audited public processed-data sources; this is a signature extension in other datasets.',fontsize=9,color=pal['ink_2']);fig.tight_layout(rect=[0,.11,1,1])
 for ext in ['png','svg']:fig.savefig(PAPER/'figures'/('isr_specificity_contexts.'+ext),dpi=220,facecolor=pal['surface'])
 plt.close(fig)
 report=['# Source ISR signature extension','','Han et al. Supplementary Table 1 supplies 129 unique mouse symbols. The exact source-symbol set is retained; a 105-gene sensitivity removes HPCS, ADI and operational-label overlap. One Ensembl identifier occurs beside two different source symbols in the PDF. That discrepancy is recorded, and symbols were not silently replaced using those identifiers.','','This separate run scores existing reusable repair/developmental, sorted neoplastic and NNK inputs using the frozen joint-depth procedure. Its gene panel changes the exact seed realization, so its operational cell groups are not assumed identical to the original HPCS run. No animal-level inferential test is built from pooled wells or sorting gates.','',f'{len(primary)} primary module/unit contrasts are eligible; {int(paired.same_direction.sum())}/{len(paired)} retain their direction in seed 29. Directions and magnitudes are provided for each unit, including negative results. An ISR RNA score does not establish pathway activation or identify the mechanism of injury.','','## Conditional Han comparator','','The author repository was pinned to e1cbecd5dc0e36cc3b22f266804597b179096514. It contains code referencing local SHH1–SHH8 10x directories and local RDS objects, but no processed matrices in its file inventory. The paper lists raw-read BioProjects; its small Supplementary Data ZIP contains an ImageJ alveolar-thickness macro. Thus the original Han developmental perturbation comparison remains unavailable from the bounded processed-data sources audited here. Raw-read alignment was outside the authorized staged plan and was not started.','','![ISR context comparisons](../../figures/isr_specificity_contexts.png)','','[Source table](han_ISR_source_table.csv), [identifier discrepancy](source_duplicate_identifiers.csv), [within-unit values](within_unit_contrasts.csv), [summary](contrast_summary.csv), [seed checks](seed_robustness.csv), [module specification](module_specification.json). Primary source: [Han et al.](https://doi.org/10.1038/s41586-023-06423-8), [author code](https://github.com/MinhoLee-DGU/2023.Han.et.al.Nature).']
 (out/'REPORT.md').write_text('\n'.join(report)+'\n',encoding='utf-8');(out/'validation.json').write_text(json.dumps(dict(status='passed',score_rows=len(s),primary_contrasts=len(primary),seed_direction_agreement=int(paired.same_direction.sum()),seed_comparisons=len(paired),finite=True),indent=2)+'\n')
 print(json.dumps(dict(primary_contrasts=len(primary),seed_agreement=int(paired.same_direction.sum()),seed_comparisons=len(paired))))

if __name__=='__main__':main()
