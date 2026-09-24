"""Validate paired human calculations and report reference-compatible niches."""
from pathlib import Path
import json
PAPER=Path(__file__).resolve().parents[1];ROOT=PAPER.parents[1]

def main():
 import numpy as np,pandas as pd,matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 from scipy.stats import false_discovery_control,spearmanr
 out=PAPER/'trials/u5_human_niche';cache=PAPER/'cache/u5_human_niche'
 pipeline=json.loads((out/'paired_pipeline_run_record.json').read_text())
 assert {'u5_human_aggregate.py','u5_human_paired_pathways.R','u5_human_compatibility.py'}<=set(pipeline['completed'])
 assert json.loads((out/'liana_run_record.json').read_text())['status']=='completed'
 annotation=json.loads((PAPER/'trials/u5_human_full/annotation_review_summary.json').read_text())
 assert annotation.get('release_decision_status')=='reviewed_for_reference_compatible_broad_niche_interpretation'
 pw=pd.read_csv(out/'camera.csv');families=0
 for _,g in pw.groupby(['config','cell_floor','correlation_setting']):assert np.allclose(false_discovery_control(g.PValue),g.FDR_global,rtol=1e-10,atol=1e-12);families+=1
 primary=pw[(pw.config=='unc20_pooled')&(pw.cell_floor==50)&pw.correlation_setting.eq('estimated')];fixed=pw[(pw.config=='unc20_pooled')&(pw.cell_floor==50)&pw.correlation_setting.eq('fixed001')]
 coverage=pd.read_csv(out/'paired_eligibility.csv');broad=coverage[(coverage.config=='unc20_pooled')&(coverage.cell_floor==50)&coverage.label.eq('__broad__')];broad.to_csv(out/'primary_broad_pair_coverage.csv',index=False)
 scores=pd.concat([pd.read_csv(f) for f in (cache/'liana').glob('*_scores.csv.gz')],ignore_index=True);assert scores.lr_probs.between(0,1).all();assert not scores.duplicated(['patient','histology','config','resource','source','target','ligand_complex','receptor_complex']).any()
 pairs=pd.read_csv(out/'liana_fixed_pairs.csv');calls=pd.read_csv(out/'liana_call_coverage.csv');rows=[];score_groups={k:g for k,g in scores.groupby(['patient','histology','resource','config'])};pair_groups={k:set(g[['source','target']].itertuples(index=False,name=None)) for k,g in pairs.groupby(['patient','histology','config'])};keys=['source','target','ligand_complex','receptor_complex']
 for call in calls[calls.config!='all_cells_primary'].itertuples():
  key=(call.patient,call.histology,call.resource);common_pairs=pair_groups.get((call.patient,call.histology,'all_cells_primary'),set())&pair_groups.get((call.patient,call.histology,call.config),set())
  row=dict(patient=call.patient,histology=call.histology,resource=call.resource,config=call.config,common_pairs=len(common_pairs),edge_jaccard=np.nan,shared_edge_spearman=np.nan,shared_edges=0)
  if common_pairs:
   def extract(config):
    g=score_groups.get((*key,config))
    if g is None:return pd.Series(dtype=float)
    g=g[[tuple(x) in common_pairs for x in g[['source','target']].itertuples(index=False,name=None)]];return g.set_index(keys).lr_probs
   a=extract('all_cells_primary');b=extract(call.config);union=set(a.index)|set(b.index);common=set(a.index)&set(b.index)
   row['edge_jaccard']=len(common)/len(union) if union else np.nan;row['shared_edges']=len(common)
   if len(common)>=3:
    ix=sorted(common);x=a.loc[ix].to_numpy();y=b.loc[ix].to_numpy()
    if np.unique(x).size>1 and np.unique(y).size>1:row['shared_edge_spearman']=spearmanr(x,y).statistic
  rows.append(row)
 robustness=pd.DataFrame(rows);robustness.to_csv(out/'liana_sensitivity_by_patient_histology.csv',index=False);robustness.groupby(['config','resource']).agg(groups=('patient','size'),comparable=('edge_jaccard','count'),median_edge_jaccard=('edge_jaccard','median'),median_shared_edge_spearman=('shared_edge_spearman','median')).reset_index().to_csv(out/'liana_sensitivity_summary.csv',index=False)
 # Every observed label/source profile is retained separately from edge scores.
 r=pd.read_csv(out/'primary_compatibility_contrasts.csv');v=pd.read_csv(out/'primary_compatibility_patient_values.csv');pal=json.loads((ROOT/'analysis/config/palette.json').read_text());contrast_order=[('AAH','normal'),('AIS','normal'),('MIA','normal'),('LUAD','normal'),('LUAD','AAH'),('LUAD','AIS'),('LUAD','MIA')]
 show=primary[primary.label=='__broad__'];row_order=sorted(set(zip(show.comp,show['set'])));fig,ax=plt.subplots(figsize=(12,10));fig.set_facecolor(pal['surface']);ax.set_facecolor(pal['surface'])
 for i,(comp,pathway) in enumerate(row_order):
  for j,(case,ref) in enumerate(contrast_order):
   z=show[(show.comp==comp)&(show['set']==pathway)&(show.case==case)&(show.reference==ref)]
   if z.empty:ax.scatter(j,i,s=32,marker='x',color=pal['muted'],linewidth=.7);continue
   assert len(z)==1;t=z.iloc[0];color=pal['categorical']['1' if t.Direction=='Up' else '2'];size=18+35*min(3,-np.log10(max(t.FDR_global,1e-300)));ax.scatter(j,i,s=size,facecolor=color if t.FDR_global<.05 else pal['surface'],edgecolor=color,lw=1)
 ax.set_xticks(range(len(contrast_order)),[a+' minus\n'+b for a,b in contrast_order]);ax.set_yticks(range(len(row_order)),[c+' | '+s.replace('HALLMARK_','').replace('REACTOME_','').replace('_',' ') for c,s in row_order]);ax.invert_yaxis();ax.tick_params(length=0,labelsize=8);ax.set_title('Human recipient pathways: paired patients, primary estimated correlation',loc='left',fontsize=13);ax.grid(color=pal['grid'],lw=.4);ax.set_axisbelow(True)
 for spine in ax.spines.values():spine.set_visible(False)
 fig.text(.02,.02,'Blue = up; orange = down. Filled circles: global q < 0.05; open circles: tested without that threshold; x: no eligible test.\nSize follows -log10(q), capped at 3. BH includes all eligible subtype/broad views and seven contrasts in the primary family.\nThe figure displays broad compartments. Source atlas labels identify reference-compatible populations and do not establish nonmalignant identity.',fontsize=9,color=pal['ink_2']);fig.tight_layout(rect=[0,.10,1,1])
 for ext in ['png','svg']:fig.savefig(PAPER/'figures'/('human_paired_recipient_pathways.'+ext),dpi=220,facecolor=pal['surface'])
 plt.close(fig)
 choices=[('IL1B','IL1R1_IL1RAP','macrophages','fibroblasts'),('TGFB1','TGFBR1_TGFBR2','macrophages','fibroblasts'),('CCL2','CCR2','fibroblasts','macrophages')];fig,axes=plt.subplots(1,3,figsize=(14,7));fig.set_facecolor(pal['surface']);plotted=[]
 for ax,(lig,rec,src,tgt) in zip(axes,choices):
  for i,(case,reference) in enumerate(contrast_order):
   a=v[(v.view=='broad')&(v.ligand==lig)&(v.receptor==rec)&(v.source==src)&(v.target==tgt)&(v.case==case)&(v.reference==reference)]
   if len(a):
    ax.scatter(a.paired_difference,i+np.linspace(-.12,.12,len(a)),s=15,color=pal['categorical']['1']);ax.scatter(a.paired_difference.mean(),i,marker='D',s=40,facecolor=pal['surface'],edgecolor=pal['ink'],lw=.8);plotted.extend(a.to_dict('records'))
   else:ax.text(.02,i,'No eligible contrast',transform=ax.get_yaxis_transform(),fontsize=7,color=pal['ink_2'])
  ax.set_yticks(range(len(contrast_order)),[a+' minus '+b for a,b in contrast_order]);ax.invert_yaxis();ax.axvline(0,color=pal['axis'],lw=.8);ax.set_title(f'{lig}/{rec}\n{src} → {tgt}',loc='left',fontsize=10);ax.set_xlabel('Paired RNA compatibility difference');ax.tick_params(labelsize=8,length=0);ax.grid(axis='x',color=pal['grid'],lw=.5);ax.set_facecolor(pal['surface'])
  for spine in ax.spines.values():spine.set_visible(False)
 fig.suptitle('Human niche RNA contrasts retain individual patients',x=.02,ha='left',fontsize=14);fig.text(.02,.02,'Dots = patient-level case-minus-reference values; diamond = mean. Every displayed contrast has at least three complete patients.\nCanonical examples are prespecified families; complete subtype/resource/sensitivity tables are retained. RNA scores do not measure communication.\nNative expression-eligibility results are separate; prior-count compatibility values can exist below native LR detection thresholds.',fontsize=9,color=pal['ink_2']);fig.tight_layout(rect=[0,.15,1,.91])
 for ext in ['png','svg']:fig.savefig(PAPER/'figures'/('human_paired_niche_RNA_contrasts.'+ext),dpi=220,facecolor=pal['surface'])
 plt.close(fig);pd.DataFrame(plotted).to_csv(out/'compatibility_figure_values.csv',index=False)
 record=dict(status='passed',primary_pathway_tests=len(primary),primary_pathway_q05=int(primary.FDR_global.lt(.05).sum()),fixed001_pathway_q05=int(fixed.FDR_global.lt(.05).sum()),BH_families_checked=families,primary_RNA_contrasts=len(r),native_LR_rows=len(scores),native_calls=int(calls.status.ne('no_eligible_pair').sum()),assay_identity_limit='reference-compatible atlas populations; cross-lineage RNA and excluded low-confidence cells remain material')
 (out/'validation.json').write_text(json.dumps(record,indent=2)+'\n')
 report=['# Human paired niche analysis','',f'All 75 libraries from 23 deposited patient labels were processed. {annotation["qc_cells"]:,} nuclei passed full-assay QC; {annotation["primary_confident_cells"]:,} met the primary reference-label confidence criterion. Repeated libraries were pooled within patient and histology. AAH, AIS, MIA and LUAD comparisons remain separate.','','## Analysis and interpretation','','Pathways use raw-count pseudobulks, filterByExpr, TMM, voom, a patient-blocked design, and official CAMERA. Global BH covers all eligible set/subtype/broad-view/histology contrasts within each declared analysis family. At least three complete patients are required. Broad and subtype views overlap and are not independent replications.','',f'The primary estimated-correlation family contains {len(primary)} tests, with {record["primary_pathway_q05"]} q < 0.05. The fixed-0.01 sensitivity has {record["fixed001_pathway_q05"]} q < 0.05; it does not replace the primary result.','',f'The descriptive RNA-compatibility arm contains {len(r)} primary resource/view contrasts, with paired patient values and omission ranges independently reconstructed. Native LIANA yielded {len(scores):,} score rows across resources and sensitivities. Pair coverage, native expression eligibility and missing edges are retained separately; absent edges are never scored as zero.','','Confidence 0.3, 30/100-cell floors, largest-library selection, prior counts 0.5/2, both LR resources, expression thresholds and cell caps remain labelled sensitivities. The full raw assay supplies library totals and normalization; the computational gene panel never replaces the assay denominator.','','## Identity and assay limitations','','These are reference-compatible candidate compartments, including AT2-like cells. The atlas does not establish nonmalignant status and does not supply a KAC or HPCS classifier. Low-confidence/unassigned cells remain in QC records and are excluded from subtype comparisons. Confidence retention by histology and full-cohort multi-marker support are reported separately in the full annotation review. Unassigned cells carry median IL1B count fractions of 52-72% across histologies; these subtype results cannot establish the dominant source across all recovered cells. Surfactant RNA is widespread outside epithelial assignments, so no cell type is assigned from SFTPC alone. Gene expression and native LR magnitude do not establish secretion, activation, physical contact or a feedback mechanism. Cross-sectional paired lesions do not prove a progression trajectory.','','The planned KAC/NF-kB joint association remains gated by the unavailable source-defined KAC measure. Human HPCS/ADI/ISR program comparisons are a separate descriptive specificity extension and are not substituted for that endpoint. Conditional ligand-target refits use their own target and expression gates.','','![Human paired pathways](../../figures/human_paired_recipient_pathways.png)','','![Human niche contrasts](../../figures/human_paired_niche_RNA_contrasts.png)','','[Pathways](camera.csv), [paired eligibility](paired_eligibility.csv), [primary RNA contrasts](primary_compatibility_contrasts.csv), [individual patients](primary_compatibility_patient_values.csv), [native-call coverage](liana_call_coverage.csv), [sensitivity summary](liana_sensitivity_summary.csv), [annotation review](../u5_human_full/ANNOTATION_REVIEW.md), [checks](validation.json).']
 (out/'REPORT.md').write_text('\n'.join(report)+'\n',encoding='utf-8');print(json.dumps(record))

if __name__=='__main__':main()
