"""Patient-level section summaries and measured maps; no inferred ROIs."""
from pathlib import Path
import json
PAPER=Path(__file__).resolve().parents[1];ROOT=PAPER.parents[1]

def main():
 import numpy as np,pandas as pd,matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 from matplotlib.colors import LinearSegmentedColormap
 out=PAPER/'trials/u5_spatial_context';cache=PAPER/'cache/u5_spatial_measurements'
 assert json.loads((out/'processing_run_record.json').read_text())['status']=='completed_measurements_require_report'
 q=pd.read_csv(out/'section_qc.csv');p=pd.read_csv(out/'whole_section_programs.csv');inventory=pd.read_csv(out/'deposited_member_inventory.csv');assert q.gsm.is_unique and len(q)==65
 candidate=inventory[inventory.name.str.contains('annot|cluster|region|roi',case=False,regex=True)];candidate.to_csv(out/'possible_pathology_annotation_members.csv',index=False);assert candidate.empty,'Review newly found annotation candidates before claiming region labels unavailable'
 assert q[q.series=='GSE307534'].coordinates_available.all() and q[q.series=='GSE307534'].position_version_disagreements.eq(0).all()
 assert not q[q.series!='GSE307534'].coordinates_available.any()
 v=p[p.eligible & p.series.eq('GSE307534')].copy();v['patient']=v.patient.astype(int)
 aggregated=v.groupby(['patient','histology','module']).agg(score=('mean_log1p_norm10000','mean'),sections=('gsm','nunique'),spots=('spots','sum')).reset_index();aggregated.to_csv(out/'patient_histology_programs.csv',index=False)
 weighted=[]
 for key,g in v.groupby(['patient','histology','module']):weighted.append(dict(zip(['patient','histology','module'],key))|dict(score=np.average(g.mean_log1p_norm10000,weights=g.spots)))
 weights=pd.DataFrame(weighted);pairs=[];summaries=[]
 for case,reference in [('LUAD','AAH'),('LUAD','AIS'),('LUAD','MIA')]:
  for mode,a in [('equal_sections',aggregated),('spot_weighted',weights)]:
   z=a[a.histology==case].merge(a[a.histology==reference],on=['patient','module'],suffixes=('_case','_reference'),validate='one_to_one');z['paired_difference']=z.score_case-z.score_reference;z['case']=case;z['reference']=reference;z['section_weighting']=mode;pairs.append(z)
   for module,g in z.groupby('module'):
    d=g.paired_difference.to_numpy();loo=(d.sum()-d)/(len(d)-1) if len(d)>1 else np.array([np.nan]);summaries.append(dict(case=case,reference=reference,module=module,section_weighting=mode,n_patients=len(g),eligible_primary_count=len(g)>=3,mean_paired_difference=d.mean(),loo_min=np.nanmin(loo) if len(d)>1 else np.nan,loo_max=np.nanmax(loo) if len(d)>1 else np.nan))
 paired=pd.concat(pairs,ignore_index=True);s=pd.DataFrame(summaries);paired.to_csv(out/'paired_patient_values.csv',index=False);s.to_csv(out/'paired_section_summary.csv',index=False)
 for row in s[s.eligible_primary_count].itertuples():
  a=paired[(paired.case==row.case)&(paired.reference==row.reference)&(paired.module==row.module)&(paired.section_weighting==row.section_weighting)];assert a.patient.is_unique;assert np.isclose(a.paired_difference.mean(),row.mean_paired_difference)
 pal=json.loads((ROOT/'analysis/config/palette.json').read_text());cmap=LinearSegmentedColormap.from_list('repo_spatial',pal['sequential_ramp']);fig,axes=plt.subplots(2,4,figsize=(13,8));fig.set_facecolor(pal['surface']);examples=[('GSM9226168','P1 | AAH'),('GSM9226169','P1 | LUAD')];genes=['IL1B','C1QA','COL1A1','KRT8'];maps={gsm:pd.read_csv(cache/gsm/'spatial_display_values.csv.gz') for gsm,_ in examples};limits=[]
 for j,gene in enumerate(genes):
  col='GENE_'+gene;all_values=np.concatenate([d[col].to_numpy() for d in maps.values()]);upper=float(np.quantile(all_values,.995));upper=upper if upper>0 else max(float(all_values.max()),1e-6);limits.append(dict(gene=gene,color_maximum=upper,scaling='shared across displayed sections; clipped at joint 99.5th percentile unless zero'))
  for i,(gsm,title) in enumerate(examples):
   a=maps[gsm];ax=axes[i,j];im=ax.scatter(a.pxl_col_in_fullres,a.pxl_row_in_fullres,c=a[col],s=3,cmap=cmap,vmin=0,vmax=upper,rasterized=True,linewidths=0);ax.set_aspect('equal');ax.invert_yaxis();ax.set_axis_off();ax.set_title(title+' | '+gene,fontsize=10);ax.set_facecolor(pal['surface'])
  bar=fig.add_axes([.045+j*.245,.16,.18,.012]);cb=fig.colorbar(im,cax=bar,orientation='horizontal');cb.set_label('log1p normalized RNA',fontsize=8);cb.ax.tick_params(labelsize=7,length=2)
 fig.suptitle('Measured spatial RNA in the first deposited patient pair',x=.02,ha='left',fontsize=14);fig.subplots_adjust(left=.02,right=.98,top=.90,bottom=.24,wspace=.18,hspace=.12)
 fig.text(.02,.02,'Every plotted point is a QC-passing tissue spot. Panels show whole sections; no gene-defined ROI or inferred cell contact is used.\nColors share a scale within each gene. The two sections are one patient, not independent replication.\nIL1B RNA, macrophage/stromal markers and KRT8 in mixed spots do not identify the producing or receiving cell.',fontsize=9,color=pal['ink_2'])
 for ext in ['png','svg']:fig.savefig(PAPER/'figures'/('human_spatial_measured_maps.'+ext),dpi=220,facecolor=pal['surface'])
 plt.close(fig);pd.DataFrame(limits).to_csv(out/'spatial_figure_color_limits.csv',index=False)
 shown=['GENE_IL1B','GENE_KRT8','HALLMARK_INFLAMMATORY_RESPONSE'];fig,ax=plt.subplots(figsize=(12,7));fig.set_facecolor(pal['surface']);ax.set_facecolor(pal['surface']);rows=[]
 for reference in ['AAH','AIS','MIA']:
  for module in shown:
   i=len(rows);a=paired[(paired.reference==reference)&(paired.module==module)&(paired.section_weighting=='equal_sections')];rows.append(f'LUAD minus {reference} | {module.removeprefix("GENE_").replace("HALLMARK_","")} | n={len(a)}')
   if len(a):
    jitter=np.linspace(-.14,.14,len(a));ax.scatter(a.paired_difference,i+jitter,s=14,color=pal['categorical']['1'],alpha=.8)
    if len(a)>=3:ax.scatter(a.paired_difference.mean(),i,s=52,marker='D',facecolor=pal['surface'],edgecolor=pal['ink'],linewidth=.9)
 ax.axvline(0,color=pal['axis'],lw=.8);ax.set_yticks(range(len(rows)),rows);ax.invert_yaxis();ax.set_xlabel('Within-patient difference in mean log1p normalized RNA');ax.set_title('Whole-section changes across all eligible deposited patient pairs',loc='left',fontsize=13);ax.grid(axis='x',color=pal['grid'],lw=.5);ax.tick_params(labelsize=8,length=0)
 for sp in ax.spines.values():sp.set_visible(False)
 fig.text(.02,.02,'Dots = patients; open diamond = mean where at least three patients are available. Repeated same-histology sections receive equal weight.\nSource lesion labels describe entire deposited samples. Cell mixture and independent-region annotation limits remain; no P values are assigned.\nSpot-weighted sensitivity and all planned pathway/gene values are in the companion tables.',fontsize=9,color=pal['ink_2']);fig.tight_layout(rect=[0,.13,1,1])
 for ext in ['png','svg']:fig.savefig(PAPER/'figures'/('human_spatial_paired_sections.'+ext),dpi=220,facecolor=pal['surface'])
 plt.close(fig)
 # Post-viral data remain individual whole-sample values without neighborhood tests.
 post=p[p.series.ne('GSE307534')].copy();post.to_csv(out/'postviral_individual_sample_programs.csv',index=False)
 post_summaries=[]
 for series,g in post.groupby('series'):
  for module,a in g[g.eligible].groupby('module'):
   case_mask=a.histology.str.startswith('COVID') if series=='GSE267226' else a.histology.str.contains('CD8')
   control_mask=a.histology.str.startswith('Control') if series=='GSE267226' else a.histology.str.contains('IgG')
   assert (case_mask|control_mask).all()
   x=a.loc[case_mask,'mean_log1p_norm10000'];y=a.loc[control_mask,'mean_log1p_norm10000']
   post_summaries.append(dict(series=series,module=module,case='PASC-PF' if series=='GSE267226' else 'anti-CD8',reference='control' if series=='GSE267226' else 'IgG',n_case=len(x),n_reference=len(y),case_mean=x.mean(),reference_mean=y.mean(),difference=x.mean()-y.mean(),inference='descriptive_only'))
 pd.DataFrame(post_summaries).to_csv(out/'postviral_descriptive_contrasts.csv',index=False)
 report=['# Spatial context: measured sections and regional inference limits','',f'Processed all {int(q.series.eq("GSE307534").sum())} deposited GSE307534 sections and all nine deposited post-viral count matrices. Human lesion sections represent {q[q.series=="GSE307534"].patient.nunique()} source patient labels, with repeated sections nested within patients. This is not additional replication independent of the RNA companion.','','## What is evaluable','','Both coordinate-file versions agree on their common barcodes in every human lesion section, and all filtered-matrix barcodes map to coordinates. QC uses the deposited tissue mask, at least 500 counts, 200 detected genes and at most 15% mitochondrial RNA. Whole-section programme means and individual paired patient differences are computed for LUAD versus AAH, AIS and MIA separately. Repeated sections are averaged equally, with spot-weighted sensitivity. These are descriptive sample-level comparisons, not an inferred progression trajectory.','','The archive inventory contains no verified independent spot-level pathology/ROI annotations. The planned pathological-region enrichment or within-region neighborhood analysis therefore remains unevaluable. Gene-score-derived regions were not substituted. Tissue spots contain mixtures, so co-expression does not assign a ligand source, receptor-bearing cell or physical contact.','','## Pathological post-viral context','','GSE267226 has three PASC-PF donors and two controls; GSE267228 has two IgG and two anti-CD8 mice. GEO archive inventories supply H5 count matrices and PNG images without coordinate tables; the H5 internal inventory also has no spatial coordinates. The author code loads complete local Space Ranger directories that are not included in these deposited files. Therefore individual whole-sample RNA-program values are available, but neighborhood or anatomical-region inference is not. The mouse intervention is anti-CD8, not anti-IL1B.','','![Measured spatial RNA](../../figures/human_spatial_measured_maps.png)','','![Whole-section paired changes](../../figures/human_spatial_paired_sections.png)','','[Section QC](section_qc.csv), [archive/H5 inventory](deposited_member_inventory.csv), [assay coverage](program_coverage.csv), [all whole-section values](whole_section_programs.csv), [patient paired values](paired_patient_values.csv), [paired summaries and weighting sensitivity](paired_section_summary.csv), [post-viral individual values](postviral_individual_sample_programs.csv). No threshold was relaxed to make a regional result pass.']
 report += ['', '## Descriptive post-viral findings', '', 'Whole-section IL1B RNA means are lower in the three deposited PASC-PF samples than in the two controls (0.041 versus 0.077 log1p-normalized units), and lower in anti-CD8 than IgG mice (0.352 versus 0.459). These small, unadjusted whole-section comparisons neither localize a pathogenic niche nor establish an IL-1 treatment response. In particular, a whole-section average is not a reproduction of a pathology-region-specific result. [All descriptive contrasts](postviral_descriptive_contrasts.csv) retain every measured gene/program and individual sample values.', '', 'Across RNA and spatial companions, all 23 RNA patient labels are present among the 25 spatial patient labels, with consistent deposited age/sex/smoking/ethnicity fields. This is metadata support for linkage, not genotype verification or proof that tissue pieces are identical. [Crosswalk check](RNA_spatial_crosswalk_validation.json). Histology-specific patient counts overlap and must not be summed as independent patients.']
 (out/'REPORT.md').write_text('\n'.join(report)+'\n',encoding='utf-8');write=dict(status='passed',sections=len(q),human_patient_labels=int(q[q.series=='GSE307534'].patient.nunique()),coordinate_versions_agree=True,unique_section_ids=True,paired_values_checked=True,independent_region_labels_available=False,postviral_coordinates_available=False);(out/'validation.json').write_text(json.dumps(write,indent=2)+'\n');print(json.dumps(write))

if __name__=='__main__':main()
