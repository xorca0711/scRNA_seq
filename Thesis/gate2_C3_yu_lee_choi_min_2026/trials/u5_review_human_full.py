"""Full-cohort label support and selection audit, separate from effects."""
from pathlib import Path
import json,sys
PAPER=Path(__file__).resolve().parents[1];ROOT=PAPER.parents[1];sys.path.insert(0,str(ROOT))
from analysis.lib.provenance import write_json_atomic,sha256_file,code_identity

def main():
 import numpy as np,pandas as pd
 out=PAPER/'trials/u5_human_full';niche=PAPER/'trials/u5_human_niche';cache=PAPER/'cache/u5_human_niche'
 assert json.loads((out/'processing_run_record.json').read_text())['status']=='completed_inputs_ready'
 assert json.loads((niche/'aggregation_run_record.json').read_text())['status']=='completed'
 q=pd.read_csv(out/'library_qc_and_annotation.csv');assert len(q)==75 and q.patient.nunique()==23 and q.gsm.is_unique
 a=q.groupby(['patient','histology']).agg(libraries=('gsm','nunique'),input_cells=('input_cells','sum'),qc_cells=('qc_cells','sum'),confident_primary=('confident_primary','sum'),confident_sensitivity=('confident_sensitivity','sum')).reset_index();a['QC_retained_fraction']=a.qc_cells/a.input_cells;a['primary_reference_fraction']=a.confident_primary/a.qc_cells;a['sensitivity_reference_fraction']=a.confident_sensitivity/a.qc_cells;a.to_csv(out/'patient_histology_annotation_retention.csv',index=False)
 markers='EPCAM KRT8 SFTPC SFTPA1 ABCA3 AGER HOPX SCGB1A1 FOXJ1 PTPRC LST1 TYROBP C1QA C1QB CD68 CSF1R PPARG FABP4 COL1A1 COL1A2 COL3A1 DCN LUM PDGFRA PI16 ACTA2 TAGLN MYH11 PECAM1 VWF MKI67 CLDN4 KRT17'.split();rows=[]
 for config in ['unc20_pooled','unc30_pooled']:
  m=pd.read_csv(niche/(config+'_units.csv'));z=np.load(cache/(config+'_all_labels.npz'));genes=z['genes'];det=z['detected_cells'];counts=z['counts'];assert len(m)==len(det);pos={g:i for i,g in enumerate(genes)}
  for key,ix in m.groupby(['histology','comp','label'],sort=True).groups.items():
   ids=list(ix);cells=int(m.loc[ids,'cells'].sum());total=int(m.loc[ids,'full_library_sum'].sum())
   for gene in markers:
    if gene not in pos:continue
    j=pos[gene];rows.append(dict(config=config,histology=key[0],comp=key[1],label=key[2],gene=gene,patients=m.loc[ids,'patient'].nunique(),cells=cells,detection_fraction=det[ids,j].sum()/cells,pseudobulk_CPM=counts[ids,j].sum()/total*1e6))
 pd.DataFrame(rows).to_csv(out/'full_cohort_marker_review.csv',index=False)
 frame=pd.DataFrame(rows);targets={'AT2':['EPCAM','SFTPC','ABCA3'],'fibroblasts':['COL1A1','DCN','LUM','PDGFRA'],'macrophages':['C1QA','C1QB','CD68','TYROBP','CSF1R']};facts=[]
 for (config,hist,comp),g in frame[frame.label=='__broad__'].groupby(['config','histology','comp']):
  wanted=g[g.gene.isin(targets[comp])];facts.append(dict(config=config,histology=hist,comp=comp,observed_markers=len(wanted),marker_detections=';'.join(f'{r.gene}:{r.detection_fraction:.3f}' for r in wanted.itertuples()),markers_detected_in_at_least_10_percent=int((wanted.detection_fraction>=.1).sum())))
 pd.DataFrame(facts).to_csv(out/'broad_lineage_marker_review.csv',index=False)
 metadata=json.loads((PAPER/'trials/u5_human_pilot/annotation_review_summary.json').read_text());assert metadata['filename_title_disagreements']==0
 report=dict(status='review_tables_ready',code=code_identity(ROOT,__file__),libraries=75,patients=23,qc_cells=int(q.qc_cells.sum()),primary_confident_cells=int(q.confident_primary.sum()),sensitivity_confident_cells=int(q.confident_sensitivity.sum()),assayed_reference_gene_fraction_min=float(q.model_gene_fraction.min()),assayed_reference_gene_fraction_max=float(q.model_gene_fraction.max()),marker_review_sha256=sha256_file(out/'full_cohort_marker_review.csv'),limitations=['healthy reference compatible candidate labels, not proof of nonmalignant identity','unassigned cells retained in QC records, excluded from subtype pseudobulks','fixed-RNA background means surfactant RNA alone cannot assign epithelial identity','fractions concern recovered assayed cells, not tissue abundance','histology labels come from deposition; no inferred normal-to-cancer lineage'],release_decision='Inspect broad_lineage_marker_review.csv and retention before interpreting effects; no automatic post-result relabelling')
 write_json_atomic(out/'annotation_review_summary.json',report);print(json.dumps({k:v for k,v in report.items() if k not in ['code','limitations','release_decision']}))

if __name__=='__main__':main()
