"""Candidate annotation marker and metadata review, before human effects."""
from pathlib import Path
import sys,json,re
PAPER=Path(__file__).resolve().parents[1];ROOT=PAPER.parents[1]

def main():
    import numpy as np,pandas as pd,anndata as ad
    out=PAPER/'trials/u5_human_pilot';cache=PAPER/'cache/u5_human_pilot'
    a=ad.read_h5ad(cache/'annotation_input.h5ad');labels=pd.read_csv(cache/'reference_candidate_annotations.csv.gz',index_col=0)
    assert labels.index.is_unique and set(labels.index)==set(a.obs_names)
    labels=labels.loc[a.obs_names];rows=[]
    genes='EPCAM KRT8 SFTPC SFTPA1 ABCA3 AGER HOPX SCGB1A1 FOXJ1 PTPRC LST1 TYROBP C1QA C1QB CD68 CSF1R PPARG FABP4 COL1A1 COL1A2 COL3A1 DCN LUM PDGFRA PI16 ACTA2 TAGLN MYH11 PECAM1 VWF MKI67 CLDN4 KRT17'.split();genes=[g for g in genes if g in a.var_names]
    x=a[:,genes].X.tocsr();total=a.obs.full_library_size.to_numpy()
    for threshold in [.2,.3]:
        confident=labels.ann_finest_level_uncertainty<=threshold
        for name in sorted(labels.ann_finest_level.unique()):
            idx=np.flatnonzero(confident & labels.ann_finest_level.eq(name));n=len(idx)
            if not n:continue
            detect=np.asarray((x[idx]>0).sum(axis=0)).ravel()/n
            means=np.asarray(x[idx].multiply((10000/np.maximum(total[idx],1))[:,None]).mean(axis=0)).ravel()
            for j,g in enumerate(genes):rows.append(dict(label=name,uncertainty_cutoff=threshold,cells=n,gene=g,detection_fraction=detect[j],mean_normalized_10000=means[j]))
    pd.DataFrame(rows).to_csv(out/'reference_marker_review.csv',index=False)
    coverage=[]
    for threshold in [.2,.3]:
        z=labels[labels.ann_finest_level_uncertainty<=threshold]
        for (patient,hist,label),n in z.groupby(['patient','histology','ann_finest_level'],observed=True).size().items():coverage.append(dict(patient=patient,histology=hist,label=label,uncertainty_cutoff=threshold,cells=n))
    pd.DataFrame(coverage).to_csv(out/'reference_confident_coverage.csv',index=False)
    # GEO titles and filenames are independent deposited fields; demographics
    # check title-derived patient identity, not pathology itself.
    source=json.loads((PAPER/'trials/u0_geo_design_audit/GSE308103.json').read_text());meta=[]
    for sample in source['samples']:
        title=sample['title'][0];patient=re.search(r'patient (\d+)',title)
        if not patient:continue
        hist='normal' if 'normal' in title.lower() else next((h for h in ['AAH','AIS','MIA','LUAD'] if h in title),None)
        chars=dict(x.split(': ',1) for x in sample['characteristics_ch1'] if ': ' in x)
        urls=[x for x in sample['supplementary_file'] if x!='NONE']
        assert len(urls)==1
        filename=urls[0].rsplit('/',1)[1];name=re.search(r'_P(\d+)_(Normal|AAH|AIS|MIA|LUAD)',filename)
        consistent=bool(name and int(name[1])==int(patient[1]) and name[2].lower()==hist.lower())
        meta.append(dict(gsm=sample['accession'],patient=int(patient[1]),histology=hist,title=title,filename=filename,title_filename_agree=consistent,age=chars.get('age'),sex=chars.get('Sex'),smoking=chars.get('smoking')))
    m=pd.DataFrame(meta);m.to_csv(out/'all_series_patient_histology_crosswalk.csv',index=False)
    consistency=[]
    for patient,g in m.groupby('patient'):
        consistency.append(dict(patient=int(patient),libraries=len(g),histologies=';'.join(sorted(g.histology.unique())),all_titles_filenames_agree=bool(g.title_filename_agree.all()),age_consistent=g.age.nunique(dropna=False)==1,sex_consistent=g.sex.nunique(dropna=False)==1,smoking_consistent=g.smoking.nunique(dropna=False)==1))
    pd.DataFrame(consistency).to_csv(out/'patient_metadata_consistency.csv',index=False)
    summary=dict(status='candidate_review_tables_ready',cells=a.n_obs,confident_primary=int((labels.ann_finest_level_uncertainty<=.2).sum()),full_series_libraries=len(m),full_series_patients=m.patient.nunique(),filename_title_disagreements=int((~m.title_filename_agree).sum()),pathology_limit='deposited histology labels; independent section-level ROI masks not supplied',decision='manual marker review still required before releasing final human subtype annotations')
    (out/'annotation_review_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary))

if __name__=='__main__':main()
