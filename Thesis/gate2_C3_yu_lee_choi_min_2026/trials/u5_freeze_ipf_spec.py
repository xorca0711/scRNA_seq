"""Freeze the declared pathway sets and existing donor-level inputs before effects."""
from pathlib import Path
import csv
import json
import sys
from datetime import datetime, timezone
PAPER=Path(__file__).resolve().parents[1]
ROOT=PAPER.parents[1]
sys.path.insert(0,str(ROOT))
from analysis.lib.provenance import sha256_file, write_json_atomic


def main():
    out=PAPER/'trials/u5_ipf_spec'
    out.mkdir(exist_ok=True)
    if (PAPER/'trials/u5_ipf_pathways/run_record.json').exists():
        raise RuntimeError('Do not refreeze the specification after effects have started')
    inputs=[]; sets={}
    resources=[ROOT/'raw_data/msigdb/h.all.v2024.1.Hs.symbols.gmt',
               PAPER/'cache/u2_sources/c2.cp.reactome.v2024.1.Hs.symbols.gmt']
    for path in resources:
        inputs.append({'path':path.relative_to(ROOT).as_posix(),'sha256':sha256_file(path)})
        for line in path.read_text().splitlines():
            fields=line.split('\t'); sets[fields[0]]=sorted(set(fields[2:]))
    shared=['HALLMARK_TNFA_SIGNALING_VIA_NFKB','HALLMARK_INFLAMMATORY_RESPONSE','HALLMARK_IL6_JAK_STAT3_SIGNALING']
    mapping={
        'macrophages':shared+['HALLMARK_GLYCOLYSIS','HALLMARK_OXIDATIVE_PHOSPHORYLATION','HALLMARK_HYPOXIA','HALLMARK_INTERFERON_GAMMA_RESPONSE'],
        'fibroblasts':shared+['HALLMARK_TGF_BETA_SIGNALING','HALLMARK_WNT_BETA_CATENIN_SIGNALING','REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION'],
        'AT2':shared+['HALLMARK_HYPOXIA','HALLMARK_GLYCOLYSIS','HALLMARK_P53_PATHWAY','HALLMARK_E2F_TARGETS','HALLMARK_G2M_CHECKPOINT']}
    selected=sorted({name for values in mapping.values() for name in values})
    assert all(name in sets for name in selected)
    with (out/'pathway_genes.tsv').open('w',newline='') as handle:
        writer=csv.writer(handle,delimiter='\t'); writer.writerow(['set','gene'])
        writer.writerows((name,gene) for name in selected for gene in sets[name])
    with (out/'compartment_sets.tsv').open('w',newline='') as handle:
        writer=csv.writer(handle,delimiter='\t'); writer.writerow(['compartment','set'])
        writer.writerows((comp,name) for comp,names in mapping.items() for name in names)
    expected={}
    for cohort in ['GSE136831','GSE135893']:
        for suffix in ['_subtypes_counts.csv.gz','_subtypes_units.csv']:
            path=ROOT/'analysis/corrections/statistics/cache'/(cohort+suffix)
            inputs.append({'path':path.relative_to(ROOT).as_posix(),'sha256':sha256_file(path),'bytes':path.stat().st_size})
        prep=json.loads((ROOT/'analysis/corrections/statistics'/(cohort+'_subtypes_preparation.json')).read_text())
        expected[cohort]=prep['total_counts']
    spec={'frozen_utc':datetime.now(timezone.utc).isoformat(),'status':'released','cohorts':list(expected),
          'expected_total_counts':expected,'primary_floor':50,'sensitivity_floors':[30,100],
          'source':'audited subtype pseudobulks from full raw matrices; macrophages included',
          'input_sources':inputs,'sets':selected,'sets_per_compartment':mapping,
          'method':'TMM_filterByExpr_voom_official_CAMERA_estimated_inter_gene_correlation',
          'sensitivity_correlation':0.01,'multiple_testing':'global_BH_across_all_compartments_subtypes_per_cohort_floor_correlation_family',
          'annotation':'all_deposited_subtypes_in_cached_compartments_no_effect_based_selection',
          'covariates':'disease_only; age/sex unavailable in reused metadata; observational',
          'assayed_fraction_floor':0.7,'minimum_tested_set_members':10,'independent_validation':False,
          'coverage_limit':'Cache retained donors with >=50 broad-compartment cells. Floor-30 subtype sensitivity is conditional on that original eligible population; it cannot recover excluded broad compartments with 30-49 cells.',
          'resource_source':'https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/c2.cp.reactome.v2024.1.Hs.symbols.gmt'}
    write_json_atomic(out/'specification.json',spec)
    with (out/'expected_counts.tsv').open('w',newline='') as handle:
        writer=csv.writer(handle,delimiter='\t');writer.writerow(['cohort','total_counts']);writer.writerows(expected.items())
    print(f'Frozen {len(selected)} distinct pathways in two IPF cohorts; resources and inputs hashed.')


if __name__=='__main__': main()
