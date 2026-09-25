"""Freeze mouse niche resources after treatment-blind lineage review."""
from pathlib import Path
import sys,json
from datetime import datetime,timezone
PAPER=Path(__file__).resolve().parents[1];ROOT=PAPER.parents[1];sys.path.insert(0,str(ROOT))
from analysis.lib.provenance import sha256_file,write_json_atomic

def read_gmt(path):
    return {a[0]:list(dict.fromkeys(a[2:])) for a in (line.rstrip().split('\t') for line in path.open(encoding='utf-8'))}

def main():
    import pandas as pd
    out=PAPER/'trials/u4_resources';cache=PAPER/'cache/u4_resources'
    hcop=cache/'human_mouse_hcop_fifteen_column.txt.gz';h=pd.read_csv(hcop,sep='\t',dtype=str)
    h=h[h.support.str.split(',').str.len()>=3][['human_symbol','mouse_symbol']].drop_duplicates()
    h=h[~h.human_symbol.isin(['-','']) & ~h.mouse_symbol.isin(['-',''])].dropna()
    h=h[~h.human_symbol.duplicated(False)&~h.mouse_symbol.duplicated(False)]
    h.to_csv(out/'strict_one_to_one_orthologs.csv',index=False);mapping=dict(zip(h.human_symbol,h.mouse_symbol))
    res=pd.read_csv(PAPER/'trials/u5_full_source_panels/GSE136831/resource_edges.csv');rows=[];omitted=[]
    for r in res.itertuples(index=False):
        genes=set(r.ligand.split('_')+r.receptor.split('_'));missing=genes-set(mapping)
        if missing:omitted.append(dict(resource=r.resource,human_ligand=r.ligand,human_receptor=r.receptor,unmapped=';'.join(sorted(missing))));continue
        rows.append(dict(resource=r.resource,ligand='_'.join(mapping[g] for g in r.ligand.split('_')),receptor='_'.join(mapping[g] for g in r.receptor.split('_')),human_ligand=r.ligand,human_receptor=r.receptor,family='core' if r.ligand in ['IL1A','IL1B','AREG','HBEGF','TGFB1'] else 'exploratory'))
    pd.DataFrame(rows).to_csv(out/'mouse_resource_edges.csv',index=False);pd.DataFrame(omitted).to_csv(out/'unmapped_resource_edges.csv',index=False)
    mh=ROOT/'raw_data/msigdb/mh.all.v2024.1.Mm.symbols.gmt';mr=cache/'m2.cp.reactome.v2024.1.Mm.symbols.gmt'
    sets=read_gmt(mh)|read_gmt(mr)
    assignments=pd.read_csv(PAPER/'trials/u5_ipf_spec/compartment_sets.tsv',sep='\t')
    assignments.compartment=assignments.compartment.map({'macrophages':'myeloid','fibroblasts':'fibroblast','AT2':'alveolar'})
    assert set(assignments['set'])<=set(sets)
    assignments.to_csv(out/'compartment_sets.tsv',sep='\t',index=False)
    pd.DataFrame([dict(set=s,gene=g) for s in sorted(set(assignments['set'])) for g in sets[s]]).to_csv(out/'pathway_genes.tsv',sep='\t',index=False)
    spec=dict(frozen_utc=datetime.now(timezone.utc).isoformat(),scope='early treatment niche analysis with independent broad-lineage annotations; does not substitute for KAC phenotype',orthology='HCOP >=3 supporting resources; strictly unique human and mouse symbols; every complex subunit must map',pathways='same declared families, version 2024.1.Mm; broad myeloid/alveolar labels replace unsupported macrophage/AT2 specificity explicitly',inputs={str(p.relative_to(ROOT)):sha256_file(p) for p in [hcop,mh,mr,PAPER/'trials/u3_lineage_annotation/frozen_cluster_annotations.csv']},statistical_scope='early endpoint association and descriptive LR; later treatment-history unresolved; primary epithelial fraction has three IgG and two treated animals at >=100 alveolar cells',minimum_units_per_arm=3,niche_floors=[30,50,100],primary_floor=50,pathway_multiplicity='BH over all evaluable early-endpoint sets and subtype views; additional sensitivity reserves equally many unobserved late-endpoint tests at P=1; no claim to completed two-endpoint primary family')
    write_json_atomic(out/'specification.json',spec);print('Frozen',len(rows),'mouse LR rows;',len(set(assignments['set'])),'pathways;',len(h),'strict ortholog pairs')

if __name__=='__main__':main()
