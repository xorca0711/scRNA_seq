"""Recover the source ISR table and strictly map frozen mouse modules."""
from pathlib import Path
import sys,json,re,csv
from datetime import datetime,timezone
PAPER=Path(__file__).resolve().parents[1];ROOT=PAPER.parents[1];sys.path.insert(0,str(ROOT))
from analysis.lib.provenance import write_json_atomic,sha256_file

def main():
 import pymupdf,pandas as pd
 cache=PAPER/'cache/u2_remaining_context_access';pdf=cache/'41586_2023_6423_MOESM1_ESM.pdf';out=PAPER/'trials/u6_isr_extension';out.mkdir(exist_ok=True)
 doc=pymupdf.open(pdf);rows=[]
 for i,page in enumerate(doc):
  for group,ensembl,symbol in re.findall(r'(gene\.clic|extended)\s+(ENSMUSG\d+)\s+(\S+)',page.get_text()):rows.append(dict(source_group=group,ensembl_id=ensembl,mouse_symbol=symbol,pdf_page=i+1))
 table=pd.DataFrame(rows);assert len(table)>=100 and table.mouse_symbol.is_unique
 table.to_csv(out/'han_ISR_source_table.csv',index=False);conflict=table[table.ensembl_id.duplicated(keep=False)];conflict.to_csv(out/'source_duplicate_identifiers.csv',index=False)
 old=json.loads((PAPER/'trials/u6_specificity/module_specification.json').read_text());old_modules=old['modules'];hpcs=set(next(m['genes'] for m in old_modules if m['name']=='HPCS_author_top100'));adi=set(next(m['genes'] for m in old_modules if m['name']=='ADI_published_holdout'))
 genes=table.mouse_symbol.tolist();operational={'Cldn4','Krt8','Sftpc','Cebpa','Slc4a11'}
 modules=[dict(name='Han_ISR_source_symbols',genes=genes,kind='source_supplementary_table_full_symbol_list',label_overlap=sorted(set(genes)&operational)),dict(name='Han_ISR_without_HPCS_ADI_or_operational_markers',genes=[g for g in genes if g not in hpcs|adi|operational],kind='specificity_sensitivity_removed_shared_genes',label_overlap=[])]
 spec=dict(frozen_utc=datetime.now(timezone.utc).isoformat(),scope='additional source-defined ISR scoring in existing reusable contexts; original Han developmental expression comparator remains conditional',source='Han et al. Nature 2023 Supplementary Table 1',doi='10.1038/s41586-023-06423-8',pdf_sha256=sha256_file(pdf),symbol_policy='exact source symbols; source Ensembl IDs not substituted for symbols because the table duplicates one ID for different symbols; no guessed alias updates',source_duplicate_identifier_rows=conflict.to_dict('records'),modules=modules,original_HPCS_ADI_spec_sha256=sha256_file(PAPER/'trials/u6_specificity/module_specification.json'),technical_sampling='same 2000-UMI joint sampling algorithm and seeds, with this frozen gene panel; exact random realizations/groups may differ from the original larger panel')
 write_json_atomic(out/'module_specification.json',spec)
 # Preserve the already completed entrypoint; derive a separately named run.
 original=PAPER/'trials/u6_score_specificity.py';text=original.read_text(encoding='utf-8').replace('u6_specificity','u6_isr_extension').replace('New HPCS/repair specificity extension','Source ISR specificity extension')
 (PAPER/'trials/u6_score_isr.py').write_text(text,encoding='utf-8')
 orth=pd.read_csv(PAPER/'trials/u4_resources/strict_one_to_one_orthologs.csv');mapping=dict(zip(orth.mouse_symbol,orth.human_symbol));human=[];audit=[]
 for m in old_modules+modules:
  source=m['genes'];mapped=[]
  for g in source:
   h=mapping.get(g);audit.append(dict(module=m['name'],mouse_symbol=g,human_symbol=h or '',mapped=bool(h)))
   if h:mapped.append(h)
  assert len(mapped)==len(set(mapped))
  human.append(dict(name=m['name'],genes=mapped,source_genes=len(source),mapped_genes=len(mapped),ortholog_fraction=len(mapped)/len(source)))
 destination=PAPER/'trials/u6_human_specificity';destination.mkdir(exist_ok=True)
 pd.DataFrame(audit).to_csv(destination/'ortholog_mapping.csv',index=False)
 write_json_atomic(destination/'module_specification.json',dict(frozen_utc=spec['frozen_utc'],modules=human,mapping='frozen HCOP >=3 evidence support and strict 1:1 in both species; original source gene count retained as denominator',ortholog_sha256=sha256_file(PAPER/'trials/u4_resources/strict_one_to_one_orthologs.csv'),primary='descriptive patient-blocked AT2 program contrasts, no source-state identity assignment',minimum_assayed_source_fraction=.7))
 pd.DataFrame([dict(module=m['name'],gene=g,source_genes=m['source_genes']) for m in human for g in m['genes']]).to_csv(destination/'human_module_genes.tsv',sep='\t',index=False)
 print(json.dumps(dict(ISR_source_genes=len(genes),ISR_reduced_genes=len(modules[1]['genes']),source_duplicate_identifier_rows=len(conflict),human_modules=len(human))))

if __name__=='__main__':main()
