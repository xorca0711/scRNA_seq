"""Freeze available marker panels and complete orthogonal pathway modules."""
from pathlib import Path
import ast
import hashlib
import json
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    target = HERE / 'modules.json'
    if target.exists():
        raise SystemExit('modules.json already exists; do not silently refreeze')
    source = ROOT / 'Thesis/gate1_02_choi_2020/choi_2020_extracts.json'
    extract = json.loads(source.read_text(encoding='utf-8'))
    m = extract['marker_sets']
    entries = []
    def add(name, genes, kind, provenance, **kw):
        entries.append(dict(name=name, genes=genes, kind=kind,
                            provenance=provenance, **kw))
    prov = {'path':str(source.relative_to(ROOT)).replace('\\','/'),
            'sha256':digest(source),'doi':'10.1016/j.stem.2020.06.020'}
    for name, key in [('DATP_marker_panel','DATP'),('AT1_early_panel','AT1_early'),
                      ('AT1_late_panel','AT1_late'),('DATP_figure7_panel','DATP_genes_in_figure_7C')]:
        add(name,m[key],'paper_marker_panel',dict(prov,json_key='marker_sets.'+key),
            label_overlap=sorted(set(m[key]) & {'Cldn4','Krt8','Sftpc'}))
    add('DATP_label_free_panel',[g for g in m['DATP'] if g not in ['Cldn4','Krt8']],
        'paper_marker_panel_minus_label_genes',dict(prov,json_key='marker_sets.DATP'),
        label_overlap=[])
    src=ROOT/'Thesis/gate1_02_choi_2020/datp_epigenetics/trials/m1_closed_or_merely_silenced.py'
    sets={}
    for n in ast.parse(src.read_text(encoding='utf-8')).body:
        if isinstance(n,ast.Assign) and isinstance(n.targets[0],ast.Name) and n.targets[0].id in ['AT2_IDENTITY','TRANSITIONAL']:
            sets[n.targets[0].id]=ast.literal_eval(n.value)
    p={'path':str(src.relative_to(ROOT)).replace('\\','/'),'sha256':digest(src)}
    add('AT2_common_panel',[g for g in sets['AT2_IDENTITY'] if g != 'Cebpa'],
        'repository_marker_panel',dict(p,variable='AT2_IDENTITY',exclusion='Cebpa removed in every genotype'),label_overlap=[])
    add('DATP_PATS_holdout_panel',sets['TRANSITIONAL'],'repository_combined_marker_panel',
        dict(p,variable='TRANSITIONAL'),label_overlap=[])
    src=HERE/'sources/strunz_2020_supplementary_data_3.xlsx'
    table=pd.read_excel(src,sheet_name='cell_types_2')
    provenance={'path':str(src.relative_to(ROOT)).replace('\\','/'),'sha256':digest(src),
      'doi':'10.1038/s41467-020-17358-3','pmcid':'PMC7366678',
      'retrieved_from':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7366678/supplementaryFiles',
      'archive_member':'41467_2020_17358_MOESM6_ESM.xlsx','sheet':'cell_types_2',
      'selection':'Every reported gene for the named author cell type; no local data-dependent selection. Each source list is capped at 400 published markers, not a whole biological program.'}
    for celltype,name in [('Krt8 ADI','ADI_published_400'),('AT2 cells','AT2_published_400'),('AT1 cells','AT1_published_400')]:
        gs=table.loc[table.cluster==celltype,'gene'].tolist()
        assert len(gs)==400 and len(set(gs))==400
        add(name,gs,'complete_published_marker_list_400',dict(provenance,cluster=celltype),
            label_overlap=sorted(set(gs)&set(['Cldn4','Krt8','Sftpc'])))
        excluded=set(['Cldn4','Krt8','Sftpc','Cebpa'])
        add(name.replace('_400','_holdout'),[g for g in gs if g not in excluded],
            'complete_published_list_minus_labels_and_Cebpa',dict(provenance,cluster=celltype,excluded=sorted(excluded)),label_overlap=[])
    gmt=ROOT/'raw_data/msigdb/mh.all.v2024.1.Mm.symbols.gmt'
    wanted={'HALLMARK_P53_PATHWAY','HALLMARK_HYPOXIA','HALLMARK_INFLAMMATORY_RESPONSE'}
    for line in gmt.read_text(encoding='utf-8').splitlines():
        name,url,*genes=line.split('\t')
        if name in wanted:
            add(name,genes,'complete_MSigDB_Hallmark_module',
                {'path':str(gmt.relative_to(ROOT)).replace('\\','/'),'sha256':digest(gmt),'url':url,
                 'version':'2024.1.Mm'},label_overlap=sorted(set(genes)&{'Cldn4','Krt8','Sftpc'}))
    assert wanted <= {x['name'] for x in entries}
    target.write_text(json.dumps({'version':'ES1-2026-09-22','modules':entries,
      'unavailable_full_state_modules':{
       'DATP':'Local authoritative extract is a short marker panel; primary PMC7487779 supplement index lists PDFs only.',
       'PATS':'No complete state module cached; primary supplementary-table download unavailable in this session.',
       'developmental_AT2_maturation':'No independently sourced complete developmental maturation set found locally; AT2 identity is not a substitute.'}},indent=2)+'\n',encoding='utf-8')
    print('Frozen',len(entries),'panels/modules;',len(set(g for e in entries for g in e['genes'])),'unique genes')

if __name__ == '__main__':
    main()
