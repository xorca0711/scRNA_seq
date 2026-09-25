"""Recover explicit CD44 identities; audit HPCS annotations and PATS resources."""
from pathlib import Path
from datetime import datetime, timezone
import ast
import csv
import gzip
import json
import re
import subprocess
import xml.etree.ElementTree as ET
import pandas as pd
from a1_robustness import sha
from a1_closure import cd44_key, annotation_audit

BASE = Path(__file__).resolve().parents[1]
ROOT = BASE.parents[1]


def main():
    out = BASE/'tables/evidence_closure'
    rp = BASE/'reports/closure_identity_run.json'
    if rp.exists() or (out.exists() and {p.name for p in out.iterdir()} != {'prior_numerical_sha256.json'}):
        raise SystemExit('Refusing to overwrite closure identity audit')
    contract_path = BASE/'config/closure_analysis_contract.json'
    cfg = json.loads(contract_path.read_text())
    assert cfg['status'] == 'frozen_before_new_numerical_results'
    out.mkdir(exist_ok=True)
    old_paths = subprocess.check_output(['git','ls-files',str(BASE.relative_to(ROOT)/'tables'),
                                        str(BASE.relative_to(ROOT)/'figures')],cwd=ROOT,text=True).splitlines()
    baseline = {p:sha(ROOT/p) for p in old_paths if (ROOT/p).is_file()}
    baseline_path = out/'prior_numerical_sha256.json'
    if baseline_path.exists():
        assert json.loads(baseline_path.read_text()) == baseline
    else:
        baseline_path.write_text(json.dumps(baseline,indent=2)+'\n')
    inputs = [contract_path, BASE/'reports/closure_source_inventory.json', BASE/'metadata/GSE273123.json']
    for item in json.loads(inputs[1].read_text())['files']:
        assert item['status'] == 'retrieved'
        assert sha(BASE/item['path']) == item['sha256']
        inputs.append(BASE/item['path'])
    meta = json.loads((BASE/'metadata/GSE273123.json').read_text())
    by_gsm = {s['accession']:s for s in meta['samples']}
    count = BASE/'cache/inputs/GSE273123/GSE273123_Count_Matrix.txt.gz'
    inputs.append(count)
    with gzip.open(count,'rt') as handle:
        columns = next(csv.reader(handle,delimiter='\t'))[1:]
    key_columns = {cd44_key(c,'column'):c for c in columns}
    assert len(key_columns) == len(columns) == 16
    rows, catalog = [], []
    for p in sorted((BASE/'cache/evidence_closure/sra').glob('*.xml')):
        root = ET.parse(p).getroot()
        exp = root.find('.//EXPERIMENT')
        gsm = root.find('.//SAMPLE').attrib['alias']
        attrs = {x.findtext('TAG'):x.findtext('VALUE') for x in root.findall('.//SAMPLE_ATTRIBUTE')}
        files = [f.attrib for f in root.findall('.//SRAFile') if f.attrib.get('supertype') == 'Original']
        catalog.append(dict(experiment=exp.attrib['accession'],gsm=gsm,library_name=root.findtext('.//LIBRARY_NAME',''),
                            attributes=json.dumps(attrs,sort_keys=True),original_files=json.dumps(files,sort_keys=True)))
        if gsm not in by_gsm:
            continue
        s = by_gsm[gsm]
        keys = {cd44_key(f['filename'],'fastq') for f in files}
        assert len(keys) == 1 and len(files) == 2
        key = next(iter(keys))
        assert key in key_columns
        title = s['title'][0]
        title_match = re.fullmatch(r'(WT|SPC Mut)([1-4]) AT2 CD44 (neg|pos)',title,re.I)
        assert title_match, title
        genotype = 'WT' if attrs['genotype']=='SP-C WT' else 'Mutant'
        assert attrs['genotype'] in ('SP-C WT','SP-C Mutant')
        assert (title_match[1]=='WT') == (genotype=='WT')
        assert key[1] == ('positive' if title_match[3]=='pos' else 'negative')
        rows.append(dict(sample_id=key_columns[key],source_alias=key[0],gsm=gsm,experiment=exp.attrib['accession'],
                         mouse=genotype+title_match[2],genotype=genotype,sort=key[1],GEO_title=title,
                         original_files=';'.join(f['filename'] for f in files),
                         evidence=f'https://www.ncbi.nlm.nih.gov/sra?term={exp.attrib["accession"]}'))
    frame = pd.DataFrame(rows).set_index('sample_id').loc[columns].reset_index()
    assert len(frame) == 16 and frame.gsm.nunique() == 16 and frame.mouse.nunique() == 8
    assert (frame.groupby('mouse')['sort'].nunique()==2).all()
    assert (frame.groupby('mouse').source_alias.nunique()==1).all()
    assert (frame.groupby('source_alias').mouse.nunique()==1).all()
    frame.to_csv(out/'cd44_sample_manifest.tsv',sep='\t',index=False)
    pd.DataFrame(catalog).to_csv(out/'SRA_explicit_identity_evidence.tsv',sep='\t',index=False)
    cells_path = BASE/'cache/evidence_closure/hpcs/tracing_depletion_analysis/01_Concatenate_cells.source.json'
    cells = json.loads(cells_path.read_text())
    author_map = ast.literal_eval(ast.parse(cells[142]['source']).body[0].value)
    assert author_map == cfg['hpcs']['mapping']
    obs_path = BASE/cfg['hpcs']['input']; inputs.extend([obs_path,cells_path])
    obs = pd.read_csv(obs_path,index_col=0,dtype=str,keep_default_na=False)
    old_cfg = json.loads((BASE/'config/hpcs_descendant_reconstruction.json').read_text())
    traced = obs[obs.batch.isin(old_cfg['select_batches'])]
    assert len(traced)==5333 and traced.index.is_unique and set(traced.Sorting_Groups)=={'traced'}
    audits = [dict(scope=name,**annotation_audit(data.to_dict('records'),author_map))
              for name,data in [('all_deposited',obs),('traced',traced)]]
    pd.DataFrame(audits).to_csv(out/'hpcs_annotation_audit.tsv',sep='\t',index=False)
    summaries = []
    for level,column in [('state','cell type'),('source','Classification'),('group','Group')]:
        for name,data in traced.groupby(column):
            row = dict(level=level,label=name,**annotation_audit(data.to_dict('records'),author_map))
            row['abstention_fraction'] = row['stringent_other']/row['cells']
            summaries.append(row)
    pd.DataFrame(summaries).to_csv(out/'hpcs_abstention_summary.tsv',sep='\t',index=False)
    pats_path = BASE/'cache/evidence_closure/pats_sra.xml'; inputs.append(pats_path)
    pats = []
    for p in ET.parse(pats_path).getroot().findall('.//EXPERIMENT_PACKAGE'):
        ex = p.find('.//EXPERIMENT'); runs=p.findall('.//RUN')
        files = [f.attrib for f in p.findall('.//SRAFile') if f.attrib.get('supertype')=='Original']
        pats.append(dict(experiment=ex.attrib['accession'],title=ex.findtext('TITLE'),
                         runs=';'.join(r.attrib['accession'] for r in runs),
                         spots=sum(int(r.attrib['total_spots']) for r in runs),
                         original_bytes=sum(int(f['size']) for f in files),
                         original_files=';'.join(f['filename'] for f in files)))
    assert len(pats)==20
    pd.DataFrame(pats).to_csv(out/'PATS_raw_resource_audit.tsv',sep='\t',index=False)
    record=dict(status='completed',utc=datetime.now(timezone.utc).isoformat(),
                cd44_libraries=16,cd44_mice=8,hpcs=audits,
                input_sha256={p.relative_to(BASE).as_posix():sha(p) for p in inputs},
                code_sha256={p.name:sha(p) for p in [Path(__file__),BASE/'scripts/a1_closure.py']},
                output_sha256={p.relative_to(BASE).as_posix():sha(p) for p in out.iterdir()})
    rp.write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(dict(cd44_libraries=16,cd44_mice=8,hpcs=audits),indent=2))


if __name__=='__main__':main()
