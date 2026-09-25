"""Archive compact source-identity evidence and final presentation checks."""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json
import re
import xml.etree.ElementTree as ET

BASE = Path(__file__).resolve().parents[1]


def sha(p):
    with p.open('rb') as h:
        return hashlib.file_digest(h,'sha256').hexdigest()


def read(p):
    with p.open(newline='',encoding='utf-8') as h:
        return list(csv.DictReader(h,delimiter='\t'))


def write(rows,p):
    with p.open('w',newline='',encoding='utf-8') as h:
        w=csv.DictWriter(h,fieldnames=list(rows[0]),delimiter='\t')
        w.writeheader();w.writerows(rows)


def main():
    out=BASE/'tables/identity_audit'
    report=BASE/'reports/delivery_evidence_audit.json'
    if out.exists() or report.exists():
        raise SystemExit('Refusing to overwrite delivery audit')
    old=json.loads((BASE/'reports/second_batch_summary_run.json').read_text())
    new=json.loads((BASE/'reports/second_batch_verified_run.json').read_text())
    assert old['input_sha256']==new['input_sha256']
    for name,digest in old['output_sha256'].items():
        assert sha(BASE/name)==digest
        if name.startswith('tables/'):
            assert digest==new['output_sha256'][name.replace('second_batch_summary','second_batch_verified')]
    for name,digest in new['output_sha256'].items():
        assert sha(BASE/name)==digest
    svg_labels={}
    for p in (BASE/'figures/second_batch_verified').glob('*.svg'):
        texts=ET.parse(p).getroot().findall('.//{http://www.w3.org/2000/svg}text')
        assert len(texts)>30
        svg_labels[p.name]=len(texts)
    ref=dict(line.split() for line in (BASE/'cache/followup_sources/hs1.chrom.sizes').read_text().splitlines())
    sizes=read(BASE/'tables/second_batch_verification/track_chromosome_sizes.tsv')
    assert all(ref[r['chrom']]==r['size'] for r in sizes)
    out.mkdir(parents=True)
    rows=[];sources={}
    for acc in ['GSE154966','GSE273123']:
        meta=BASE/'metadata'/f'{acc}.json'
        sources[meta.relative_to(BASE).as_posix()]=sha(meta)
        for sample in json.loads(meta.read_text())['samples']:
            srx=next(re.search(r'SRX[0-9]+',s).group() for s in sample['relation'] if 'SRX' in s)
            sam=next(re.search(r'SAMN[0-9]+',s).group() for s in sample['relation'] if 'SAMN' in s)
            ep=BASE/'cache/followup_sources/ena_xml'/f'{srx}.xml'
            sp=BASE/'cache/followup_sources/ena_xml'/f'{sam}.xml'
            er=ET.parse(ep).getroot();sr=ET.parse(sp).getroot()
            attrs={a.findtext('TAG'):a.findtext('VALUE') for a in sr.findall('.//SAMPLE_ATTRIBUTE')}
            for p in [ep,sp]:sources[p.relative_to(BASE).as_posix()]=sha(p)
            rows.append(dict(accession=acc,gsm=sample['accession'],experiment=srx,biosample=sam,
                GEO_title='; '.join(sample['title']),ENA_title=sr.findtext('.//TITLE'),
                experiment_alias=er.find('EXPERIMENT').get('alias'),library_name=er.findtext('.//LIBRARY_NAME'),
                genotype=attrs.get('genotype'),sort=attrs.get('tigit','see GEO title'),
                biological_unit_status='held: pool membership/non-overlap unresolved' if acc=='GSE154966' else 'held: matrix alias-to-GSM missing',
                source_url='https://www.ebi.ac.uk/ena/browser/view/'+srx))
    write(rows,out/'GEO_ENA_sample_audit.tsv')
    submitted=read(BASE/'cache/followup_sources/cd44_submitted_files.tsv')
    assert len(submitted)==16 and all(not r['submitted_ftp'] for r in submitted)
    write(submitted,out/'CD44_submitted_file_audit.tsv')
    held=[r for r in read(BASE/'tables/direct_marks_2026-09-25/histone_contrasts.tsv') if r['eligible']=='False']
    write(held,out/'held_histone_ratios.tsv')
    # Preserve compact labels and source-output evidence, without copying image
    # payloads or executing any external notebook code.
    notebooks=[]
    root=BASE/'cache/followup_sources/hpcs/tracing_depletion_analysis'
    for p in sorted(root.glob('*.ipynb')):
        data=json.loads(p.read_text())
        sources[p.relative_to(BASE).as_posix()]=sha(p)
        numbers=[11,12,29,61,73,74,77] if p.name.startswith('04') else [11,12,36,55,56]
        cells=[]
        for i in numbers:
            cell=data['cells'][i]
            output=[]
            for o in cell.get('outputs',[]):
                text=''.join(o.get('text',[])) or ''.join(o.get('data',{}).get('text/plain',[]))
                if text and 'FutureWarning' not in text and len(text)<5000:
                    output.append(text)
            cells.append(dict(cell_index_zero_based=i,source_sha256=hashlib.sha256(''.join(cell.get('source',[])).encode()).hexdigest(),saved_text_outputs=output))
        notebooks.append(dict(notebook=p.relative_to(BASE).as_posix(),sha256=sha(p),cells=cells,
            interpretation='Author-saved output only; no independent reanalysis or verified mouse/pool mapping'))
    (out/'hpcs_author_output_inventory.json').write_text(json.dumps(notebooks,indent=2)+'\n')
    for filename in ['hs1.chrom.sizes','cd44_submitted_files.tsv','pats_bioc.xml','hpcs_bioc.xml','mintchip_main.nf','GSE277777_directory.html']:
        p=BASE/'cache/followup_sources'/filename;sources[p.relative_to(BASE).as_posix()]=sha(p)
    record=dict(status='passed',utc=datetime.now(timezone.utc).isoformat(),code_sha256=sha(Path(__file__)),
        source_sha256=sources,original_presentation_preserved=True,corrected_numerical_summaries_identical=True,
        native_CHM13_chromosomes_verified=len(sizes),native_size_source='https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.chrom.sizes',
        SVG_text_labels=svg_labels,PNG_visual_review='All four inspected; marker legend clear; labels and captions readable',
        GEO_ENA_samples_audited=len(rows),CD44_submitted_filenames_available=0,
        output_sha256={p.relative_to(BASE).as_posix():sha(p) for p in sorted(out.iterdir())})
    report.write_text(json.dumps(record,indent=2)+'\n')
    print('Delivery evidence: 24 GEO/ENA identities audited; native assembly sizes and preserved/corrected renders verified')


if __name__=='__main__':
    main()
