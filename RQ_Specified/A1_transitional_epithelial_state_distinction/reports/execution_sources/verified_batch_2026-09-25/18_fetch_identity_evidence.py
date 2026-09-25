"""Fetch bounded public metadata/code for the outstanding identity audit.

Downloads are evidence only: no downloaded code is executed. The cached tree
pins author code to a commit. Existing inventory/output is never overwritten.
"""
from pathlib import Path
from datetime import datetime, timezone
from concurrent.futures import ThreadPoolExecutor
from urllib.parse import quote
import hashlib
import json
import argparse
import re
import urllib.request

BASE = Path(__file__).resolve().parents[1]
CACHE = BASE/'cache/followup_sources'


def fetch(item):
    name, url, expected_bytes = item
    p = CACHE/name
    try:
        if not p.exists():
            req = urllib.request.Request(url, headers={'User-Agent': 'scRNA_seq public research'})
            with urllib.request.urlopen(req, timeout=45) as r:
                data = r.read(20*1024**2+1)
            if len(data) > 20*1024**2:
                raise ValueError('Exceeds frozen 20-MiB source-document ceiling')
            if expected_bytes is not None and len(data) != expected_bytes:
                raise ValueError('Byte count differs from pinned Git tree')
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_bytes(data)
        data = p.read_bytes()
        return dict(path=p.relative_to(BASE).as_posix(),url=url,bytes=len(data),
                    sha256=hashlib.sha256(data).hexdigest())
    except Exception as exc:
        return dict(name=name,url=url,error=str(exc))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-id', default='identity_evidence_inventory')
    parser.add_argument('--metadata-only', action='store_true', help='Audit remaining sequence identities and assembly index only')
    parser.add_argument('--remaining-metadata', action='store_true', help='Compact submitted-file aliases and native chromosome sizes')
    parser.add_argument('--hpcs-sample-notebooks', action='store_true', help='Read source hash-to-mouse/pool assignments without executing code')
    args = parser.parse_args()
    if not re.fullmatch(r'[a-zA-Z0-9_-]+',args.run_id):
        parser.error('Invalid run ID')
    report = BASE/'reports'/f'{args.run_id}.json'
    if report.exists():
        raise SystemExit('Refusing to overwrite source inventory')
    tree = json.loads((CACHE/'hpcs_code_tree.json').read_text())
    chosen = [r for r in tree['tree'] if r['type'] == 'blob' and (
        r['path'] in ['README.md','tracing/README.md','tracing_depletion_analysis/README.md',
            'tracing_depletion_analysis/04_Evaluate_traced_cells_Slc4a11.ipynb',
            'tracing_depletion_analysis/05_Evaluate_traced_cells_Hopx.ipynb'] or
        (r['path'].startswith('tracing/') and (r['path'].endswith('README.md') or
            (r['path'].endswith('.ipynb') and r['size'] < 1000000))))]
    items = [('hpcs/'+r['path'],'https://raw.githubusercontent.com/dbetel/HPCS_LUAD/'+tree['sha']+'/'+quote(r['path']),r['size']) for r in chosen]
    items += [
        ('chm13v2.0.fa.gz.fai','https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz.fai',None),
        ('zenodo_17662770.json','https://zenodo.org/api/records/17662770',None),
    ]
    if args.metadata_only:
        items = [items[-2], ('GSE277777_directory.html','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE277nnn/GSE277777/suppl/',None)]
        for acc in ['GSE154966','GSE273123']:
            metadata = json.loads((BASE/'metadata'/f'{acc}.json').read_text())
            for sample in metadata['samples']:
                for relation in sample['relation']:
                    match = re.search(r'(SRX[0-9]+|SAMN[0-9]+)',relation)
                    if match:
                        accession = match.group()
                        items.append((f'ena_xml/{accession}.xml','https://www.ebi.ac.uk/ena/browser/api/xml/'+accession,None))
    if args.remaining_metadata:
        items = [
            ('hs1.chrom.sizes','https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.chrom.sizes',None),
            ('cd44_submitted_files.tsv','https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1140210&result=read_run&fields=run_accession,experiment_accession,sample_accession,sample_title,submitted_ftp&format=tsv',None),
        ]
    if args.hpcs_sample_notebooks:
        chosen = [r for r in tree['tree'] if r['type']=='blob' and r['path'].startswith('tracing/')
                  and r['path'].endswith('.ipynb') and r['size']>1000000]
        items = [('hpcs/'+r['path'],'https://raw.githubusercontent.com/dbetel/HPCS_LUAD/'+tree['sha']+'/'+quote(r['path']),r['size']) for r in chosen]
    with ThreadPoolExecutor(max_workers=3) as pool:
        records = list(pool.map(fetch, items))
    report.write_text(json.dumps(dict(utc=datetime.now(timezone.utc).isoformat(),
        author_code_commit=tree['sha'],code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        file_ceiling_bytes=20*1024**2,files=records),indent=2)+'\n')
    for r in records:
        print(r.get('path',r.get('name')),r.get('bytes',r.get('error')))


if __name__ == '__main__':
    main()
