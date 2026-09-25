"""Retrieve bounded public source metadata and author code; never execute it."""
from __future__ import annotations
import concurrent.futures
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from analysis.lib.provenance import archive_existing_record, write_json_atomic

PAPER = Path(__file__).resolve().parents[1]
CACHE = PAPER / 'cache' / 'u2_sources'
OUT = PAPER / 'trials' / 'u2_source_audit'
CAP = 12_000_000


def fetch(name, url):
    path = CACHE / name
    t = time.monotonic()
    try:
        if not path.exists():
            req = urllib.request.Request(url, headers={'User-Agent': 'Public-research-reanalysis/1.0'})
            with urllib.request.urlopen(req, timeout=40) as response:
                data = response.read(CAP + 1)
            if len(data) > CAP:
                raise ValueError('Source exceeds bounded metadata cap')
            path.write_bytes(data)
        data = path.read_bytes()
        if name.endswith('.xml'):
            root = ET.fromstring(data)
            if root.tag not in {'article', 'collection'}:
                raise ValueError('Response is not a full-text article or BioC collection')
        return {'name': name, 'url': url, 'bytes': len(data),
                'sha256': hashlib.sha256(data).hexdigest(), 'seconds': round(time.monotonic()-t, 3), 'status': 'retrieved'}
    except Exception as exc:
        return {'name': name, 'url': url, 'status': 'unavailable', 'error': str(exc)}


def main():
    CACHE.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    sources = [
        ('author_tree.json', 'https://api.github.com/repos/FuduanPeng/LungPCA_Code/git/trees/main?recursive=1'),
        ('zenodo_17148540.json', 'https://zenodo.org/api/records/17148540'),
        ('zenodo_17172149.json', 'https://zenodo.org/api/records/17172149'),
        ('peng_fulltext.xml', 'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12980502/fullTextXML'),
        ('peng_bioc.xml', 'https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_xml/PMC12980502/unicode'),
    ]
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        records = list(pool.map(lambda x: fetch(*x), sources))
    tree_path = CACHE / 'author_tree.json'
    if tree_path.exists():
        tree = json.loads(tree_path.read_text(encoding='utf-8'))
        sha = tree['sha']
        code_sources = [(x['path'], 'https://raw.githubusercontent.com/FuduanPeng/LungPCA_Code/'+sha+'/'+urllib.parse.quote(x['path']))
                        for x in tree['tree'] if x['type'] == 'blob' and x.get('size', CAP+1) < 1_000_000]
        with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
            records += list(pool.map(lambda x: fetch(*x), code_sources))
    inventories = []
    for record_id in ['17148540', '17172149']:
        path = CACHE / ('zenodo_'+record_id+'.json')
        if path.exists():
            obj = json.loads(path.read_text(encoding='utf-8'))
            inventories.append({'record': record_id, 'title': obj.get('metadata', {}).get('title'),
                                'access_right': obj.get('metadata', {}).get('access_right'),
                                'files': [{'key': f['key'], 'size': f['size'], 'checksum': f.get('checksum'), 'links': f.get('links')} for f in obj.get('files', [])]})
    result = {'run_utc': datetime.now(timezone.utc).isoformat(), 'sources': records,
              'zenodo_inventories': inventories, 'external_code_executed': False}
    archive_existing_record(OUT / 'run_record.json')
    write_json_atomic(OUT / 'run_record.json', result)
    print(json.dumps({'source_status': [(r['name'], r['status']) for r in records], 'zenodo_inventories': inventories}, indent=2), flush=True)


if __name__ == '__main__':
    main()
