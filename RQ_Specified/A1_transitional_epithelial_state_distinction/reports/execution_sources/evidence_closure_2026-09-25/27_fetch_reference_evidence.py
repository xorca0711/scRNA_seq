"""Small primary references/catalog records for adaptive analysis eligibility."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import urllib.request
import time

BASE = Path(__file__).resolve().parents[1]


def main():
    report = BASE / 'reports/closure_reference_inventory.json'
    if report.exists():
        raise SystemExit('Refusing to overwrite reference inventory')
    items = []
    for name, pmc in [('cd44', 'PMC12484947'), ('ire1_2025', 'PMC12520674'),
                      ('sage_perturb', 'PMC13015491'), ('zhou', 'PMC8684104')]:
        items.append((name + '_bioc.xml',
                      f'https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_xml/{pmc}/unicode'))
    meta = json.loads((BASE / 'metadata/GSE141635.json').read_text())
    exps = [v.split('term=')[1] for s in meta['samples'] for v in s['relation'] if 'term=SRX' in v]
    assert len(exps) == len(set(exps)) == 20
    items.append(('pats_sra.xml', 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id='
                  + ','.join(exps) + '&rettype=full&retmode=xml'))
    records = []
    for name, url in items:
        path = BASE / 'cache/evidence_closure' / name
        record = dict(name=name, url=url, cap_bytes=20*1024**2)
        try:
            if path.exists():
                data = path.read_bytes()
                record['cache_reused'] = True
            else:
                with urllib.request.urlopen(url, timeout=90) as response:
                    data = response.read(record['cap_bytes']+1)
                assert len(data) <= record['cap_bytes']
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(data)
            record.update(path=path.relative_to(BASE).as_posix(), bytes=len(data),
                          sha256=hashlib.sha256(data).hexdigest(), status='retrieved')
        except Exception as exc:
            record.update(status='retrieval_failed', error=str(exc))
        records.append(record)
        print(name, record['status'], record.get('bytes', record.get('error')), flush=True)
        time.sleep(.4)
    report.write_text(json.dumps(dict(utc=datetime.now(timezone.utc).isoformat(), files=records,
                                     code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()), indent=2)+'\n')


if __name__ == '__main__':
    main()
