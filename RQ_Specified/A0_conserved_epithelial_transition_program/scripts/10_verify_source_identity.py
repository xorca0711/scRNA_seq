"""Check the author donor/capture crosswalk against GEO, without expression effects."""
import csv
import gzip
import hashlib
import json
import re
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
CACHE = BASE / 'cache/continuation_v1'
OUT = BASE / 'tables/pilot_v1'


def main():
    meta_path = CACHE / 'Sountoulidis_lungdev_meta.tsv'
    geo_path = CACHE / 'GSE215898_family.soft.gz'
    with meta_path.open() as f:
        captures = sorted({(r['donor'], r['age'], r['name'])
                           for r in csv.DictReader(f, delimiter='\t')})
    with gzip.open(geo_path, 'rt') as f:
        soft = f.read()
    geo = {}
    for block in soft.split('^SAMPLE = ')[1:]:
        title = re.search(r'^!Sample_title = (.+)', block, re.M)
        age = re.search(r'^!Sample_characteristics_ch1 = age: PCW([0-9.]+)', block, re.M)
        if title and age:
            name = re.search(r'\[(10[Xx]\d+w\d+)\]', title[1])
            if name:
                key = name[1].lower().replace('w', '_')
                assert key not in geo
                geo[key] = (block.splitlines()[0], float(age[1]))
    assert len(captures) == len(geo) == 39
    assert len({c[0] for c in captures}) == 17
    assert len({c[2] for c in captures}) == len(captures)
    rows = []
    for donor, age, capture in captures:
        accession, geo_age = geo[capture.lower()]
        assert float(age) == geo_age
        rows.append(dict(donor=donor, capture=capture, age_pcw=age,
                         GEO_sample=accession, GEO_age_pcw=geo_age))
    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / 'D2_source_crosswalk.tsv').open('x', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    record = dict(status='PASS', captures=39, author_donors=17, age_conflicts=0,
                  interpretation='Author donor IDs define units; GEO captures do not add replication. Matching ages check source identity, not donor independence by themselves.',
                  inputs={p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in [meta_path, geo_path]})
    with (OUT / 'source_identity.json').open('x') as f:
        json.dump(record, f, indent=2)
        f.write('\n')
    print(json.dumps({k: record[k] for k in ['status', 'captures', 'author_donors', 'age_conflicts']}))


if __name__ == '__main__':
    main()
