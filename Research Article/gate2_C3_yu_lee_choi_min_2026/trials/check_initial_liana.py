"""Check initial LR output integrity without rerunning expression analysis."""
import csv
import hashlib
import json
import math
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

OUT = Path(__file__).resolve().parent / 'u5_ipf_liana'


def rows(path):
    with path.open(newline='', encoding='utf-8') as handle:
        return list(csv.DictReader(handle))


def main():
    record = json.loads((OUT / 'run_record.json').read_text())
    spec = json.loads((OUT / 'specification.json').read_text())
    assert record['status'] == 'completed_descriptive_LR_requires_interpretation'
    assert record['input_count_parity'] == 'passed'
    assert sorted(record['completed_donors']) == sorted(spec['donors'])
    pairs = rows(OUT / 'fixed_donor_pairs.csv')
    allowed = {(r['donor'], r['source'], r['target']) for r in pairs}
    diseases = {r['donor']: r['disease'] for r in pairs}
    assert len(diseases) == len(spec['donors'])
    groups = rows(OUT / 'selected_groups.csv')
    assert all(int(r['cells']) >= spec['cell_floor'] for r in groups)
    files = []
    for resource in spec['resources']:
        resource_edges = {(r['ligand'], r['receptor'])
                          for r in rows(OUT / ('resource_' + resource + '.csv'))}
        for donor in spec['donors']:
            path = OUT / (resource + '_' + donor + '.csv')
            result = rows(path)
            keys = set()
            for r in result:
                assert r['donor'] == donor and r['resource'] == resource
                assert r['disease'] == diseases[donor]
                assert (donor, r['source'], r['target']) in allowed
                assert (r['ligand_complex'], r['receptor_complex']) in resource_edges
                key = (r['source'], r['target'], r['ligand_complex'], r['receptor_complex'])
                assert key not in keys
                keys.add(key)
                for column in ['ligand_props', 'receptor_props', 'lr_probs']:
                    value = float(r[column])
                    assert math.isfinite(value) and 0 <= value <= 1
                assert float(r['ligand_props']) >= spec['expression_fraction'] - 1e-9
                assert float(r['receptor_props']) >= spec['expression_fraction'] - 1e-9
                assert not r.get('cellchat_pvals'), 'Cell-permutation P values were disabled'
            files.append({'file': path.name, 'rows': len(result),
                          'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
    report = {
        'checked_utc': datetime.now(timezone.utc).isoformat(),
        'status': 'passed',
        'scope': 'file completeness, fixed-pair/resource membership, row uniqueness, '
                 'score bounds, detection eligibility and disabled cell P values; '
                 'not biological validation or missing-edge interpretation',
        'donors': len(spec['donors']), 'donors_by_disease': dict(Counter(diseases.values())),
        'output_files': len(files), 'empty_files': [f['file'] for f in files if not f['rows']],
        'rows_by_resource': {resource: sum(f['rows'] for f in files
                                         if f['file'].startswith(resource + '_'))
                             for resource in spec['resources']},
        'spec_sha256': hashlib.sha256((OUT / 'specification.json').read_bytes()).hexdigest(),
        'files': files,
    }
    (OUT / 'validation.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({k: v for k, v in report.items() if k != 'files'}))


if __name__ == '__main__':
    main()
