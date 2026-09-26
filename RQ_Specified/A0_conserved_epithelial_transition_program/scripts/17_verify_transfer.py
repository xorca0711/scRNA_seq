"""Verify transfer aggregation and independently calculate tied-rank scores from counts."""
import csv
import gzip
import hashlib
import json
import math
import statistics
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

BASE = Path(__file__).resolve().parents[1]
OUT = BASE / 'tables/pilot_v1'
WORK = BASE / 'processed/pilot_v1'


def read_rows(path):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(4*1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()


def verify():
    run = json.loads((OUT / 'transfer_run.json').read_text())
    checks = 0
    for key, path in [('script_sha256', BASE/'scripts/13_transfer_pilot.py'),
                      ('configuration_sha256', BASE/'config/pilot_v1.json'),
                      ('module_sha256', OUT/'frozen_programme.json'),
                      ('preparation_sha256', OUT/'preparation.json'),
                      ('execution_plan_sha256', BASE/'TRANSFER_EXECUTION.md')]:
        assert run[key] == sha(path); checks += 1
    for name, digest in run['outputs'].items():
        assert sha(OUT/name) == digest; checks += 1
    cells = read_rows(OUT/'transfer_cell_scores.tsv.gz')
    grouped = defaultdict(list)
    for cell in cells:
        score = float(cell['programme_score'])
        assert 0 < score <= 1
        grouped[(cell['role'], cell['unit'], cell['state'])].append(score)
    means = {k: math.fsum(v)/len(v) for k,v in grouped.items()}
    max_aggregation_error = 0
    for row in read_rows(OUT/'transfer_unit_scores.tsv'):
        key = (row['role'], row['unit'], row['state'])
        assert len(grouped[key]) == int(row['cells'])
        error = abs(means[key]-float(row['programme_score']))
        max_aggregation_error = max(max_aggregation_error, error)
        assert error < 1e-12; checks += 2
    differences = defaultdict(dict)
    for row in read_rows(OUT/'transfer_unit_differences.tsv'):
        role, unit, endpoint = row['role'], row['unit'], row['endpoint']
        expected = means[(role,unit,'intermediate')]-means[(role,unit,endpoint)]
        assert abs(expected-float(row['difference'])) < 1e-12; checks += 1
        differences[(role,endpoint)][unit] = expected
    v1_pass = True
    for row in read_rows(OUT/'transfer_summary.tsv'):
        values = differences[(row['role'],row['endpoint'])]
        med = statistics.median(values.values())
        fraction = sum(v>0 for v in values.values())/len(values)
        omitted = [statistics.median(v for u,v in values.items() if u!=omit) for omit in values]
        assert len(values) == int(row['units'])
        assert abs(med-float(row['median_difference'])) < 1e-12
        assert abs(fraction-float(row['positive_fraction'])) < 1e-12
        assert abs(min(omitted)-float(row['minimum_leave_one_out_median'])) < 1e-12
        passed = med>0 and fraction>=2/3 and min(omitted)>0
        assert passed == (row['criterion_met']=='True'); checks += 5
        if row['role']=='V1':
            v1_pass &= passed
    for row in read_rows(OUT/'transfer_leave_one_out.tsv'):
        values = differences[(row['role'],row['endpoint'])]
        med = statistics.median(v for u,v in values.items() if u!=row['omitted_unit'])
        assert abs(med-float(row['median_difference'])) < 1e-12; checks += 1
    days = {r['unit']:int(r['day']) for r in cells if r['role']=='D1'}
    for row in read_rows(OUT/'repair_day_sensitivity.tsv'):
        values = differences[('D1',row['endpoint'])]
        if row['sensitivity']=='equal_day_mean':
            per_day = defaultdict(list)
            for u,v in values.items():
                per_day[days[u]].append(v)
            expected = statistics.mean(statistics.mean(v) for v in per_day.values())
        else:
            expected = statistics.median(v for u,v in values.items() if days[u]!=int(row['omitted_day']))
        assert abs(expected-float(row['effect'])) < 1e-12; checks += 1
    assert v1_pass == run['V1_primary_rule_passed']; checks += 1
    if not v1_pass:
        assert run['status']=='STOP_PRIMARY_TRANSFER_NOT_SUPPORTED'
        assert run['specificity']=='pruned_after_primary_transfer_failure_as_declared_before_V1_scores'
    # A separate score calculation counts values below/equal to each programme
    # gene; it does not call scipy.rankdata or import the scoring implementation.
    ortho = read_rows(OUT/'ortholog_universe.tsv')
    module = json.loads((OUT/'frozen_programme.json').read_text())
    module_rows = [i for i,r in enumerate(ortho) if r['human'] in module['human_genes']]
    probes = 0; max_score_error = 0
    for role in ['D1','D2','V1']:
        native = json.loads((WORK/f'{role}_genes.json').read_text())
        index = {g:i for i,g in enumerate(native)}
        metadata = read_rows(WORK/f'{role}_cells.tsv')
        col_index = {r['cell_id']:i for i,r in enumerate(metadata)}
        native_rows = [index[r['human' if role=='D2' else 'mouse']] for r in ortho]
        x = np.memmap(WORK/f'{role}_counts.bin',mode='r',dtype=np.int32,shape=(len(native),len(metadata)))
        reduced = np.asarray(x[native_rows,:],dtype=np.int32); del x
        seen = set()
        for cell in cells:
            key = (cell['unit'],cell['state'])
            if cell['role'] != role or key in seen:
                continue
            seen.add(key)
            counts = reduced[:,col_index[cell['cell_id']]]
            ranks = [(int(np.count_nonzero(counts < counts[i])) +
                      (int(np.count_nonzero(counts == counts[i]))+1)/2)/len(ortho) for i in module_rows]
            score = math.fsum(ranks)/len(ranks)
            error = abs(score-float(cell['programme_score']))
            assert error < 1e-12
            max_score_error = max(max_score_error,error); probes += 1
        del reduced
    return dict(status='PASS',completed_utc=datetime.now(timezone.utc).isoformat(),
                verifier_sha256=sha(Path(__file__)),transfer_run_sha256=sha(OUT/'transfer_run.json'),
                arithmetic_checks=checks,cells_aggregated=len(cells),raw_score_probes=probes,
                maximum_aggregation_error=max_aggregation_error,maximum_raw_score_error=max_score_error,
                V1_primary_rule_passed=v1_pass)


if __name__ == '__main__':
    destination = OUT/'transfer_verification.json'
    if destination.exists():
        raise SystemExit('Refusing to overwrite verification')
    record = verify()
    destination.write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record))
