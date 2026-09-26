"""Apply the committed programme once; decide whether specificity testing is warranted."""
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

BASE = Path(__file__).resolve().parents[1]
ROOT = BASE.parents[1]
OUT = BASE / 'tables/pilot_v1'
WORK = BASE / 'processed/pilot_v1'


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(4*1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    targets = ['transfer_cell_scores.tsv.gz', 'transfer_unit_scores.tsv',
               'transfer_unit_differences.tsv', 'transfer_summary.tsv',
               'transfer_leave_one_out.tsv', 'repair_day_sensitivity.tsv', 'transfer_run.json']
    if any((OUT / name).exists() for name in targets):
        raise SystemExit('Refusing to overwrite transfer')
    cfg = json.loads((BASE / 'config/pilot_v1.json').read_text())
    module = json.loads((OUT / 'frozen_programme.json').read_text())
    assert module['status'] == 'FROZEN_FOR_TRANSFER' and 20 <= module['selected_genes'] <= 50
    assert module['config_sha256'] == sha(BASE / 'config/pilot_v1.json')
    relative = (OUT / 'frozen_programme.json').relative_to(ROOT).as_posix()
    assert subprocess.check_output(['git', 'show', 'HEAD:'+relative], cwd=ROOT) == (OUT / 'frozen_programme.json').read_bytes()
    prep = json.loads((OUT / 'preparation.json').read_text())
    for name, record in prep['prepared'].items():
        assert sha(WORK / name) == record['sha256'], name
    universe = pd.read_csv(OUT / 'ortholog_universe.tsv', sep='\t')
    assert sha(OUT / 'ortholog_universe.tsv') == module['ortholog_universe_sha256']
    selected = np.flatnonzero(universe.human.isin(module['human_genes']))
    assert len(selected) == module['selected_genes']
    all_cells = []
    for role in ['D1', 'D2', 'V1']:
        native = json.loads((WORK / f'{role}_genes.json').read_text())
        cells = pd.read_csv(WORK / f'{role}_cells.tsv', sep='\t')
        index = {g: i for i, g in enumerate(native)}
        rows = [index[g] for g in universe['human' if role == 'D2' else 'mouse']]
        x = np.memmap(WORK / f'{role}_counts.bin', mode='r', dtype=np.int32, shape=(len(native), len(cells)))
        # Read the reduced matrix once, avoiding repeated network-backed mmap scans.
        counts = np.asarray(x[rows, :], dtype=np.int32)
        del x
        score = np.empty(len(cells))
        for begin in range(0, len(cells), 128):
            ranks = rankdata(counts[:, begin:begin+128], method='average', axis=0) / len(universe)
            score[begin:begin+ranks.shape[1]] = ranks[selected, :].mean(axis=0)
        del counts
        cells['role'] = role
        cells['programme_score'] = score
        fields = ['role', 'cell_id', 'unit', 'state', 'library_total', 'programme_score']
        if role == 'D1':
            fields.append('day')
        all_cells.append(cells[fields])
        print(f'{role}: scored {len(cells)} cells', flush=True)
    cells = pd.concat(all_cells, ignore_index=True)
    unit_scores = cells.groupby(['role', 'unit', 'state'], as_index=False).agg(
        cells=('programme_score', 'size'), programme_score=('programme_score', 'mean'))
    differences = []
    for (role, unit), group in unit_scores.groupby(['role', 'unit']):
        group = group.set_index('state')
        for endpoint in ['start', 'destination']:
            assert min(group.loc['intermediate', 'cells'], group.loc[endpoint, 'cells']) >= cfg['cell_floor']
            differences.append(dict(role=role, unit=unit, endpoint=endpoint,
                                    difference=float(group.loc['intermediate', 'programme_score'] - group.loc[endpoint, 'programme_score'])))
    differences = pd.DataFrame(differences)
    summaries = []
    leave_out = []
    for (role, endpoint), group in differences.groupby(['role', 'endpoint']):
        assert len(group) >= cfg['unit_floor']
        med = float(group.difference.median())
        fraction = float((group.difference > 0).mean())
        omitted = []
        for unit in group.unit:
            value = float(group.loc[group.unit != unit, 'difference'].median())
            omitted.append(value)
            leave_out.append(dict(role=role, endpoint=endpoint, omitted_unit=unit, median_difference=value))
        passed = med > 0 and fraction >= 2/3 and min(omitted) > 0
        summaries.append(dict(role=role, endpoint=endpoint, units=len(group), median_difference=med,
                              positive_units=int((group.difference > 0).sum()), positive_fraction=fraction,
                              minimum_leave_one_out_median=min(omitted), criterion_met=passed))
    summary = pd.DataFrame(summaries)
    days = cells[cells.role == 'D1'][['unit', 'day']].drop_duplicates()
    assert days.unit.is_unique
    day_differences = differences[differences.role == 'D1'].merge(days, on='unit', validate='many_to_one')
    day_rows = []
    for endpoint, group in day_differences.groupby('endpoint'):
        day_rows.append(dict(endpoint=endpoint, sensitivity='equal_day_mean', omitted_day='',
                             effect=float(group.groupby('day').difference.mean().mean())))
        for day in sorted(group.day.unique()):
            day_rows.append(dict(endpoint=endpoint, sensitivity='leave_one_day_out_median', omitted_day=int(day),
                                 effect=float(group.loc[group.day != day, 'difference'].median())))
    cells.to_csv(OUT / targets[0], sep='\t', index=False, float_format='%.17g', compression={'method':'gzip','mtime':0})
    unit_scores.to_csv(OUT / targets[1], sep='\t', index=False, float_format='%.17g')
    differences.to_csv(OUT / targets[2], sep='\t', index=False, float_format='%.17g')
    summary.to_csv(OUT / targets[3], sep='\t', index=False, float_format='%.17g')
    pd.DataFrame(leave_out).to_csv(OUT / targets[4], sep='\t', index=False, float_format='%.17g')
    pd.DataFrame(day_rows).to_csv(OUT / targets[5], sep='\t', index=False, float_format='%.17g')
    passed = bool(summary.loc[summary.role == 'V1', 'criterion_met'].all())
    record = dict(completed_utc=datetime.now(timezone.utc).isoformat(),
                  git_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                  script_sha256=sha(Path(__file__)), configuration_sha256=sha(BASE/'config/pilot_v1.json'),
                  execution_plan_sha256=sha(BASE/'TRANSFER_EXECUTION.md'), module_sha256=sha(OUT/'frozen_programme.json'),
                  preparation_sha256=sha(OUT/'preparation.json'), V1_primary_rule_passed=passed,
                  status='TRANSFER_PASSED_SPECIFICITY_PENDING' if passed else 'STOP_PRIMARY_TRANSFER_NOT_SUPPORTED',
                  specificity='pending' if passed else 'pruned_after_primary_transfer_failure_as_declared_before_V1_scores',
                  outputs={name:sha(OUT/name) for name in targets if name!='transfer_run.json'})
    (OUT / targets[-1]).write_text(json.dumps(record, indent=2)+'\n')
    print(json.dumps({'status':record['status'],'summaries':summaries}))


if __name__ == '__main__':
    main()
