"""Audit A10 plate/target/imaging allocation and all GEO sample metadata; no fit."""
from __future__ import annotations

import argparse
import collections
import datetime as dt
import gzip
import hashlib
import itertools
import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
OUT = HERE / 'tables' / 'followup_v1'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def soft_records(path):
    records, current = [], None
    with gzip.open(path, 'rt', encoding='utf-8') as handle:
        for line in handle:
            if line.startswith('^SAMPLE = '):
                current = collections.defaultdict(list)
                current['accession'] = [line.strip().split(' = ', 1)[1]]
                records.append(current)
            elif line.startswith('^'):
                current = None
            elif current is not None and line.startswith('!Sample_') and ' = ' in line:
                key, value = line.rstrip('\r\n').split(' = ', 1)
                current[key].append(value)
    return records


def main():
    import numpy as np
    import pandas as pd

    parser = argparse.ArgumentParser()
    parser.add_argument('--source-cache', type=Path, required=True)
    args = parser.parse_args()
    if (OUT / 'diagnostic_run.json').exists() or list(OUT.glob('diagnostic_*.tsv')):
        raise SystemExit('Refusing to overwrite diagnostic outputs')
    prior = json.loads((HERE / 'tables/audit_run.json').read_text())
    inputs = {}
    for rel, meta in prior['inputs'].items():
        path = args.source_cache / Path(rel).name
        assert sha(path) == meta['sha256'], path
        inputs[path.name] = {'bytes': path.stat().st_size, 'sha256': sha(path)}
    soft = HERE / 'cache/GSE307112_family.soft.gz'
    retrieval = json.loads((HERE / 'cache/geo_followup_retrieval.json').read_text())
    assert sha(soft) == retrieval['sha256']
    inputs[soft.name] = retrieval
    for name in ('module_scores.tsv', 'well_join.tsv', 'revised_per_unit.tsv'):
        inputs[name] = {'sha256': sha(HERE / 'tables' / name)}

    data = pd.read_csv(HERE / 'tables/module_scores.tsv', sep='\t')
    joined = pd.read_csv(HERE / 'tables/well_join.tsv', sep='\t')
    qc = pd.read_csv(args.source_cache / 'GSE307112_xenome_stats.csv.gz')
    imaging = pd.read_csv(args.source_cache / 'GSE307112_imaging_outputs.csv.gz')
    design = pd.read_csv(args.source_cache / 'GSE307112_plate_design.csv.gz')
    assert len(data) == len(joined) == len(qc) == 886
    assert data['library name'].is_unique and qc['library name'].is_unique
    assert set(data['library name']) == set(qc['library name']) == set(joined['library name'])
    data['plate'] = data.unit.str.split('-rep').str[0]
    imaging['unit'] = [f'{p}-rep{str(r).split("-")[-1]}' for p, r in zip(imaging.plate, imaging.plate_replicate)]
    assert not imaging.duplicated(['unit', 'well', 'day']).any()
    area = imaging.pivot(index=['unit', 'well'], columns='day', values='organoids_area_mean')
    area.columns = ['raw_' + str(c) for c in area.columns]
    check = data.merge(area.reset_index(), left_on=['unit', 'well_key'], right_on=['unit', 'well'], validate='one_to_one')
    for day in ('day07', 'day14'):
        assert np.allclose(check['organoids_area_mean_' + day], check['raw_' + day], equal_nan=True)
    qcols = [c for c in qc if c.startswith('xenome_numReads') and c != 'xenome_numReadsInput']
    q = qc[['library name', 'xenome_numReadsInput'] + qcols].copy()
    for col in qcols:
        q[col + '_fraction'] = q[col] / q.xenome_numReadsInput
    data = data.merge(q[['library name'] + [c + '_fraction' for c in qcols]], on='library name', validate='one_to_one')
    data['eligible'] = data.organoids_area_mean_day07.gt(0) & data.organoids_area_mean_day14.gt(0) & data.epi_total.gt(0)
    assert int(data.eligible.sum()) == 885
    data['log2_area_change'] = np.log2(data.organoids_area_mean_day14 + 1) - np.log2(data.organoids_area_mean_day07 + 1)

    output = {}
    output['diagnostic_target_plate.tsv'] = pd.crosstab(data.target, data.plate).reset_index()
    top = data.target.value_counts().idxmax()
    overlap = []
    for scope, frame in [('all', data), ('drop_top_and_TDTOMATO', data[~data.target.isin([top, 'TDTOMATO'])])]:
        sets = {p: set(g.target) for p, g in frame.groupby('plate')}
        for p, q in itertools.combinations(sorted(sets), 2):
            common = sorted(sets[p] & sets[q])
            overlap.append({'scope': scope, 'plate_a': p, 'plate_b': q, 'targets_a': len(sets[p]),
                            'targets_b': len(sets[q]), 'common_targets': len(common), 'common_names': '|'.join(common)})
    output['diagnostic_target_overlap.tsv'] = pd.DataFrame(overlap)

    summaries = []
    metrics = ['epi_total', 'fib_total', 'epi_fraction', 'xenome_numReadsInput',
               'organoids_area_mean_day07', 'organoids_area_mean_day14',
               'organoids_count_day07', 'organoids_count_day14', 'log2_area_change'] + [c + '_fraction' for c in qcols]
    for scope in ('plate', 'unit'):
        for key, group in data.groupby(scope):
            row = {'scope': scope, 'group': key, 'libraries': len(group), 'eligible': int(group.eligible.sum()),
                   'targets': group.target.nunique(), 'TDTOMATO': int(group.target.eq('TDTOMATO').sum()),
                   'top_target': top, 'top_target_wells': int(group.target.eq(top).sum())}
            for name in metrics:
                values = group[name].dropna()
                for label, value in [('min', values.min()), ('median', values.median()), ('max', values.max())]:
                    row[name + '_' + label] = value
            summaries.append(row)
    output['diagnostic_group_summary.tsv'] = pd.DataFrame(summaries)
    output['diagnostic_segmentation.tsv'] = imaging.groupby(['plate', 'day', 'segmentation_model'], dropna=False).size().rename('imaging_rows').reset_index()
    imrows = []
    for (plate, day), group in imaging.groupby(['plate', 'day']):
        # Diagnostic only: equality to an imaged area requires untransformed mean area.
        implied = group.organoids_area_mean * group.organoids_count / group.organoids_area_prop.replace(0, np.nan)
        imrows.append({'plate': plate, 'day': day, 'rows': len(group), 'positive_count': int(group.organoids_count.gt(0).sum()),
                       'area_prop_median': group.organoids_area_prop.median(),
                       'area_times_count_over_coverage_min': implied.min(),
                       'area_times_count_over_coverage_median': implied.median(),
                       'area_times_count_over_coverage_max': implied.max()})
    output['diagnostic_imaging_arithmetic.tsv'] = pd.DataFrame(imrows)

    guide_rows = []
    for (plate, gene), group in design.groupby(['plate', 'gene_name']):
        guides = {str(v) for c in ('sgrna1', 'sgrna2', 'sgrna3') for v in group[c].dropna()}
        guide_rows.append({'plate': plate, 'gene_name': gene, 'layout_rows': len(group),
                           'unique_guide_sequences': len(guides), 'guide_tuple_count': len(group[['sgrna1', 'sgrna2', 'sgrna3']].drop_duplicates())})
    output['diagnostic_guide_layout.tsv'] = pd.DataFrame(guide_rows)

    records = soft_records(soft)
    samples, catalogue = [], collections.Counter()
    for record in records:
        chars = {}
        for value in record.get('!Sample_characteristics_ch1', []):
            key, sep, val = value.partition(': ')
            chars[key] = val if sep else ''
        library = [s.removeprefix('Library name: ') for s in record.get('!Sample_description', []) if s.startswith('Library name: ')]
        assert len(library) == 1, record['accession']
        sample = {'gsm': record['accession'][0], 'library_name': library[0],
                  'title': '|'.join(record.get('!Sample_title', [])), **chars}
        samples.append(sample)
        for field, values in record.items():
            if field.startswith('!Sample_') and any(x in field for x in ('protocol', 'processing', 'characteristics', 'source_name', 'description')):
                for value in values:
                    catalogue[(field, value)] += 1
    assert len(records) == 886 and len({s['gsm'] for s in samples}) == 886
    sample_table = pd.DataFrame(samples)
    assert sample_table.library_name.is_unique
    assert set(sample_table.library_name) == set(data['library name']), 'GEO descriptions differ from libraries'
    output['diagnostic_geo_sample_fields.tsv'] = sample_table
    output['diagnostic_geo_field_catalogue.tsv'] = pd.DataFrame([
        {'field': k, 'value': v, 'sample_count': n} for (k, v), n in sorted(catalogue.items())])
    structured_identity = [c for c in sample_table if re.search(r'donor|isolation|animal|mouse.?id|lot|preparation', c, re.I)]
    multi_plate = output['diagnostic_target_plate.tsv'].set_index('target').gt(0).sum(axis=1)
    record = {'stage': 'followup_design_audit', 'completed_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
              'script_sha256': sha(__file__), 'plan_sha256': sha(HERE / 'FOLLOWUP_PLAN.md'),
              'inputs': inputs, 'libraries': len(data), 'eligible': int(data.eligible.sum()),
              'geo_sample_records': len(records), 'geo_description_library_ids_match': True,
              'technical_attempt': 'Initial assertion expected titles to be library IDs and stopped before output. Corrected to exact Library name: description fields; no cohort or criterion changed.',
              'structured_sample_fields': list(sample_table.columns), 'structured_preparation_fields': structured_identity,
              'targets': int(data.target.nunique()), 'targets_on_multiple_plates': int(multi_plate.gt(1).sum()),
              'shared_target_names': list(multi_plate[multi_plate.gt(1)].index), 'most_represented_target': top,
              'target_overlap_after_removals': [r for r in overlap if r['scope'] != 'all'],
              'segmentation_models': sorted(imaging.segmentation_model.astype(str).unique().tolist()),
              'no_model_fitted': True, 'python_libraries': {'numpy': np.__version__, 'pandas': pd.__version__}}
    OUT.mkdir(parents=True, exist_ok=True)
    for name, frame in output.items():
        assert not (OUT / name).exists()
        frame.to_csv(OUT / name, sep='\t', index=False)
    record['outputs'] = {name: sha(OUT / name) for name in output}
    (OUT / 'diagnostic_run.json').write_text(json.dumps(record, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({k: record[k] for k in ('libraries', 'eligible', 'geo_sample_records', 'structured_sample_fields',
                     'structured_preparation_fields', 'targets_on_multiple_plates', 'shared_target_names', 'most_represented_target', 'segmentation_models')}))


if __name__ == '__main__':
    main()
