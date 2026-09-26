"""Frozen A10 growth-block decomposition and joint plate/target shift stress test."""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
OUT = HERE / 'tables/followup_v1'
SPEC = HERE / 'config/a10_followup_models.json'
OUTPUTS = ['model_scores.tsv', 'model_predictions.tsv', 'model_fold_errors.tsv',
           'model_comparisons.tsv', 'model_fold_comparisons.tsv', 'model_run.json']


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def train_predict(train_x, train_y, test_x):
    """No test outcome is accepted, and feature transforms are training-only."""
    import numpy as np
    mean, sd = train_x.mean(axis=0), train_x.std(axis=0)
    keep = sd > 0
    a = np.column_stack([np.ones(len(train_x)), (train_x[:, keep] - mean[keep]) / sd[keep]])
    b = np.column_stack([np.ones(len(test_x)), (test_x[:, keep] - mean[keep]) / sd[keep]])
    beta, _, rank, _ = np.linalg.lstsq(a, train_y, rcond=None)
    assert rank == a.shape[1], 'Training design is rank deficient'
    return b @ beta


def centre(values, groups):
    import pandas as pd
    frame = pd.DataFrame(values)
    return (frame - frame.groupby(pd.Series(groups)).transform('mean')).to_numpy()


def legacy_r2(x, y, groups):
    import numpy as np
    pred = np.full(len(y), np.nan)
    for group in np.unique(groups):
        tr, te = groups != group, groups == group
        keep = x[tr].std(axis=0) > 0
        keep[0] = True
        beta, *_ = np.linalg.lstsq(x[tr][:, keep], y[tr], rcond=None)
        pred[te] = x[te][:, keep] @ beta
    return 1 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum()


def reproduce_grid(data, growth):
    import numpy as np
    import pandas as pd
    expected = pd.read_csv(HERE / 'tables/revised_grid.tsv', sep='\t')
    checks = []
    groups = data.unit.to_numpy()
    for programme, blocks in [('frozen_repair_modules', {'mouse': [c for c in data if c.startswith('epi__')],
                                                      'human': [c for c in data if c.startswith('fib__')]}),
                              ('growth_hallmark', growth)]:
        for baseline, within in [('pooled', False), ('within_unit', True)]:
            y = np.log2(data.organoids_area_mean_day14.to_numpy() + 1)
            cont = np.column_stack([np.log2(data.organoids_area_mean_day07 + 1),
                                    data.epi_fraction, np.log2(data.epi_total + 1)])
            if within:
                y = centre(y, groups).ravel()
                cont = centre(cont, groups)
            else:
                cont = np.column_stack([cont, pd.get_dummies(data.plate, drop_first=True).to_numpy(float)])
            base = np.column_stack([np.ones(len(data)), cont])
            matrices = []
            for species in ('mouse', 'human'):
                m = data[blocks[species]].to_numpy(float)
                if within:
                    m = centre(m, groups)
                sd = m.std(axis=0)
                matrices.append((m - m.mean(axis=0)) / np.where(sd > 0, sd, 1))
            b = legacy_r2(base, y, groups)
            e = legacy_r2(np.column_stack([base, matrices[0]]), y, groups)
            f = legacy_r2(np.column_stack([base, *matrices]), y, groups)
            row = expected[(expected.programme_set == programme) & (expected.baseline == baseline)].iloc[0]
            for metric, observed in [('baseline_r2', b), ('epithelial_increment', e-b), ('fibroblast_increment', f-e)]:
                error = abs(float(row[metric]) - float(observed))
                assert error < 0.00005, (programme, baseline, metric, error)
                checks.append({'programme': programme, 'baseline': baseline, 'metric': metric,
                               'expected_rounded': float(row[metric]), 'observed': float(observed), 'absolute_error': error})
    return checks


def main():
    import numpy as np
    import pandas as pd

    parser = argparse.ArgumentParser()
    parser.add_argument('--source-cache', type=Path, required=True)
    parser.add_argument('--data-root', type=Path, required=True)
    args = parser.parse_args()
    assert not any((OUT / name).exists() for name in OUTPUTS), 'Refusing to overwrite'
    spec = json.loads(SPEC.read_text())
    extraction = json.loads((HERE / 'tables/hallmark_extract_run.json').read_text())
    cache = args.source_cache / 'a10_hallmark_counts.npz'
    assert sha(cache) == extraction['cache']['sha256']
    npz = np.load(cache, allow_pickle=True)
    wanted = [s for group in spec['programme_blocks'].values() for s in group]
    scores = pd.DataFrame({'library name': [str(x) for x in npz['libraries']]})
    assert scores['library name'].is_unique
    input_hashes = {'hallmark_npz': sha(cache), 'original_module_scores': sha(HERE / 'tables/module_scores.tsv')}
    growth = {}
    coverage = {}
    for species, prefix, suffix in [('mouse', 'mh', 'Mm'), ('human', 'h', 'Hs')]:
        gmt = args.data_root / f'raw_data/msigdb/{prefix}.all.v2024.1.{suffix}.symbols.gmt'
        assert sha(gmt) == extraction['gmt'][species]['sha256']
        input_hashes[species + '_gmt'] = sha(gmt)
        genes = {str(g): i for i, g in enumerate(npz[species + '_genes'])}
        values = np.log2(npz[species + '_counts'] / np.maximum(npz[species + '_totals'], 1) * 1e6 + 1)
        sets = {p[0]: p[2:] for line in gmt.read_text().splitlines() if (p := line.split('\t'))[0] in wanted}
        assert set(sets) == set(wanted)
        growth[species] = []
        for name in wanted:
            indices = [genes[g] for g in sets[name] if g in genes]
            assert len(indices) / len(sets[name]) >= 0.7
            col = species + '__' + name
            scores[col] = values[indices].mean(axis=0)
            growth[species].append(col)
            coverage[col] = {'assayed': len(indices), 'source': len(sets[name])}
    data = pd.read_csv(HERE / 'tables/module_scores.tsv', sep='\t').merge(scores, on='library name', validate='one_to_one')
    assert len(data) == 886
    data = data[data.organoids_area_mean_day07.gt(0) & data.organoids_area_mean_day14.gt(0) & data.epi_total.gt(0)].reset_index(drop=True)
    assert len(data) == 885
    data['plate'] = data.unit.str.split('-rep').str[0]
    instrument = reproduce_grid(data, growth)

    basecols = ['day7', 'epi_fraction', 'epi_size']
    pcols = ['mouse__' + n for n in spec['programme_blocks']['proliferation']]
    rcols = ['mouse__' + n for n in spec['programme_blocks']['remaining_growth']]
    models = {'baseline': basecols, 'baseline+proliferation': basecols + pcols,
              'baseline+remaining_growth': basecols + rcols,
              'baseline+proliferation+remaining_growth': basecols + pcols + rcols}
    predictions, errors, comparisons, fold_comparisons = [], [], [], []
    for setting in spec['sensitivity_settings']:
        sub = data.copy()
        if setting in ('drop_TIGIT', 'drop_both'):
            sub = sub[sub.target != 'TIGIT']
        if setting in ('drop_TDTOMATO', 'drop_both'):
            sub = sub[sub.target != 'TDTOMATO']
        sub = sub.reset_index(drop=True)
        for scale in spec['outcome_scales']:
            frame = sub.copy()
            y = frame.organoids_area_mean_day14.to_numpy(float)
            frame['day7'] = frame.organoids_area_mean_day07
            if scale == 'inherited_log2':
                y = np.log2(y + 1)
                frame['day7'] = np.log2(frame.day7 + 1)
            frame['epi_size'] = np.log2(frame.epi_total + 1)
            allcols = basecols + pcols + rcols
            raw_x = frame[allcols].to_numpy(float)
            assert np.isfinite(raw_x).all() and np.isfinite(y).all()
            for evaluation in spec['evaluation']:
                if evaluation not in ('within_group', 'plate_shift'):
                    continue
                folds = frame.unit.to_numpy() if evaluation == 'within_group' else frame.plate.to_numpy()
                xx = centre(raw_x, frame.unit.to_numpy()) if evaluation == 'within_group' else raw_x
                yy = centre(y, frame.unit.to_numpy()).ravel() if evaluation == 'within_group' else y.copy()
                pred_frame = frame[['library name', 'unit', 'plate', 'target']].copy()
                pred_frame['setting'], pred_frame['scale'], pred_frame['evaluation'] = setting, scale, evaluation
                pred_frame['outcome'] = yy
                fold_errors = {}
                for model, columns in models.items():
                    x = xx[:, [allcols.index(c) for c in columns]]
                    pred = np.full(len(frame), np.nan)
                    for fold in sorted(set(folds)):
                        tr, te = folds != fold, folds == fold
                        pred[te] = train_predict(x[tr], yy[tr], x[te])
                        residual = yy[te] - pred[te]
                        sse = float(residual @ residual)
                        fold_errors[(model, fold)] = sse
                        unseen = set(frame.loc[te, 'target']) - set(frame.loc[tr, 'target'])
                        errors.append({'setting': setting, 'scale': scale, 'evaluation': evaluation, 'fold': fold,
                                       'model': model, 'wells': int(te.sum()), 'sse': sse, 'mse': sse / te.sum(),
                                       'heldout_targets': int(frame.loc[te, 'target'].nunique()), 'unseen_targets': len(unseen),
                                       'wells_with_unseen_target': int(frame.loc[te, 'target'].isin(unseen).sum())})
                    assert np.isfinite(pred).all()
                    pred_frame[model] = pred
                predictions.append(pred_frame)
                for item in spec['comparisons']:
                    fold_gains = []
                    ref_total = cand_total = 0.0
                    for fold in sorted(set(folds)):
                        ref = fold_errors[(item['reference'], fold)]
                        cand = fold_errors[(item['candidate'], fold)]
                        assert ref > 0
                        gain = 1 - cand / ref
                        ref_total += ref
                        cand_total += cand
                        fold_gains.append(gain)
                        fold_comparisons.append({'setting': setting, 'scale': scale, 'evaluation': evaluation,
                                                 'comparison': item['name'], 'fold': fold, 'relative_error_reduction': gain})
                    comparisons.append({'setting': setting, 'scale': scale, 'evaluation': evaluation,
                                        'comparison': item['name'], 'wells': len(frame), 'folds': len(fold_gains),
                                        'reference_sse': ref_total, 'candidate_sse': cand_total,
                                        'relative_error_reduction': 1-cand_total/ref_total,
                                        'equal_fold_mean_reduction': float(np.mean(fold_gains)),
                                        'positive_folds': int(sum(g > 0 for g in fold_gains)),
                                        'minimum_fold_reduction': float(min(fold_gains)), 'maximum_fold_reduction': float(max(fold_gains))})
    summaries = pd.DataFrame(comparisons)
    margin = spec['decisions']['descriptive_margin_relative_error_reduction']
    verdicts = {}
    for item in spec['comparisons']:
        rows = summaries[(summaries.comparison == item['name']) & summaries.setting.isin(['primary', 'drop_both'])]
        within = rows[rows.evaluation == 'within_group']
        plates = rows[rows.evaluation == 'plate_shift']
        assert len(within) == len(plates) == 4
        verdicts[item['name']] = {
            'within_group_margin_met': bool(within.relative_error_reduction.ge(margin).all()),
            'consistent_plate_shift_gain': bool(plates.relative_error_reduction.ge(margin).all() and plates.minimum_fold_reduction.gt(0).all())}
    frames = {'model_scores.tsv': scores, 'model_predictions.tsv': pd.concat(predictions, ignore_index=True),
              'model_fold_errors.tsv': pd.DataFrame(errors), 'model_comparisons.tsv': summaries,
              'model_fold_comparisons.tsv': pd.DataFrame(fold_comparisons)}
    record = {'stage': 'followup_models', 'completed_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
              'script_sha256': sha(__file__), 'spec_sha256': sha(SPEC), 'input_hashes': input_hashes,
              'coverage': coverage, 'instrument_reproduction': instrument, 'verdicts': verdicts,
              'numpy': np.__version__, 'pandas': pd.__version__,
              'interpretation': 'same-screen descriptive comparison; biological units unresolved; plate shift also changes target mix'}
    for name, frame in frames.items():
        frame.to_csv(OUT / name, sep='\t', index=False, float_format='%.12g')
    record['outputs'] = {name: sha(OUT / name) for name in frames}
    (OUT / 'model_run.json').write_text(json.dumps(record, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({'instrument_checks': len(instrument), 'verdicts': verdicts}))
    selected = summaries[(summaries.setting.isin(['primary', 'drop_both'])) & (summaries.scale == 'inherited_log2')]
    print(selected[['setting', 'evaluation', 'comparison', 'relative_error_reduction', 'positive_folds', 'folds']].to_string(index=False))


if __name__ == '__main__':
    main()
