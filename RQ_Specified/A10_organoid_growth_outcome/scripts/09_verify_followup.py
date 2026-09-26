"""Verify saved arithmetic and reproduce plate predictions with independent QR fits."""
from __future__ import annotations

import datetime as dt
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
OUT = HERE / 'tables/followup_v1'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    import numpy as np
    import pandas as pd
    from scipy.linalg import lstsq

    assert not (OUT / 'verification.json').exists(), 'Refusing to overwrite verification'
    assert not (OUT / 'model_absolute_performance.tsv').exists(), 'Refusing to overwrite diagnostics'
    run = json.loads((OUT / 'model_run.json').read_text())
    audit = json.loads((OUT / 'diagnostic_run.json').read_text())
    spec = json.loads((HERE / 'config/a10_followup_models.json').read_text())
    checks = []

    def require(label, condition):
        assert condition, label
        checks.append(label)

    for kind, record in [('models', run), ('diagnostics', audit)]:
        for name, digest in record['outputs'].items():
            require(kind + ' hash ' + name, sha(OUT / name) == digest)
    require('frozen model specification', sha(HERE / 'config/a10_followup_models.json') == run['spec_sha256'])
    require('executed model script', sha(HERE / 'scripts/08_fit_followup.py') == run['script_sha256'])
    require('original score input', sha(HERE / 'tables/module_scores.tsv') == run['input_hashes']['original_module_scores'])
    require('12 previous instrument values reproduced', len(run['instrument_reproduction']) == 12 and
            all(r['absolute_error'] < .00005 for r in run['instrument_reproduction']))
    pred = pd.read_csv(OUT / 'model_predictions.tsv', sep='\t')
    folds = pd.read_csv(OUT / 'model_fold_errors.tsv', sep='\t')
    summary = pd.read_csv(OUT / 'model_comparisons.tsv', sep='\t')
    fc = pd.read_csv(OUT / 'model_fold_comparisons.tsv', sep='\t')
    require('one prediction per library and evaluation', not pred.duplicated(['library name','setting','scale','evaluation']).any())
    require('prediction settings and size', len(pred) == 13232 and set(pred.setting) == set(spec['sensitivity_settings']))
    absolute = []
    models = spec['models']
    for (setting, scale, evaluation), part in pred.groupby(['setting','scale','evaluation']):
        fold_column = 'unit' if evaluation == 'within_group' else 'plate'
        for fold, group in part.groupby(fold_column):
            ff = folds[(folds.setting == setting) & (folds.scale == scale) & (folds.evaluation == evaluation) & (folds.fold == fold)]
            ss = {}
            for model in models:
                ss[model] = float(np.square(group.outcome - group[model]).sum())
                old = ff[ff.model == model].iloc[0]
                require(f'error arithmetic {setting}/{scale}/{evaluation}/{fold}/{model}', np.isclose(ss[model], old.sse, rtol=1e-8, atol=1e-9))
                if evaluation == 'plate_shift':
                    sst = float(np.square(group.outcome - group.outcome.mean()).sum())
                    absolute.append({'setting':setting,'scale':scale,'plate':fold,'model':model,'wells':len(group),
                                     'sse':ss[model],'rmse':float(np.sqrt(ss[model]/len(group))),
                                     'r2_against_heldout_mean':1-ss[model]/sst})
            for item in spec['comparisons']:
                gain = 1-ss[item['candidate']]/ss[item['reference']]
                old = fc[(fc.setting == setting) & (fc.scale == scale) & (fc.evaluation == evaluation) &
                         (fc.fold == fold) & (fc.comparison == item['name'])].iloc[0]
                require(f'fold gain {setting}/{scale}/{evaluation}/{fold}/{item["name"]}', np.isclose(gain, old.relative_error_reduction, rtol=1e-7, atol=1e-8))
        for item in spec['comparisons']:
            reference = np.square(part.outcome - part[item['reference']]).sum()
            candidate = np.square(part.outcome - part[item['candidate']]).sum()
            old = summary[(summary.setting == setting) & (summary.scale == scale) & (summary.evaluation == evaluation) & (summary.comparison == item['name'])].iloc[0]
            require(f'pooled gain {setting}/{scale}/{evaluation}/{item["name"]}', np.isclose(1-candidate/reference, old.relative_error_reduction, rtol=1e-7, atol=1e-8))
    for name, verdict in run['verdicts'].items():
        x = summary[(summary.comparison == name) & summary.setting.isin(['primary','drop_both'])]
        w, p = x[x.evaluation == 'within_group'], x[x.evaluation == 'plate_shift']
        passed = bool(w.relative_error_reduction.ge(.02).all())
        transferred = bool(passed and p.relative_error_reduction.ge(.02).all() and p.minimum_fold_reduction.gt(0).all())
        require('decision ' + name, verdict == {'within_group_margin_met':passed, 'consistent_plate_shift_gain':transferred})

    # Independent raw-design QR fit, with no call into script 08 or test outcome use.
    source = pd.read_csv(HERE / 'tables/module_scores.tsv', sep='\t').merge(
        pd.read_csv(OUT / 'model_scores.tsv', sep='\t'), on='library name', validate='one_to_one')
    source = source[source.organoids_area_mean_day07.gt(0) & source.organoids_area_mean_day14.gt(0) & source.epi_total.gt(0)].copy()
    source['plate'] = source.unit.str.split('-rep').str[0]
    growth = ['mouse__'+v for vals in spec['programme_blocks'].values() for v in vals]
    qr_errors = []
    for setting in ('primary','drop_both'):
        data = source.copy()
        if setting == 'drop_both':
            data = data[~data.target.isin(['TIGIT','TDTOMATO'])]
        for scale in spec['outcome_scales']:
            day7 = data.organoids_area_mean_day07.to_numpy()
            y = data.organoids_area_mean_day14.to_numpy()
            if scale == 'inherited_log2':
                day7, y = np.log2(day7+1), np.log2(y+1)
            base = np.column_stack([np.ones(len(data)),day7,data.epi_fraction,np.log2(data.epi_total+1)])
            for model, design in [('baseline',base),('baseline+proliferation+remaining_growth',np.column_stack([base,data[growth]]))]:
                for plate in sorted(data.plate.unique()):
                    tr, te = data.plate.ne(plate).to_numpy(), data.plate.eq(plate).to_numpy()
                    beta, _, rank, _ = lstsq(design[tr], y[tr], lapack_driver='gelsy')
                    require(f'QR full rank {setting}/{scale}/{plate}/{model}', rank == design.shape[1])
                    predicted = design[te] @ beta
                    recorded = pred[(pred.setting == setting)&(pred.scale == scale)&(pred.evaluation == 'plate_shift')].set_index('library name').loc[data.loc[te,'library name'],model].to_numpy()
                    error = float(np.max(np.abs(predicted-recorded)))
                    require(f'training-only QR prediction {setting}/{scale}/{plate}/{model}', error < 1e-8)
                    # Changing held-out outcomes leaves training design/response identical.
                    changed_y = y.copy()
                    changed_y[te] = 1e6 + np.arange(int(te.sum()))
                    changed_beta, *_ = lstsq(design[tr], changed_y[tr], lapack_driver='gelsy')
                    require(f'heldout-outcome exclusion {setting}/{scale}/{plate}/{model}', np.array_equal(beta, changed_beta))
                    qr_errors.append(error)
    absolute_frame = pd.DataFrame(absolute)
    absolute_frame.to_csv(OUT / 'model_absolute_performance.tsv',sep='\t',index=False,float_format='%.12g')
    record = {'completed_utc':dt.datetime.now(dt.timezone.utc).isoformat(),'checks_passed':len(checks),
              'checks':checks,'independent_qr_max_prediction_error':max(qr_errors),'independent_qr_fits':len(qr_errors),
              'script_sha256':sha(__file__),'absolute_diagnostic_sha256':sha(OUT/'model_absolute_performance.tsv'),
              'interpretation':'Held-out mean is used only for the diagnostic R2 denominator, never for plate-shift predictions.'}
    (OUT/'verification.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:record[k] for k in ['checks_passed','independent_qr_max_prediction_error','independent_qr_fits']}))
    print(absolute_frame[(absolute_frame.setting=='primary')&(absolute_frame.scale=='inherited_log2')&absolute_frame.model.isin(['baseline','baseline+proliferation+remaining_growth'])][['plate','model','rmse','r2_against_heldout_mean']].to_string(index=False))


if __name__ == '__main__':
    main()
