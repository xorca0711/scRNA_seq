"""Stage 3c: sensitivity checks on the fitted increments. Does not replace the primary.

Each check asks whether an increment survives a specific rival that stage 1 flagged
or that the design invites: fibroblast library size as an unmodelled composition term,
the one target with extreme replication, the control libraries, and dependence on a
single unit. Also tests the assay validation against chance.

Refuses to overwrite its outputs.
"""
from __future__ import annotations

import datetime as dt
import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
OUT = HERE / 'tables'
CONFIG = HERE / 'config/a10_outcome_contract.json'
OUTPUTS = ['sensitivity.tsv', 'per_unit_increments.tsv', 'module_contributions.tsv', 'robustness_run.json']
EPS = 1.0


def sha256(p: Path) -> str:
    d = hashlib.sha256()
    with p.open('rb') as h:
        for b in iter(lambda: h.read(1 << 24), b''):
            d.update(b)
    return d.hexdigest()


def rel(p: Path) -> str:
    return str(p.relative_to(ROOT)).replace('\\', '/')


def held_out(design, y, folds):
    import numpy as np
    pred = np.full(len(y), np.nan)
    for unit in np.unique(folds):
        tr, te = folds != unit, folds == unit
        Xtr = design[tr]
        keep = Xtr.std(axis=0) > 0
        keep[0] = True
        beta, *_ = np.linalg.lstsq(Xtr[:, keep], y[tr], rcond=None)
        pred[te] = design[te][:, keep] @ beta
    ss_res = float(((y - pred) ** 2).sum())
    ss_tot = float(((y - y.mean()) ** 2).sum())
    return 1 - ss_res / ss_tot, pred


def main() -> None:
    import numpy as np
    import pandas as pd
    from scipy import stats

    OUT.mkdir(exist_ok=True)
    existing = [n for n in OUTPUTS if (OUT / n).exists()]
    if existing:
        raise SystemExit(f'Refusing to overwrite: {existing}')
    fit = json.loads((OUT / 'fit_run.json').read_text(encoding='utf-8'))
    margin = fit['declared_margin_held_out_r2']
    data = pd.read_csv(OUT / 'module_scores.tsv', sep='\t')

    d = data.dropna(subset=['organoids_area_mean_day07', 'organoids_area_mean_day14']).copy()
    d = d[(d.organoids_area_mean_day07 > 0) & (d.organoids_area_mean_day14 > 0) & (d.epi_total > 0)].reset_index(drop=True)
    if len(d) != fit['wells_fitted']:
        raise SystemExit(f'Row count {len(d)} differs from the fitted {fit["wells_fitted"]}; refusing to compare')

    epi_cols = [c for c in d.columns if c.startswith('epi__')]
    fib_cols = [c for c in d.columns if c.startswith('fib__')]

    def build(frame, extra_baseline=()):
        y = np.log2(frame.organoids_area_mean_day14.to_numpy() + EPS)
        folds = frame.unit.to_numpy()
        plate = pd.get_dummies(frame.unit.str.split('-rep').str[0], drop_first=True).to_numpy(float)
        cols = [np.ones(len(frame)), np.log2(frame.organoids_area_mean_day07.to_numpy() + EPS),
                frame.epi_fraction.to_numpy(), np.log2(frame.epi_total.to_numpy() + 1.0)]
        for name in extra_baseline:
            cols.append(np.log2(frame[name].to_numpy() + 1.0))
        base = np.column_stack(cols + [plate])

        def std(block):
            m = frame[block].to_numpy(float)
            s = np.nanstd(m, axis=0)
            return np.nan_to_num((m - np.nanmean(m, axis=0)) / np.where(s > 0, s, 1))
        return y, folds, base, std(epi_cols), std(fib_cols)

    rows = []

    def check(label, frame, extra=()):
        y, folds, base, eb, fb = build(frame, extra)
        r_base, _ = held_out(base, y, folds)
        r_epi, _ = held_out(np.column_stack([base, eb]), y, folds)
        r_fib, _ = held_out(np.column_stack([base, eb, fb]), y, folds)
        rows.append({'check': label, 'wells': len(frame), 'units': int(pd.Series(folds).nunique()),
                     'baseline_r2': round(r_base, 4),
                     'epithelial_increment': round(r_epi - r_base, 4),
                     'fibroblast_increment': round(r_fib - r_epi, 4),
                     'epithelial_clears_margin': bool(r_epi - r_base >= margin),
                     'fibroblast_clears_margin': bool(r_fib - r_epi >= margin)})

    check('primary, as fitted', d)
    check('add fibroblast library size to baseline', d, extra=('fib_total',))
    top = d.target.value_counts().idxmax()
    check(f'drop the extreme-replication target ({top})', d[d.target != top])
    check('drop control-like libraries (TDTOMATO)', d[d.target != 'TDTOMATO'])
    check('drop both', d[(d.target != top) & (d.target != 'TDTOMATO')])
    pd.DataFrame(rows).to_csv(OUT / 'sensitivity.tsv', sep='\t', index=False)

    # Per-unit: recompute the increments holding out each unit from the evaluation only.
    y, folds, base, eb, fb = build(d)
    _, p_base = held_out(base, y, folds)
    _, p_epi = held_out(np.column_stack([base, eb]), y, folds)
    _, p_fib = held_out(np.column_stack([base, eb, fb]), y, folds)
    urows = []
    for unit in sorted(set(folds)):
        m = folds == unit
        sst = float(((y[m] - y[m].mean()) ** 2).sum())
        if sst <= 0:
            continue
        r = lambda p: 1 - float(((y[m] - p[m]) ** 2).sum()) / sst
        urows.append({'unit': unit, 'wells': int(m.sum()), 'baseline_r2': round(r(p_base), 4),
                      'epithelial_increment': round(r(p_epi) - r(p_base), 4),
                      'fibroblast_increment': round(r(p_fib) - r(p_epi), 4)})
    pd.DataFrame(urows).to_csv(OUT / 'per_unit_increments.tsv', sep='\t', index=False)
    fib_pos = sum(1 for u in urows if u['fibroblast_increment'] > 0)
    epi_pos = sum(1 for u in urows if u['epithelial_increment'] > 0)

    # Exploratory: which single module, added alone, moves the held-out fit most.
    mrows = []
    r_base, _ = held_out(base, y, folds)
    for col in epi_cols + fib_cols:
        v = d[col].to_numpy(float)
        s = np.nanstd(v)
        x = np.nan_to_num((v - np.nanmean(v)) / (s if s > 0 else 1))
        r, _ = held_out(np.column_stack([base, x]), y, folds)
        mrows.append({'feature': col, 'compartment': 'epithelial' if col.startswith('epi__') else 'fibroblast',
                      'module': col.split('__', 1)[1], 'single_feature_increment': round(r - r_base, 4)})
    mrows.sort(key=lambda r: -r['single_feature_increment'])
    pd.DataFrame(mrows).to_csv(OUT / 'module_contributions.tsv', sep='\t', index=False)

    ko = pd.read_csv(OUT / 'knockout_validation.tsv', sep='\t')
    n, k = len(ko), int((ko.difference < 0).sum())
    binom = float(stats.binomtest(k, n, 0.5, alternative='greater').pvalue)

    record = {
        'stage': '3c', 'completed_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
        'script': rel(Path(__file__).resolve()), 'script_sha256': sha256(Path(__file__).resolve()),
        'declared_margin': margin, 'fit_run_sha256': sha256(OUT / 'fit_run.json'),
        'sensitivity': rows,
        'per_unit': {'units': len(urows), 'epithelial_increment_positive_in': epi_pos,
                     'fibroblast_increment_positive_in': fib_pos},
        'assay_validation': {'genes_testable': n, 'reduced_when_targeted': k,
                             'binomial_one_sided_p': binom,
                             'reading': 'supports that editing reduces the targeted transcript and that the mouse '
                                        'sheet is the perturbed compartment'},
        'top_single_features': mrows[:6],
        'python': sys.version.split()[0], 'scipy': __import__('scipy').__version__,
        'outputs': {n_: sha256(OUT / n_) for n_ in OUTPUTS if n_ != 'robustness_run.json'},
    }
    (OUT / 'robustness_run.json').write_text(json.dumps(record, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({'sensitivity': rows, 'per_unit': record['per_unit'],
                      'assay_validation': record['assay_validation'],
                      'top_single_features': mrows[:6]}, indent=2))


if __name__ == '__main__':
    main()
