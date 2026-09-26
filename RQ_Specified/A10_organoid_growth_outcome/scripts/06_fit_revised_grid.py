"""Stage 4b: fit the declared four-cell grid, so baseline and programme set are separable.

Two baselines, pooled and within-unit, crossed with two programme sets, the frozen
repair modules and the growth-nominated Hallmark blocks. Every cell is evaluated by
leaving out one whole unit, with the mandatory sensitivity checks and per-unit
increments the specification requires.

Refuses to overwrite.
"""
from __future__ import annotations

import datetime as dt
import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
CACHE = HERE / 'cache'
OUT = HERE / 'tables'
SPEC = HERE / 'config/a10_revised_specification.json'
OUTPUTS = ['revised_grid.tsv', 'revised_sensitivity.tsv', 'revised_per_unit.tsv', 'revised_run.json']
GMT_RELATIVE = {'mouse': 'raw_data/msigdb/mh.all.v2024.1.Mm.symbols.gmt',
                'human': 'raw_data/msigdb/h.all.v2024.1.Hs.symbols.gmt'}
EPS = 1.0
PRIOR = 1.0


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
    import argparse
    import numpy as np
    import pandas as pd

    parser = argparse.ArgumentParser()
    parser.add_argument('--data-root', type=Path, default=ROOT,
                        help='checkout holding raw_data; defaults to the repository root')
    data_root = parser.parse_args().data_root.resolve()

    OUT.mkdir(exist_ok=True)
    existing = [n for n in OUTPUTS if (OUT / n).exists()]
    if existing:
        raise SystemExit(f'Refusing to overwrite: {existing}')
    spec = json.loads(SPEC.read_text(encoding='utf-8'))
    margin = spec['tests']['declared_margin_held_out_r2']
    blocks = {k.replace('primary_block_', '').replace('secondary_block_', ''): v['sets']
              for k, v in spec['growth_nominated_programmes'].items() if isinstance(v, dict) and 'sets' in v}
    hall_run = json.loads((OUT / 'hallmark_extract_run.json').read_text(encoding='utf-8'))
    hz = np.load(CACHE / 'a10_hallmark_counts.npz', allow_pickle=True)
    if sha256(CACHE / 'a10_hallmark_counts.npz') != hall_run['cache']['sha256']:
        raise SystemExit('Hallmark cache differs from its run record; refusing to fit')

    gmt = {}
    for species, meta in hall_run['gmt'].items():
        path = data_root / GMT_RELATIVE[species]
        if not path.exists():
            raise SystemExit(f'Missing {species} gene set file: {path}')
        if sha256(path) != meta['sha256']:
            raise SystemExit(f'{species} gene set file differs from the extraction record; refusing to fit')
        sets = {}
        for line in path.read_text(encoding='utf-8').splitlines():
            p = line.split('\t')
            if p and p[0] in {s for v in blocks.values() for s in v}:
                sets[p[0]] = [g for g in p[2:] if g]
        gmt[species] = sets

    def hallmark_scores(species):
        genes = [str(g) for g in hz[f'{species}_genes']]
        idx = {g: i for i, g in enumerate(genes)}
        cpm = np.log2((hz[f'{species}_counts'] / np.maximum(hz[f'{species}_totals'], 1)) * 1e6 + PRIOR)
        out = {}
        for block, names in blocks.items():
            for s in names:
                rows = [idx[g] for g in gmt[species][s] if g in idx]
                out[f'{species[:3]}__{block}__{s}'] = cpm[rows].mean(axis=0) if rows else np.full(cpm.shape[1], np.nan)
        return out

    base_data = pd.read_csv(OUT / 'module_scores.tsv', sep='\t')
    libs = [str(x) for x in hz['libraries']]
    hall = pd.DataFrame({'library name': libs})
    for species in ['mouse', 'human']:
        for k, v in hallmark_scores(species).items():
            hall[k] = v
    data = base_data.merge(hall, on='library name', how='inner', validate='one_to_one')

    d = data.dropna(subset=['organoids_area_mean_day07', 'organoids_area_mean_day14']).copy()
    d = d[(d.organoids_area_mean_day07 > 0) & (d.organoids_area_mean_day14 > 0) & (d.epi_total > 0)].reset_index(drop=True)

    frozen = {'epithelial': [c for c in d.columns if c.startswith('epi__')],
              'fibroblast': [c for c in d.columns if c.startswith('fib__')]}
    growth = {'epithelial': [c for c in d.columns if c.startswith('mou__growth__')],
              'fibroblast': [c for c in d.columns if c.startswith('hum__growth__')]}
    extra = {compartment: {block: [c for c in d.columns
                                   if c.startswith(f'{"mou" if compartment == "epithelial" else "hum"}__{block}__')]
                           for block in ('niche_signalling', 'environment')}
             for compartment in ('epithelial', 'fibroblast')}

    def cells(frame, within):
        y = np.log2(frame.organoids_area_mean_day14.to_numpy() + EPS)
        cont = {'day7': np.log2(frame.organoids_area_mean_day07.to_numpy() + EPS),
                'epi_fraction': frame.epi_fraction.to_numpy(float),
                'epi_size': np.log2(frame.epi_total.to_numpy() + 1.0)}
        folds = frame.unit.to_numpy()
        if within:
            def centre(v):
                s = pd.Series(v)
                return (s - s.groupby(pd.Series(folds)).transform('mean')).to_numpy()
            y = centre(y)
            cont = {k: centre(v) for k, v in cont.items()}
            base = np.column_stack([np.ones(len(frame))] + list(cont.values()))
        else:
            plate = pd.get_dummies(pd.Series(folds).str.split('-rep').str[0], drop_first=True).to_numpy(float)
            base = np.column_stack([np.ones(len(frame))] + list(cont.values()) + [plate])

        def block(cols):
            m = frame[cols].to_numpy(float)
            if within:
                m = np.column_stack([(pd.Series(m[:, j]) - pd.Series(m[:, j]).groupby(pd.Series(folds)).transform('mean')).to_numpy()
                                     for j in range(m.shape[1])])
            s = np.nanstd(m, axis=0)
            return np.nan_to_num((m - np.nanmean(m, axis=0)) / np.where(s > 0, s, 1))
        return y, folds, base, block

    rows, per_unit, sens = [], [], []
    for set_name, sets in [('frozen_repair_modules', frozen), ('growth_hallmark', growth)]:
        for baseline_name, within in [('pooled', False), ('within_unit', True)]:
            y, folds, base, block = cells(d, within)
            eb, fb = block(sets['epithelial']), block(sets['fibroblast'])
            r_base, p_base = held_out(base, y, folds)
            r_epi, p_epi = held_out(np.column_stack([base, eb]), y, folds)
            r_fib, p_fib = held_out(np.column_stack([base, eb, fb]), y, folds)
            rows.append({'programme_set': set_name, 'baseline': baseline_name, 'wells': len(d),
                         'units': int(pd.Series(folds).nunique()), 'baseline_r2': round(r_base, 4),
                         'epithelial_increment': round(r_epi - r_base, 4),
                         'fibroblast_increment': round(r_fib - r_epi, 4),
                         'epithelial_clears_margin': bool(r_epi - r_base >= margin),
                         'fibroblast_clears_margin': bool(r_fib - r_epi >= margin)})
            for unit in sorted(set(folds)):
                m = folds == unit
                sst = float(((y[m] - y[m].mean()) ** 2).sum())
                if sst <= 0:
                    continue
                r = lambda p: 1 - float(((y[m] - p[m]) ** 2).sum()) / sst
                per_unit.append({'programme_set': set_name, 'baseline': baseline_name, 'unit': unit,
                                 'wells': int(m.sum()), 'baseline_r2': round(r(p_base), 4),
                                 'epithelial_increment': round(r(p_epi) - r(p_base), 4),
                                 'fibroblast_increment': round(r(p_fib) - r(p_epi), 4)})
            top = d.target.value_counts().idxmax()
            for label, sub in [('primary', d), (f'drop {top}', d[d.target != top]),
                               ('drop TDTOMATO', d[d.target != 'TDTOMATO']),
                               ('drop both', d[(d.target != top) & (d.target != 'TDTOMATO')])]:
                ys, fs, bs, bl = cells(sub.reset_index(drop=True), within)
                e2, f2 = bl(sets['epithelial']), bl(sets['fibroblast'])
                rb, _ = held_out(bs, ys, fs)
                re_, _ = held_out(np.column_stack([bs, e2]), ys, fs)
                rf, _ = held_out(np.column_stack([bs, e2, f2]), ys, fs)
                sens.append({'programme_set': set_name, 'baseline': baseline_name, 'check': label,
                             'wells': len(sub), 'baseline_r2': round(rb, 4),
                             'epithelial_increment': round(re_ - rb, 4), 'fibroblast_increment': round(rf - re_, 4),
                             'epithelial_clears_margin': bool(re_ - rb >= margin),
                             'fibroblast_clears_margin': bool(rf - re_ >= margin)})
            if set_name == 'growth_hallmark':
                for compartment, cols in extra.items():
                    prev = np.column_stack([base, eb, fb])
                    r_prev = r_fib
                    for bname, cc in cols.items():
                        if not cc:
                            continue
                        r_new, _ = held_out(np.column_stack([prev, block(cc)]), y, folds)
                        rows.append({'programme_set': f'secondary:{compartment}:{bname}', 'baseline': baseline_name,
                                     'wells': len(d), 'units': int(pd.Series(folds).nunique()),
                                     'baseline_r2': round(r_prev, 4), 'epithelial_increment': None,
                                     'fibroblast_increment': round(r_new - r_prev, 4),
                                     'epithelial_clears_margin': None,
                                     'fibroblast_clears_margin': bool(r_new - r_prev >= margin)})
    for name, table in [('revised_grid.tsv', rows), ('revised_sensitivity.tsv', sens), ('revised_per_unit.tsv', per_unit)]:
        pd.DataFrame(table).to_csv(OUT / name, sep='\t', index=False)

    def verdict(set_name, baseline_name, which):
        primary = next(r for r in sens if r['programme_set'] == set_name and r['baseline'] == baseline_name
                       and r['check'] == 'primary')
        both = next(r for r in sens if r['programme_set'] == set_name and r['baseline'] == baseline_name
                    and r['check'] == 'drop both')
        key = f'{which}_clears_margin'
        if primary[key] and both[key]:
            return 'supported'
        if primary[key]:
            return 'fragile'
        return 'not_supported'

    verdicts = {f'{s}|{b}|{w}': verdict(s, b, w)
                for s in ('frozen_repair_modules', 'growth_hallmark')
                for b in ('pooled', 'within_unit') for w in ('epithelial', 'fibroblast')}
    record = {
        'stage': '4b', 'completed_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
        'script': rel(Path(__file__).resolve()), 'script_sha256': sha256(Path(__file__).resolve()),
        'spec_sha256': sha256(SPEC), 'hallmark_cache_sha256': hall_run['cache']['sha256'],
        'declared_margin': margin, 'wells': len(d), 'grid': rows, 'verdicts': verdicts,
        'interpretation_limit': 'within-screen descriptive association; the biological unit is unresolved',
        'within_unit_caveat': spec['within_unit_baseline']['honest_limitation'],
        'python': sys.version.split()[0], 'numpy': np.__version__, 'pandas': pd.__version__,
        'outputs': {n: sha256(OUT / n) for n in OUTPUTS if n != 'revised_run.json'},
    }
    (OUT / 'revised_run.json').write_text(json.dumps(record, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({'grid': [r for r in rows if not str(r['programme_set']).startswith('secondary')],
                      'verdicts': verdicts}, indent=2))


if __name__ == '__main__':
    main()
