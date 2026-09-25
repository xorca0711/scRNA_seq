"""Stage 3b: nested models for organoid growth, evaluated by leaving out whole units.

Outcome is day-14 mean organoid area conditional on day-7 mean area. Baseline holds
day-7 size, plate and well composition. Two block comparisons follow in fixed order:
the epithelial programme block, then the fibroblast block. Evaluation leaves out one
whole plate-replicate unit at a time, never sibling wells.

Also runs an assay validation that needs no extra data: genes that are both module
members and perturbation targets should fall in the libraries where they were targeted.

Refuses to overwrite its outputs.
"""
from __future__ import annotations

import datetime as dt
import hashlib
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
CACHE = HERE / 'cache'
OUT = HERE / 'tables'
CONFIG = HERE / 'config/a10_outcome_contract.json'
NPZ = CACHE / 'a10_selected_counts.npz'
CONTRACT = ROOT / 'RQ_Specified/A5_A11_shared_component_contract/tables'
ES1 = ROOT / 'Research Article/epithelial_state_specificity/modules.json'
ORTHOLOGS = ROOT / 'Research Article/gate2_C3_yu_lee_choi_min_2026/trials/u4_resources/strict_one_to_one_orthologs.csv'
OUTPUTS = ['model_comparison.tsv', 'module_scores.tsv', 'knockout_validation.tsv', 'fit_run.json']
IDENTITY = ['AT1_published_400', 'AT2_published_400']
CONTROLS = ['HALLMARK_P53_PATHWAY', 'HALLMARK_HYPOXIA', 'HALLMARK_INFLAMMATORY_RESPONSE']
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


def module_definitions():
    import pandas as pd
    frozen = json.loads((CONTRACT / 'frozen_modules.json').read_text(encoding='utf-8'))
    member = pd.read_csv(CONTRACT / 'module_membership.tsv', sep='\t')
    es1 = {m['name']: m['genes'] for m in json.loads(ES1.read_text(encoding='utf-8'))['modules']}
    shared = member[member.module == 'shared_remodelling']
    mods = {k: list(v['genes']) for k, v in frozen['modules'].items()}
    mods['shared_D_and_I'] = sorted(shared[shared.in_D & shared.in_I].gene)
    mods['shared_I_and_L'] = sorted(shared[shared.in_I & shared.in_L].gene)
    for n in IDENTITY + CONTROLS:
        mods[n] = list(es1[n])
    return mods


def score_modules(counts, genes, totals, modules, mapper=None):
    """Mean log2 CPM over each module's present genes, per library."""
    import numpy as np
    index = {g: i for i, g in enumerate(genes)}
    cpm = np.log2((counts / np.maximum(totals, 1)) * 1e6 + PRIOR)
    out, coverage = {}, {}
    for name, members in modules.items():
        wanted = [mapper[g] for g in members if mapper is None or g in mapper] if mapper else members
        rows = [index[g] for g in wanted if g in index]
        coverage[name] = len(rows) / len(members)
        out[name] = cpm[rows].mean(axis=0) if rows else np.full(counts.shape[1], np.nan)
    return out, coverage


def held_out_r2(design, y, folds):
    """Leave-one-unit-out predictive R squared, fitted by least squares on the training folds."""
    import numpy as np
    pred = np.full(len(y), np.nan)
    for unit in np.unique(folds):
        tr, te = folds != unit, folds == unit
        Xtr, Xte = design[tr], design[te]
        keep = Xtr.std(axis=0) > 0
        keep[0] = True
        beta, *_ = np.linalg.lstsq(Xtr[:, keep], y[tr], rcond=None)
        pred[te] = Xte[:, keep] @ beta
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    return 1 - ss_res / ss_tot, pred


def main() -> None:
    import numpy as np
    import pandas as pd

    OUT.mkdir(exist_ok=True)
    existing = [n for n in OUTPUTS if (OUT / n).exists()]
    if existing:
        raise SystemExit(f'Refusing to overwrite: {existing}')
    config = json.loads(CONFIG.read_text(encoding='utf-8'))
    margin = config['tests']['declared_margin_held_out_r2']
    extract = json.loads((OUT / 'extract_run.json').read_text(encoding='utf-8'))
    if sha256(NPZ) != extract['cache']['sha256']:
        raise SystemExit('The cached counts differ from the extraction record; refusing to fit')

    z = np.load(NPZ, allow_pickle=True)
    libraries = [str(x) for x in z['libraries']]
    modules = module_definitions()
    orth = pd.read_csv(ORTHOLOGS)
    m2h = dict(zip(orth.mouse_symbol, orth.human_symbol))

    epi, epi_cov = score_modules(z['mouse_counts'], [str(g) for g in z['mouse_genes']], z['mouse_totals'], modules)
    fib, fib_cov = score_modules(z['human_counts'], [str(g) for g in z['human_genes']], z['human_totals'], modules, m2h)

    species = pd.read_csv(CACHE / 'GSE307112_xenome_stats.csv.gz')
    imaging = pd.read_csv(CACHE / 'GSE307112_imaging_outputs.csv.gz')
    imaging['well_key'] = imaging.well.str.upper()
    imaging['unit'] = [f'{p}-rep{str(r).split("-")[-1]}' for p, r in zip(imaging.plate, imaging.plate_replicate)]
    wide = imaging.pivot_table(index=['unit', 'well_key'], columns='day',
                               values=['organoids_area_mean', 'organoids_count']).reset_index()
    wide.columns = [f'{a}_{b}' if b else a for a, b in wide.columns]

    lib = pd.DataFrame({'library name': libraries})
    lib['unit'] = [f'plate{n.split("-")[0]}-rep{n.split("-")[1].split("_")[0]}' for n in libraries]
    lib['well_key'] = [re.match(r'^\d-\d_([A-Z]\d{2})_', n).group(1) for n in libraries]
    lib['target'] = [n.split('_', 2)[2] for n in libraries]
    lib['epi_total'] = z['mouse_totals']
    lib['fib_total'] = z['human_totals']
    lib['epi_fraction'] = lib.epi_total / (lib.epi_total + lib.fib_total)
    for name in modules:
        lib[f'epi__{name}'] = epi[name]
        lib[f'fib__{name}'] = fib[name]
    data = lib.merge(wide, on=['unit', 'well_key'], how='left').merge(
        species[['library name', 'xenome_numReadsInput']], on='library name', how='left')
    data.to_csv(OUT / 'module_scores.tsv', sep='\t', index=False)

    d = data.dropna(subset=['organoids_area_mean_day07', 'organoids_area_mean_day14']).copy()
    d = d[(d.organoids_area_mean_day07 > 0) & (d.organoids_area_mean_day14 > 0)]
    d = d[d.epi_total > 0]
    y = np.log2(d.organoids_area_mean_day14.to_numpy() + EPS)
    folds = d.unit.to_numpy()
    plate = pd.get_dummies(d.unit.str.split('-rep').str[0], drop_first=True).to_numpy(float)
    base_cols = [np.ones(len(d)), np.log2(d.organoids_area_mean_day07.to_numpy() + EPS),
                 d.epi_fraction.to_numpy(), np.log2(d.epi_total.to_numpy() + 1.0)]
    baseline = np.column_stack(base_cols + [plate])

    def block(prefix):
        cols = [c for c in d.columns if c.startswith(prefix)]
        m = d[cols].to_numpy(float)
        m = (m - np.nanmean(m, axis=0)) / np.where(np.nanstd(m, axis=0) > 0, np.nanstd(m, axis=0), 1)
        return np.nan_to_num(m), cols

    epi_block, epi_cols = block('epi__')
    fib_block, fib_cols = block('fib__')
    models = {'baseline': baseline,
              'baseline_plus_epithelial': np.column_stack([baseline, epi_block]),
              'baseline_plus_epithelial_plus_fibroblast': np.column_stack([baseline, epi_block, fib_block])}
    rows, r2 = [], {}
    for name, X in models.items():
        value, _ = held_out_r2(X, y, folds)
        r2[name] = value
        rows.append({'model': name, 'features': X.shape[1], 'wells': len(d),
                     'units': int(pd.Series(folds).nunique()), 'held_out_r2': round(value, 4)})
    epi_gain = r2['baseline_plus_epithelial'] - r2['baseline']
    fib_gain = r2['baseline_plus_epithelial_plus_fibroblast'] - r2['baseline_plus_epithelial']
    for label, gain in [('epithelial_increment', epi_gain), ('fibroblast_increment', fib_gain)]:
        rows.append({'model': label, 'features': None, 'wells': len(d), 'units': None,
                     'held_out_r2': round(gain, 4)})
    pd.DataFrame(rows).to_csv(OUT / 'model_comparison.tsv', sep='\t', index=False)

    # Assay validation: module genes that are also perturbation targets should fall when targeted.
    mouse_genes = [str(g) for g in z['mouse_genes']]
    gi = {g: i for i, g in enumerate(mouse_genes)}
    cpm = np.log2((z['mouse_counts'] / np.maximum(z['mouse_totals'], 1)) * 1e6 + PRIOR)
    targets = pd.Series([n.split('_', 2)[2] for n in libraries])
    krows = []
    for gene, idx in sorted(gi.items()):
        hit = targets.str.upper() == gene.upper()
        if int(hit.sum()) < 3:
            continue
        v = cpm[idx]
        krows.append({'gene': gene, 'libraries_targeted': int(hit.sum()),
                      'mean_log2cpm_targeted': round(float(v[hit.to_numpy()].mean()), 3),
                      'mean_log2cpm_other': round(float(v[~hit.to_numpy()].mean()), 3),
                      'difference': round(float(v[hit.to_numpy()].mean() - v[~hit.to_numpy()].mean()), 3)})
    pd.DataFrame(krows).to_csv(OUT / 'knockout_validation.tsv', sep='\t', index=False)
    reduced = sum(1 for r in krows if r['difference'] < 0)

    verdict = ('supported' if epi_gain >= margin else
               'below_declared_margin' if epi_gain < margin else 'inconclusive')
    record = {
        'stage': '3b', 'completed_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
        'script': rel(Path(__file__).resolve()), 'script_sha256': sha256(Path(__file__).resolve()),
        'config_sha256': sha256(CONFIG), 'cache_sha256': extract['cache']['sha256'],
        'declared_margin_held_out_r2': margin,
        'wells_fitted': len(d), 'units': int(pd.Series(folds).nunique()),
        'held_out_r2': {k: round(v, 4) for k, v in r2.items()},
        'epithelial_increment': round(epi_gain, 4), 'fibroblast_increment': round(fib_gain, 4),
        'primary_verdict_epithelial': verdict,
        'module_coverage': {'epithelial': {k: round(v, 3) for k, v in epi_cov.items()},
                            'fibroblast': {k: round(v, 3) for k, v in fib_cov.items()}},
        'knockout_validation': {'genes_testable': len(krows), 'genes_reduced_when_targeted': reduced},
        'evaluation': 'leave one whole plate-replicate unit out; units are batches, not verified independent preparations',
        'interpretation_limit': 'within-screen descriptive association; expression and outcome are both day 14',
        'python': sys.version.split()[0], 'pandas': pd.__version__, 'numpy': np.__version__,
        'outputs': {n: sha256(OUT / n) for n in OUTPUTS if n != 'fit_run.json'},
    }
    (OUT / 'fit_run.json').write_text(json.dumps(record, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({k: record[k] for k in ['wells_fitted', 'units', 'held_out_r2', 'epithelial_increment',
                                             'fibroblast_increment', 'declared_margin_held_out_r2',
                                             'primary_verdict_epithelial', 'knockout_validation']}, indent=2))


if __name__ == '__main__':
    main()
