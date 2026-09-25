"""Stage 3: report module eligibility before any expression score is computed.

Reads gene names only from each target dataset, never counts. Two gates:
gate 1, the assayed source fraction (strict ortholog mapping for the human arm,
then presence in the target gene index; denominator is the module's own length);
gate 2, the complete-unit floor, cited from tracked records rather than recomputed.

Fails closed: before reporting any new module, the script must reproduce the
coverage values of the two completed runs exactly. If either check fails it
writes the failed check and stops, without a coverage table.

Usage: 02_report_coverage.py [--data-root PATH]
  --data-root  directory holding raw_data/ and the ignored paper caches;
               defaults to the repository root.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import platform
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
OUT = HERE / 'tables'
OUTPUTS = ['instrument_check.tsv', 'gene_universes.tsv', 'coverage_gate.tsv',
           'coverage_gene_detail.tsv', 'unit_floor.tsv', 'coverage_run.json']
THRESHOLD = 0.7
UNIT_FLOOR = 3

ES1 = ROOT / 'Research Article/epithelial_state_specificity'
PAPER = ROOT / 'Research Article/gate2_C3_yu_lee_choi_min_2026'
ORTHOLOGS = PAPER / 'trials/u4_resources/strict_one_to_one_orthologs.csv'
U6_SPEC = PAPER / 'trials/u6_specificity/module_specification.json'
U6_SUMMARY = PAPER / 'trials/u6_human_specificity/paired_program_summary.csv'
FROZEN = HERE / 'tables/frozen_modules.json'
MEMBERSHIP = HERE / 'tables/module_membership.tsv'
FREEZE_RUN = HERE / 'tables/freeze_run.json'

MOUSE_CHECK = ['ADI_published_400', 'ADI_published_holdout', 'AT1_published_400', 'AT1_published_holdout',
               'AT2_published_400', 'AT2_published_holdout', 'HALLMARK_HYPOXIA',
               'HALLMARK_INFLAMMATORY_RESPONSE', 'HALLMARK_P53_PATHWAY']
HUMAN_CHECK = ['ADI_published_holdout', 'HPCS_author_top100', 'HPCS_without_ADI_or_operational_markers',
               'AT1_published_holdout', 'AT2_published_holdout', 'DATP_PATS_holdout_panel']
RELEVANCE = {
    ('mouse', 'shared_remodelling'): 'A5 primary', ('mouse', 'development_specific'): 'A5 primary',
    ('mouse', 'shared_D_and_I'): 'A5 declared secondary',
    ('human', 'shared_remodelling'): 'A11 primary', ('human', 'lesion_specific'): 'A11 primary',
    ('human', 'shared_I_and_L'): 'A11 declared secondary',
}


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1 << 24), b''):
            digest.update(block)
    return digest.hexdigest()


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT)).replace('\\', '/')
    except ValueError:
        return str(path).replace('\\', '/')


POSTQC = 'Research Article/gate1_01_niethamer_2025/GSE262927/processed/postQC.h5ad'


def postqc_gene_names(path: Path) -> tuple[str, ...]:
    """Variable names of the post-QC object the precedent scored, not the raw deposit."""
    import anndata as ad
    import h5py
    with h5py.File(path, 'r') as h:
        var = ad.io.read_elem(h['var'])
    return tuple(var.index.astype(str))


def verify_postqc(path: Path) -> str:
    """The object must be byte-identical to the input the precedent recorded."""
    record = json.loads((ES1 / 'results/run_record.json').read_text(encoding='utf-8'))
    entries = [e for e in record['inputs'] if str(e['path']).replace('\\', '/').endswith('GSE262927/processed/postQC.h5ad')]
    if len(entries) != 1:
        raise SystemExit(f'Expected one recorded post-QC input in the precedent run record, found {len(entries)}')
    if path.stat().st_size != entries[0]['bytes']:
        raise SystemExit('Post-QC object size differs from the precedent input; refusing to use it')
    got = sha256(path)
    if got != entries[0]['sha256']:
        raise SystemExit(f'Post-QC object hash {got} differs from the precedent input; refusing to use it')
    return got


def tenx_gene_names(path: Path) -> tuple[str, ...]:
    import h5py
    with h5py.File(path, 'r') as h:
        names = [x.decode() for x in h['matrix/features/name'][:]]
        types = [x.decode() for x in h['matrix/features/feature_type'][:]]
    return tuple(n for n, t in zip(names, types) if t == 'Gene Expression')


def single_universe(label: str, files: list[Path], reader) -> tuple[frozenset, dict]:
    if not files:
        raise SystemExit(f'No files found for {label}; refusing to report')
    indexes = {rel(p): reader(p) for p in files}
    distinct = set(indexes.values())
    if len(distinct) != 1:
        raise SystemExit(f'{label}: {len(distinct)} different gene indexes across {len(files)} files; refusing to pick one')
    names = next(iter(distinct))
    digest = sha256_bytes('\n'.join(sorted(set(names))).encode())
    return frozenset(names), {'dataset': label, 'n_files': len(files), 'n_features': len(names),
                              'n_unique_names': len(set(names)), 'sorted_names_sha256': digest}


def write_tsv(path: Path, rows: list[dict]) -> None:
    header = list(rows[0].keys())
    lines = ['\t'.join(header)] + ['\t'.join('' if r[h] is None else str(r[h]) for h in header) for r in rows]
    path.write_text('\n'.join(lines) + '\n', encoding='utf-8')


def main() -> None:
    import h5py
    import pandas as pd

    parser = argparse.ArgumentParser()
    parser.add_argument('--data-root', type=Path, default=ROOT)
    data_root = parser.parse_args().data_root.resolve()

    OUT.mkdir(exist_ok=True)
    existing = [n for n in OUTPUTS if (OUT / n).exists()]
    if existing:
        raise SystemExit(f'Refusing to overwrite: {existing}')

    freeze = json.loads(FREEZE_RUN.read_text(encoding='utf-8'))
    for name in ['frozen_modules.json', 'module_membership.tsv']:
        if sha256(OUT / name) != freeze['outputs'][name]:
            raise SystemExit(f'{name} no longer matches the stage 2 run record; refusing to report')

    raw = data_root / 'raw_data'
    human_cache = data_root / 'Research Article/gate2_C3_yu_lee_choi_min_2026/cache/u5_human_full'
    universes, universe_rows = {}, []
    for label, files in [
        ('GSE247130', sorted((raw / 'GSE247130').glob('GSE247130_Aggregate*_filtered_feature_bc_matrix.h5'))),
        ('GSE310539', sorted((raw / 'GSE310539').glob('*_filtered_feature_bc_matrix.h5'))),
    ]:
        universes[label], row = single_universe(label, files, tenx_gene_names)
        universe_rows.append(dict(row, arm='mouse', universe_source='raw 10x Gene Expression feature names'))
    postqc = data_root / POSTQC
    postqc_sha = verify_postqc(postqc)
    universes['GSE262927'], row = single_universe('GSE262927', [postqc], postqc_gene_names)
    universe_rows.append(dict(row, arm='mouse',
                              universe_source=f'post-QC object variable names, sha256 {postqc_sha}, as scored by the precedent'))

    def human_reader(p: Path) -> tuple[str, ...]:
        with h5py.File(p, 'r') as h:
            return tuple(x.decode() if isinstance(x, bytes) else str(x) for x in h['genes'][:])
    universes['GSE308103'], row = single_universe(
        'GSE308103', sorted(human_cache.glob('*/raw_counts.csc.h5')), human_reader)
    universe_rows.append(dict(row, arm='human', universe_source='processed per-library gene index shared by all libraries'))

    orth = pd.read_csv(ORTHOLOGS)
    m2h = dict(zip(orth.mouse_symbol, orth.human_symbol))
    es1_modules = {m['name']: m['genes'] for m in json.loads((ES1 / 'modules.json').read_text(encoding='utf-8'))['modules']}
    u6_modules = {m['name']: m['genes'] for m in json.loads(U6_SPEC.read_text(encoding='utf-8'))['modules']}

    def mouse_count(genes, universe):
        return sum(g in universe for g in genes)

    def human_counts(genes, universe):
        mapped = [m2h[g] for g in genes if g in m2h]
        if len(mapped) != len(set(mapped)):
            raise SystemExit('Strict one-to-one mapping produced a collision; refusing to report')
        return len(mapped), sum(h in universe for h in mapped)

    checks, ok = [], True
    scores = pd.read_csv(ES1 / 'results/module_scores.csv')
    expected = scores[['cohort', 'module', 'n_genes', 'total_genes']].drop_duplicates()
    for _, r in expected[expected.module.isin(MOUSE_CHECK)].iterrows():
        genes = es1_modules[r.module]
        observed = mouse_count(genes, universes[r.cohort])
        match = observed == int(r.n_genes) and len(genes) == int(r.total_genes)
        ok &= match
        checks.append({'arm': 'mouse', 'dataset': r.cohort, 'module': r.module, 'expected': int(r.n_genes),
                       'observed': observed, 'denominator': len(genes), 'match': match,
                       'expected_source': rel(ES1 / 'results/module_scores.csv')})
    summary = pd.read_csv(U6_SUMMARY)
    human_expected = summary[['module', 'source_coverage']].drop_duplicates()
    for _, r in human_expected[human_expected.module.isin(HUMAN_CHECK)].iterrows():
        genes = u6_modules[r.module]
        _, assayed = human_counts(genes, universes['GSE308103'])
        observed = assayed / len(genes)
        match = abs(observed - float(r.source_coverage)) < 1e-6
        ok &= match
        checks.append({'arm': 'human', 'dataset': 'GSE308103', 'module': r.module,
                       'expected': round(float(r.source_coverage), 6), 'observed': round(observed, 6),
                       'denominator': len(genes), 'match': match, 'expected_source': rel(U6_SUMMARY)})
    write_tsv(OUT / 'instrument_check.tsv', checks)
    write_tsv(OUT / 'gene_universes.tsv', universe_rows)
    if not ok or len(checks) < len(MOUSE_CHECK) * 3 + len(HUMAN_CHECK):
        failed = [c for c in checks if not c['match']]
        (OUT / 'coverage_run.json').write_text(json.dumps(
            {'status': 'refused_instrument_check_failed', 'failed': failed, 'n_checks': len(checks)}, indent=2) + '\n',
            encoding='utf-8')
        raise SystemExit(f'Instrument check failed or incomplete ({len(failed)} failed, {len(checks)} run); '
                         'no coverage reported')

    frozen = json.loads(FROZEN.read_text(encoding='utf-8'))
    member = pd.read_csv(MEMBERSHIP, sep='\t')
    shared = member[member.module == 'shared_remodelling']
    modules = {name: list(frozen['modules'][name]['genes']) for name in
               ['shared_remodelling', 'development_specific', 'lesion_specific', 'injury_residual']}
    modules['shared_D_and_I'] = sorted(shared[shared.in_D & shared.in_I].gene)
    modules['shared_I_and_L'] = sorted(shared[shared.in_I & shared.in_L].gene)
    modules['shared_D_and_L'] = sorted(shared[shared.in_D & shared.in_L].gene)
    for key in 'DIL':
        modules[f'source_{key}'] = sorted(member[member[f'in_{key}']].gene)
    for name in ['AT1_published_400', 'AT2_published_400', 'HALLMARK_P53_PATHWAY', 'HALLMARK_HYPOXIA',
                 'HALLMARK_INFLAMMATORY_RESPONSE']:
        modules[name] = list(es1_modules[name])
    role = {'shared_remodelling': 'contract', 'development_specific': 'contract', 'lesion_specific': 'contract',
            'injury_residual': 'reference', 'shared_D_and_I': 'declared secondary',
            'shared_I_and_L': 'declared secondary', 'shared_D_and_L': 'declared secondary',
            'source_D': 'source list', 'source_I': 'source list', 'source_L': 'source list',
            'AT1_published_400': 'identity axis', 'AT2_published_400': 'identity axis',
            'HALLMARK_P53_PATHWAY': 'control axis', 'HALLMARK_HYPOXIA': 'control axis',
            'HALLMARK_INFLAMMATORY_RESPONSE': 'control axis'}

    rows = []
    for name, genes in modules.items():
        n = len(genes)
        for arm, dataset in [('mouse', 'GSE247130'), ('mouse', 'GSE310539'), ('mouse', 'GSE262927'),
                             ('human', 'GSE308103')]:
            base = {'module': name, 'role': role[name], 'relevant_to': RELEVANCE.get((arm, name), ''),
                    'arm': arm, 'dataset': dataset, 'n_genes': n}
            if n == 0:
                rows.append(dict(base, n_orthologs=None, ortholog_fraction=None, n_assayed=0,
                                 assayed_fraction=None, gate1='empty_module'))
                continue
            if arm == 'mouse':
                assayed = mouse_count(genes, universes[dataset])
                n_orth, orth_frac = None, None
            else:
                n_orth, assayed = human_counts(genes, universes[dataset])
                orth_frac = round(n_orth / n, 4)
            frac = assayed / n
            rows.append(dict(base, n_orthologs=n_orth, ortholog_fraction=orth_frac, n_assayed=assayed,
                             assayed_fraction=round(frac, 4),
                             gate1='pass' if frac >= THRESHOLD else 'assayed_source_gene_fraction_below_0.7'))
    write_tsv(OUT / 'coverage_gate.tsv', rows)

    detail = []
    for name in ['shared_remodelling', 'development_specific', 'lesion_specific']:
        for g in modules[name]:
            h = m2h.get(g)
            detail.append({'module': name, 'mouse_symbol': g,
                           **{f'in_{d}': g in universes[d] for d in ['GSE247130', 'GSE310539', 'GSE262927']},
                           'human_ortholog': h or '', 'in_GSE308103': bool(h) and h in universes['GSE308103']})
    write_tsv(OUT / 'coverage_gene_detail.tsv', detail)

    primary = summary[(summary.config == 'unc20_pooled') & (summary.label == 'AT2') & (summary.cell_floor == 50)]
    per_contrast = primary[['case', 'reference', 'n_patients']].drop_duplicates()
    if per_contrast.duplicated(['case', 'reference']).any():
        raise SystemExit('Primary configuration gives more than one patient count per contrast; refusing to cite')
    units = [{'arm': 'human', 'question': 'A11', 'unit': 'patient',
              'contrast': f'{r.case} minus {r.reference}', 'n_complete_units': int(r.n_patients),
              'floor': UNIT_FLOOR, 'gate2': 'pass' if r.n_patients >= UNIT_FLOOR else 'fewer_than_3_complete_units',
              'source': f'{rel(U6_SUMMARY)}; config unc20_pooled, label AT2, cell floor 50'}
             for r in per_contrast.itertuples()]
    es1_summary = (ES1 / 'results/SUMMARY.md').read_text(encoding='utf-8')
    found = re.search(r'External sample check: (\d+) animals scored; (\d+) meet both 30-cell group floors', es1_summary)
    es1_run = json.loads((ES1 / 'results/run_record.json').read_text(encoding='utf-8'))
    if not found or int(found.group(1)) != es1_run['external']['animal_count']:
        raise SystemExit('Cannot confirm the external animal count in the tracked specificity record; refusing to cite')
    if 'no neonatal arm' not in es1_run['external']['caution']:
        raise SystemExit('The tracked external caution no longer states the absence of a neonatal arm')
    retention = pd.read_csv(ES1 / 'results/retention.csv')
    multiome = retention[retention.cohort.isin(['GSE247130', 'GSE310539'])]
    pooled_only = bool((multiome.unit_type == 'pooled_library').all())
    units.append({'arm': 'mouse', 'question': 'A5', 'unit': 'animal',
                  'contrast': 'adult injury, external animals meeting both group floors',
                  'n_complete_units': int(found.group(2)), 'floor': UNIT_FLOOR,
                  'gate2': 'pass' if int(found.group(2)) >= UNIT_FLOOR else 'fewer_than_3_complete_units',
                  'source': f"{rel(ES1 / 'results/SUMMARY.md')}; {found.group(1)} animals scored"})
    units.append({'arm': 'mouse', 'question': 'A5', 'unit': 'animal',
                  'contrast': 'neonatal development, independent animals',
                  'n_complete_units': 0, 'floor': UNIT_FLOOR, 'gate2': 'fewer_than_3_complete_units',
                  'source': ('external cohort caution states no neonatal arm; every multiome neonatal unit is a '
                             f'pooled library ({pooled_only}) in {rel(ES1 / "results/retention.csv")}')})
    write_tsv(OUT / 'unit_floor.tsv', units)

    record = {
        'stage': 3, 'status': 'reported', 'completed_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
        'script': rel(Path(__file__).resolve()), 'script_sha256': sha256(Path(__file__).resolve()),
        'data_root': str(data_root).replace('\\', '/'),
        'python': sys.version.split()[0], 'platform': platform.platform(), 'pandas': pd.__version__,
        'h5py': h5py.__version__,
        'tracked_inputs': {rel(p): sha256(p) for p in [FROZEN, MEMBERSHIP, FREEZE_RUN, ORTHOLOGS, U6_SPEC, U6_SUMMARY,
                                                       ES1 / 'modules.json', ES1 / 'results/module_scores.csv',
                                                       ES1 / 'results/SUMMARY.md', ES1 / 'results/run_record.json',
                                                       ES1 / 'results/retention.csv']},
        'gene_universes': universe_rows,
        'instrument_checks': {'n': len(checks), 'all_match': ok},
        'expression_counts_read': False,
        'outputs': {n: sha256(OUT / n) for n in OUTPUTS if n != 'coverage_run.json'},
    }
    (OUT / 'coverage_run.json').write_text(json.dumps(record, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({'instrument_checks': record['instrument_checks'], 'universes': universe_rows}, indent=2))


if __name__ == '__main__':
    main()
