"""Stage 4a: second streaming pass for the growth-nominated Hallmark genes.

The first pass kept only the frozen repair modules. The declared revision needs
Hallmark blocks in each species' own annotation, so both sheets are streamed again.
Library totals are recomputed and checked against the first pass, which validates that
the two passes read the same workbook the same way.

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
WORKBOOK = CACHE / 'GSE307112_gene_counts.xlsx'
GMT_RELATIVE = {'mouse': 'raw_data/msigdb/mh.all.v2024.1.Mm.symbols.gmt',
                'human': 'raw_data/msigdb/h.all.v2024.1.Hs.symbols.gmt'}
NPZ = CACHE / 'a10_hallmark_counts.npz'
OUTPUTS = ['hallmark_extract_run.json', 'hallmark_coverage.tsv']

sys.path.insert(0, str(HERE / 'scripts'))


def sha256(p: Path) -> str:
    d = hashlib.sha256()
    with p.open('rb') as h:
        for b in iter(lambda: h.read(1 << 24), b''):
            d.update(b)
    return d.hexdigest()


def rel(p: Path) -> str:
    return str(p.relative_to(ROOT)).replace('\\', '/')


def read_gmt(path: Path, wanted: set[str]) -> dict:
    sets = {}
    for line in path.read_text(encoding='utf-8').splitlines():
        parts = line.split('\t')
        if parts and parts[0] in wanted:
            sets[parts[0]] = [g for g in parts[2:] if g]
    missing = wanted - set(sets)
    if missing:
        raise SystemExit(f'Declared sets absent from {path.name}: {sorted(missing)}')
    return sets


def main() -> None:
    import argparse
    import numpy as np
    from importlib.machinery import SourceFileLoader
    mod = SourceFileLoader('pass1', str(HERE / 'scripts/02_extract_counts.py')).load_module()

    parser = argparse.ArgumentParser()
    parser.add_argument('--data-root', type=Path, default=ROOT,
                       help='checkout holding raw_data; defaults to the repository root')
    data_root = parser.parse_args().data_root.resolve()
    gmt_paths = {k: data_root / v for k, v in GMT_RELATIVE.items()}
    for species, path in gmt_paths.items():
        if not path.exists():
            raise SystemExit(f'Missing {species} gene set file: {path}')

    OUT.mkdir(exist_ok=True)
    existing = [n for n in OUTPUTS if (OUT / n).exists()] + ([NPZ.name] if NPZ.exists() else [])
    if existing:
        raise SystemExit(f'Refusing to overwrite: {existing}')
    spec = json.loads(SPEC.read_text(encoding='utf-8'))
    blocks = {k: v['sets'] for k, v in spec['growth_nominated_programmes'].items() if isinstance(v, dict) and 'sets' in v}
    declared = {g for s in blocks.values() for g in s}

    first = json.loads((OUT / 'extract_run.json').read_text(encoding='utf-8'))
    if sha256(WORKBOOK) != first['workbook_sha256']:
        raise SystemExit('Workbook differs from the first pass; refusing to proceed')

    import zipfile
    zf = zipfile.ZipFile(WORKBOOK)
    arrays, meta, cov_rows = {}, {}, []
    prior = np.load(CACHE / 'a10_selected_counts.npz', allow_pickle=True)
    for species, sheet in mod.SHEETS.items():
        sets = read_gmt(gmt_paths[species], declared)
        wanted = {g for genes in sets.values() for g in genes}
        print(f'scanning {species}: {len(wanted)} declared genes', flush=True)
        r = mod.scan(zf, sheet, wanted)
        cols = sorted(r['order'])
        names = [r['order'][i] for i in cols]
        if list(prior['libraries']) != [str(x) for x in names]:
            raise SystemExit(f'{species} library order differs from the first pass')
        totals = np.array([r['totals'][c] for c in cols], dtype=np.int64)
        if not np.array_equal(totals, prior[f'{species}_totals']):
            raise SystemExit(f'{species} library totals differ from the first pass; the passes disagree')
        genes = sorted(r['keep'])
        mat = np.zeros((len(genes), len(cols)), dtype=np.int64)
        index = {c: j for j, c in enumerate(cols)}
        for gi, g in enumerate(genes):
            for c, v in r['keep'][g].items():
                mat[gi, index[c]] = v
        arrays[f'{species}_counts'] = mat
        arrays[f'{species}_genes'] = np.array(genes, dtype=object)
        arrays[f'{species}_totals'] = totals
        found = set(genes)
        for block, members in blocks.items():
            for s in members:
                n = len(sets[s])
                k = sum(1 for g in sets[s] if g in found)
                cov_rows.append({'species': species, 'block': block, 'set': s, 'set_genes': n,
                                 'genes_assayed': k, 'fraction': round(k / n, 4)})
        meta[species] = {'rows_scanned': r['rows'], 'declared_genes': len(wanted), 'genes_found': len(genes),
                         'sets': len(sets)}
        print(f'   {species}: {len(genes)} of {len(wanted)} declared genes found', flush=True)
    arrays['libraries'] = prior['libraries']
    np.savez_compressed(NPZ, **arrays)

    header = list(cov_rows[0].keys())
    (OUT / 'hallmark_coverage.tsv').write_text(
        '\n'.join(['\t'.join(header)] + ['\t'.join(str(r[h]) for h in header) for r in cov_rows]) + '\n',
        encoding='utf-8')

    record = {
        'stage': '4a', 'completed_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
        'script': rel(Path(__file__).resolve()), 'script_sha256': sha256(Path(__file__).resolve()),
        'spec_sha256': sha256(SPEC), 'workbook_sha256': first['workbook_sha256'],
        'gmt': {k: {'path': v.name, 'sha256': sha256(v)} for k, v in gmt_paths.items()},
        'blocks': {k: v for k, v in blocks.items()},
        'per_species': meta,
        'library_totals_match_first_pass': True,
        'cache': {'path': rel(NPZ), 'sha256': sha256(NPZ), 'bytes': NPZ.stat().st_size},
        'outputs': {n: sha256(OUT / n) for n in OUTPUTS if n != 'hallmark_extract_run.json'},
    }
    (OUT / 'hallmark_extract_run.json').write_text(json.dumps(record, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({'per_species': meta, 'totals_match_first_pass': True}, indent=2))


if __name__ == '__main__':
    main()
