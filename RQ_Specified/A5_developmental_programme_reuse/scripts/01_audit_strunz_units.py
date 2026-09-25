"""Audit the Strunz 2020 high-resolution time course as an A5 test cohort.

Reads only metadata: the per-cell label table and the gene list. Never reads counts.
Answers two questions before any A5 scoring is designed: does each sample reach the
unit floor for a within-mouse transitional-versus-type-2 contrast, and are the A5
modules covered by the deposited gene list. Refuses to overwrite its outputs.
"""
from __future__ import annotations

import datetime as dt
import gzip
import hashlib
import json
import platform
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
CACHE = HERE / 'cache'
OUT = HERE / 'tables'
OUTPUTS = ['strunz_sample_units.tsv', 'strunz_module_coverage.tsv', 'audit_run.json']
CONTRACT = ROOT / 'RQ_Specified/A5_A11_shared_component_contract/tables'
FLOOR = 30
UNIT_FLOOR = 3
SOURCES = {
    'GSE141259_HighResolution_cellinfo.csv.gz': {
        'url': 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE141nnn/GSE141259/suppl/GSE141259_HighResolution_cellinfo.csv.gz',
        'bytes': 1424472,
        'sha256': 'fe82e4a9c86720186c03c602148ea52118458e50a0847b869ce375b50a284da8'},
    'GSE141259_HighResolution_genes.txt.gz': {
        'url': 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE141nnn/GSE141259/suppl/GSE141259_HighResolution_genes.txt.gz',
        'bytes': 70000,
        'sha256': 'e7c86817507d191a030350c7f595661e283a0ff10a4f7d5b6fcf52b455892c04'},
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT)).replace('\\', '/')


def main() -> None:
    import pandas as pd

    OUT.mkdir(exist_ok=True)
    existing = [n for n in OUTPUTS if (OUT / n).exists()]
    if existing:
        raise SystemExit(f'Refusing to overwrite: {existing}')
    for name, meta in SOURCES.items():
        path = CACHE / name
        if not path.exists():
            raise SystemExit(f'Missing {name}; fetch it from {meta["url"]} into {rel(CACHE)} first')
        if path.stat().st_size != meta['bytes'] or sha256(path) != meta['sha256']:
            raise SystemExit(f'{name} differs from the audited file; refusing to use it')

    cells = pd.read_csv(CACHE / 'GSE141259_HighResolution_cellinfo.csv.gz', sep='\t')
    if cells.groupby('identifier').sample_id.nunique().max() != 1 or cells.groupby('sample_id').identifier.nunique().max() != 1:
        raise SystemExit('Sample identifiers are not one-to-one; refusing to treat samples as units')
    table = cells.pivot_table(index=['identifier', 'sample_id', 'time_point'], columns='cell_type',
                              values='cell_barcode', aggfunc='count', fill_value=0).reset_index()
    table['day'] = table.time_point.str.extract(r'(\d+)').astype(int)
    for label in ['Krt8+ ADI', 'AT2', 'AT2 activated']:
        if label not in table.columns:
            raise SystemExit(f'Expected author label {label!r} is absent')
    table['meets_floor'] = (table['Krt8+ ADI'] >= FLOOR) & (table['AT2'] >= FLOOR)
    table['control_sample'] = table.sample_id.str.startswith('NC')
    units = table[['identifier', 'sample_id', 'time_point', 'day', 'control_sample', 'Krt8+ ADI', 'AT2',
                   'AT2 activated', 'meets_floor']].sort_values(['day', 'sample_id'])
    units.to_csv(OUT / 'strunz_sample_units.tsv', sep='\t', index=False)

    genes = {line.strip() for line in gzip.open(CACHE / 'GSE141259_HighResolution_genes.txt.gz', 'rt', encoding='utf-8')}
    frozen = json.loads((CONTRACT / 'frozen_modules.json').read_text(encoding='utf-8'))
    member = pd.read_csv(CONTRACT / 'module_membership.tsv', sep='\t')
    shared = member[member.module == 'shared_remodelling']
    development = member[member.module == 'development_specific']
    modules = {
        'development_specific': frozen['modules']['development_specific']['genes'],
        'development_specific_minus_identity': sorted(
            development[~development.in_AT1_published_400 & ~development.in_AT2_published_400].gene),
        'shared_remodelling': frozen['modules']['shared_remodelling']['genes'],
        'shared_D_and_I': sorted(shared[shared.in_D & shared.in_I].gene),
    }
    rows = []
    for name, module in modules.items():
        present = [g for g in module if g in genes]
        rows.append({'module': name, 'n_genes': len(module), 'n_present': len(present),
                     'fraction': round(len(present) / len(module), 4),
                     'gate1': 'pass' if len(present) / len(module) >= 0.7 else 'fail',
                     'missing': ';'.join(g for g in module if g not in genes) or 'none'})
    pd.DataFrame(rows).to_csv(OUT / 'strunz_module_coverage.tsv', sep='\t', index=False)

    injured = units[~units.control_sample]
    record = {
        'purpose': 'A5 cohort audit; metadata only, no counts read',
        'completed_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
        'script': rel(Path(__file__).resolve()), 'script_sha256': sha256(Path(__file__).resolve()),
        'python': sys.version.split()[0], 'platform': platform.platform(), 'pandas': pd.__version__,
        'sources': {name: dict(meta, cached_as=rel(CACHE / name)) for name, meta in SOURCES.items()},
        'contract_inputs': {rel(CONTRACT / n): sha256(CONTRACT / n) for n in ['frozen_modules.json', 'module_membership.tsv']},
        'samples': int(len(units)), 'control_samples': int(units.control_sample.sum()),
        'injured_samples_meeting_floor': int(injured.meets_floor.sum()),
        'days_meeting_floor': sorted(int(d) for d in injured[injured.meets_floor].day.unique()),
        'floor_per_group': FLOOR, 'unit_floor': UNIT_FLOOR,
        'gene_list_size': len(genes),
        'outputs': {n: sha256(OUT / n) for n in OUTPUTS if n != 'audit_run.json'},
        'counts_read': False,
    }
    (OUT / 'audit_run.json').write_text(json.dumps(record, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({k: record[k] for k in ['samples', 'control_samples', 'injured_samples_meeting_floor',
                                             'days_meeting_floor', 'gene_list_size']}, indent=2))


if __name__ == '__main__':
    main()
