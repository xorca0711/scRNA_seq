"""Eligibility gates for the pre-registered Kim 2020 test. No module score is computed.

Reads the deposit's gene names, author cell labels and patient table. The raw count
matrix is streamed only to collect its row names; no count value is kept or used.
Gate 1: assayed source fraction of each module, strict ortholog mapping then
presence in the Kim gene index, denominator the module's own mouse length.
Gate 2: patients with at least 50 cells in both arms, at least three required.
Also fixes the membership of the one pre-registered derived module. Refuses to
overwrite its outputs.

Usage: 01_eligibility_gates.py [--data-root PATH]
"""
from __future__ import annotations

import argparse
import datetime as dt
import gzip
import hashlib
import json
import platform
import re
import sys
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
OUT = HERE / 'tables'
OUTPUTS = ['kim_pairing.tsv', 'kim_cell_counts.tsv', 'gate1_coverage.tsv', 'derived_modules.json', 'gates_run.json']
CONTRACT = ROOT / 'RQ_Specified/A5_A11_shared_component_contract/tables'
ORTHOLOGS = ROOT / 'Research Article/gate2_C3_yu_lee_choi_min_2026/trials/u4_resources/strict_one_to_one_orthologs.csv'
ES1_MODULES = ROOT / 'Research Article/epithelial_state_specificity/modules.json'
CONFIG = HERE / 'config/kim2020_test_contract.json'
KIM = 'raw_data/GSE131907'
NORMAL_LABELS = {'AT2'}
LESION_LABELS = {'tS1', 'tS2', 'tS3'}
CELL_FLOOR = 50
UNIT_FLOOR = 3
THRESHOLD = 0.7
CONTROLS = ['HALLMARK_P53_PATHWAY', 'HALLMARK_HYPOXIA', 'HALLMARK_INFLAMMATORY_RESPONSE']


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


def patient_table(xlsx: Path) -> list[dict]:
    """Read the author sample table from sheet XML; the workbook has a malformed property openpyxl rejects."""
    ns = {'m': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    book = zipfile.ZipFile(xlsx)
    shared = [''.join(t.text or '' for t in si.iter('{%s}t' % ns['m']))
              for si in ET.fromstring(book.read('xl/sharedStrings.xml')).findall('m:si', ns)]
    rows = []
    for row in ET.fromstring(book.read('xl/worksheets/sheet1.xml')).iter('{%s}row' % ns['m']):
        values = {}
        for cell in row.findall('m:c', ns):
            column = re.match(r'[A-Z]+', cell.get('r')).group(0)
            v = cell.find('m:v', ns)
            values[column] = '' if v is None else (shared[int(v.text)] if cell.get('t') == 's' else v.text)
        rows.append(values)
    header = next(r for r in rows if 'Patient id' in r.values())
    keys = {v: k for k, v in header.items()}
    table = []
    for r in rows:
        sample = r.get(keys['Samples'], '')
        if re.fullmatch(r'LUNG_[NT]\d+', sample or ''):
            table.append({'patient': r[keys['Patient id']], 'sample': sample, 'origin': r[keys['Tissue origins']],
                          'histology': r.get(keys['Histology'], ''), 'stage': r.get(keys['Stages'], ''),
                          'egfr': r.get(keys['EGFR'], '')})
    return table


def kim_gene_index(matrix: Path) -> list[str]:
    """Row names of the raw count matrix. Each line is read to its first tab only; values are discarded."""
    names = []
    with gzip.open(matrix, 'rb') as handle:
        handle.readline()
        for line in handle:
            names.append(line.split(b'\t', 1)[0].decode().strip('"'))
    return names


def main() -> None:
    import pandas as pd

    parser = argparse.ArgumentParser()
    parser.add_argument('--data-root', type=Path, default=ROOT)
    data_root = parser.parse_args().data_root.resolve()
    OUT.mkdir(exist_ok=True)
    existing = [n for n in OUTPUTS if (OUT / n).exists()]
    if existing:
        raise SystemExit(f'Refusing to overwrite: {existing}')
    config = json.loads(CONFIG.read_text(encoding='utf-8'))

    kim = data_root / KIM
    xlsx = kim / 'GSE131907_Lung_Cancer_Feature_Summary.xlsx'
    annotation = kim / 'GSE131907_Lung_Cancer_cell_annotation.txt.gz'
    matrix = kim / 'GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz'

    samples = patient_table(xlsx)
    by_patient = {}
    for s in samples:
        by_patient.setdefault(s['patient'], {})[s['origin']] = s
    pairs = []
    for patient, arms in sorted(by_patient.items()):
        if {'nLung', 'tLung'} <= set(arms):
            pairs.append({'patient': patient, 'normal_sample': arms['nLung']['sample'],
                          'tumour_sample': arms['tLung']['sample'], 'histology': arms['tLung']['histology'],
                          'stage': arms['tLung']['stage'], 'egfr': arms['tLung']['egfr']})
    expected = config['cohort']['verified_paired_patients']
    if sorted(p['patient'] for p in pairs) != sorted(expected):
        raise SystemExit('Paired patients in the deposit table differ from the pre-registered list; refusing to continue')
    pd.DataFrame(pairs).to_csv(OUT / 'kim_pairing.tsv', sep='\t', index=False)

    cells = pd.read_csv(annotation, sep='\t', usecols=['Sample', 'Cell_subtype'])
    counts = []
    for p in pairs:
        n_normal = int(((cells.Sample == p['normal_sample']) & cells.Cell_subtype.isin(NORMAL_LABELS)).sum())
        n_lesion = int(((cells.Sample == p['tumour_sample']) & cells.Cell_subtype.isin(LESION_LABELS)).sum())
        counts.append({'patient': p['patient'], 'normal_type2_cells': n_normal, 'lesion_epithelial_cells': n_lesion,
                       'eligible': n_normal >= CELL_FLOOR and n_lesion >= CELL_FLOOR})
    counts = pd.DataFrame(counts)
    counts.to_csv(OUT / 'kim_cell_counts.tsv', sep='\t', index=False)
    eligible_patients = counts[counts.eligible].patient.tolist()

    genes = kim_gene_index(matrix)
    if len(genes) != len(set(genes)):
        raise SystemExit('Duplicate gene names in the Kim matrix; refusing to define presence')
    universe = set(genes)

    orth = pd.read_csv(ORTHOLOGS)
    m2h = dict(zip(orth.mouse_symbol, orth.human_symbol))
    frozen = json.loads((CONTRACT / 'frozen_modules.json').read_text(encoding='utf-8'))
    member = pd.read_csv(CONTRACT / 'module_membership.tsv', sep='\t')
    es1 = {m['name']: m['genes'] for m in json.loads(ES1_MODULES.read_text(encoding='utf-8'))['modules']}
    shared = member[member.module == 'shared_remodelling']
    lesion = member[member.module == 'lesion_specific']
    in_control = lesion[[f'in_{c}' for c in CONTROLS]].any(axis=1)
    stress_excluded = sorted(lesion[~in_control].gene)
    derived = {'lesion_specific_stress_excluded': {
        'definition': 'lesion_specific minus every gene in the three frozen mouse Hallmark control sets',
        'removed': sorted(lesion[in_control].gene), 'genes': stress_excluded, 'n': len(stress_excluded)}}
    (OUT / 'derived_modules.json').write_text(json.dumps(derived, indent=2) + '\n', encoding='utf-8')

    modules = {
        'lesion_specific': ('primary', frozen['modules']['lesion_specific']['genes']),
        'lesion_specific_stress_excluded': ('secondary', stress_excluded),
        'shared_I_and_L': ('secondary', sorted(shared[shared.in_I & shared.in_L].gene)),
        'shared_remodelling': ('secondary baseline and descriptive', frozen['modules']['shared_remodelling']['genes']),
        'AT1_published_400': ('descriptive identity axis', es1['AT1_published_400']),
        'AT2_published_400': ('descriptive identity axis', es1['AT2_published_400']),
        **{c: ('descriptive control axis', es1[c]) for c in CONTROLS},
    }
    rows = []
    for name, (role, module) in modules.items():
        mapped = [m2h[g] for g in module if g in m2h]
        if len(mapped) != len(set(mapped)):
            raise SystemExit(f'Ortholog collision in {name}')
        assayed = [h for h in mapped if h in universe]
        frac = len(assayed) / len(module)
        rows.append({'module': name, 'role': role, 'n_genes': len(module), 'n_orthologs': len(mapped),
                     'n_assayed': len(assayed), 'assayed_fraction': round(frac, 4),
                     'gate1': 'pass' if frac >= THRESHOLD else 'assayed_source_gene_fraction_below_0.7'})
    pd.DataFrame(rows).to_csv(OUT / 'gate1_coverage.tsv', sep='\t', index=False)

    primary = next(r for r in rows if r['module'] == 'lesion_specific')
    record = {
        'purpose': 'eligibility gates for the pre-registered Kim 2020 test; no module score computed',
        'completed_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
        'script': rel(Path(__file__).resolve()), 'script_sha256': sha256(Path(__file__).resolve()),
        'config_sha256': sha256(CONFIG), 'data_root': str(data_root).replace('\\', '/'),
        'python': sys.version.split()[0], 'platform': platform.platform(), 'pandas': pd.__version__,
        'inputs': {rel(p): {'bytes': p.stat().st_size, 'sha256': sha256(p)} for p in [xlsx, annotation, matrix]},
        'tracked_inputs': {rel(p): sha256(p) for p in [CONTRACT / 'frozen_modules.json', CONTRACT / 'module_membership.tsv',
                                                       ORTHOLOGS, ES1_MODULES, CONFIG]},
        'kim_gene_index': {'n_genes': len(genes), 'sorted_names_sha256': hashlib.sha256('\n'.join(sorted(universe)).encode()).hexdigest()},
        'gate1_primary': primary['gate1'], 'gate1_primary_fraction': primary['assayed_fraction'],
        'gate2': {'eligible_patients': eligible_patients, 'n_eligible': len(eligible_patients),
                  'cell_floor': CELL_FLOOR, 'unit_floor': UNIT_FLOOR,
                  'verdict': 'pass' if len(eligible_patients) >= UNIT_FLOOR else 'fewer_than_3_complete_units'},
        'counts_used': False,
        'outputs': {n: sha256(OUT / n) for n in OUTPUTS if n != 'gates_run.json'},
    }
    (OUT / 'gates_run.json').write_text(json.dumps(record, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({'gate1': rows, 'gate2': record['gate2'], 'kim_genes': len(genes)}, indent=2))


if __name__ == '__main__':
    main()
