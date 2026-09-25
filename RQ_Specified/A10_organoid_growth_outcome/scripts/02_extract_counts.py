"""Stage 3a: stream the species-split count workbook and cache only what the plan needs.

The workbook has two sheets, mouse and human, each 886 libraries wide with seven
annotation columns. Sheets are about 1.3 GB of XML once decompressed, so rows are
streamed rather than loaded. One pass per sheet computes each library's total counts,
which is the normalization denominator, and retains counts for a declared gene set.

Writes a compact cache plus a tracked run record. Refuses to overwrite.
"""
from __future__ import annotations

import datetime as dt
import hashlib
import json
import re
import sys
import time
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
CACHE = HERE / 'cache'
OUT = HERE / 'tables'
CONFIG = HERE / 'config/a10_outcome_contract.json'
WORKBOOK = CACHE / 'GSE307112_gene_counts.xlsx'
WORKBOOK_SHA = '60736e5c153a4fb92908f19148ebb5ef2b59f7ade97fb1056deee7d0205acc33'
CONTRACT = ROOT / 'RQ_Specified/A5_A11_shared_component_contract/tables'
ES1 = ROOT / 'Research Article/epithelial_state_specificity/modules.json'
ORTHOLOGS = ROOT / 'Research Article/gate2_C3_yu_lee_choi_min_2026/trials/u4_resources/strict_one_to_one_orthologs.csv'
SHEETS = {'mouse': 'xl/worksheets/sheet1.xml', 'human': 'xl/worksheets/sheet2.xml'}
NPZ = CACHE / 'a10_selected_counts.npz'
OUTPUTS = ['extract_run.json', 'library_totals.tsv']
CELL = re.compile(rb'<c r="([A-Z]+)\d+"[^>]*?>(?:<v>([^<]*)</v>|<is><t>([^<]*)</t></is>)</c>')
ROW = re.compile(rb'<row [^>]*?>(.*?)</row>', re.S)
IDENTITY = ['AT1_published_400', 'AT2_published_400']
CONTROLS = ['HALLMARK_P53_PATHWAY', 'HALLMARK_HYPOXIA', 'HALLMARK_INFLAMMATORY_RESPONSE']


def sha256(path: Path) -> str:
    d = hashlib.sha256()
    with path.open('rb') as h:
        for b in iter(lambda: h.read(1 << 24), b''):
            d.update(b)
    return d.hexdigest()


def rel(p: Path) -> str:
    return str(p.relative_to(ROOT)).replace('\\', '/')


def col_index(ref: bytes) -> int:
    n = 0
    for ch in ref:
        n = n * 26 + (ch - 64)
    return n - 1


def needed_genes() -> tuple[dict, dict]:
    import pandas as pd
    frozen = json.loads((CONTRACT / 'frozen_modules.json').read_text(encoding='utf-8'))
    member = pd.read_csv(CONTRACT / 'module_membership.tsv', sep='\t')
    es1 = {m['name']: m['genes'] for m in json.loads(ES1.read_text(encoding='utf-8'))['modules']}
    shared = member[member.module == 'shared_remodelling']
    modules = {k: list(v['genes']) for k, v in frozen['modules'].items()}
    modules['shared_D_and_I'] = sorted(shared[shared.in_D & shared.in_I].gene)
    modules['shared_I_and_L'] = sorted(shared[shared.in_I & shared.in_L].gene)
    for name in IDENTITY + CONTROLS:
        modules[name] = list(es1[name])
    orth = pd.read_csv(ORTHOLOGS)
    m2h = dict(zip(orth.mouse_symbol, orth.human_symbol))
    mouse = sorted({g for genes in modules.values() for g in genes})
    human = sorted({m2h[g] for g in mouse if g in m2h})
    return modules, {'mouse': mouse, 'human': human, 'ortholog_map': m2h}


def scan(zf: zipfile.ZipFile, sheet: str, wanted: set[str]) -> dict:
    keep, totals, header, sym_col = {}, None, None, None
    seen = detected = 0
    start = time.time()
    with zf.open(sheet) as handle:
        buffer = b''
        tail = False
        while True:
            chunk = handle.read(1 << 23)
            if not chunk:
                tail = True
            buffer += chunk
            cut = buffer.rfind(b'</row>')
            block, buffer = (buffer, b'') if tail else ((buffer[:cut + 6], buffer[cut + 6:]) if cut >= 0 else (b'', buffer))
            for m in ROW.finditer(block):
                cells = CELL.findall(m.group(1))
                if header is None:
                    header = {col_index(c[0]): (c[2] or c[1]).decode() for c in cells}
                    sym_col = next(i for i, v in header.items() if v == 'symbol')
                    libs = {i: v for i, v in header.items() if re.match(r'^\d-\d_[A-Z]\d{2}_', v)}
                    totals = {i: 0 for i in libs}
                    order = libs
                    continue
                seen += 1
                symbol = None
                values = {}
                for ref, num, txt in cells:
                    i = col_index(ref)
                    if i == sym_col:
                        symbol = txt.decode()
                    elif i in totals and num:
                        v = int(num) if b'.' not in num else int(float(num))
                        totals[i] += v
                        if v:
                            values[i] = v
                if symbol in wanted:
                    keep.setdefault(symbol, {}).update(values)
                    detected += 1
                if seen % 20000 == 0:
                    print(f'   {sheet}: {seen} rows, {detected} wanted, {time.time()-start:.0f}s', flush=True)
            if tail:
                break
    return {'header': header, 'order': order, 'totals': totals, 'keep': keep, 'rows': seen}


def main() -> None:
    import numpy as np

    OUT.mkdir(exist_ok=True)
    existing = [n for n in OUTPUTS if (OUT / n).exists()] + ([NPZ.name] if NPZ.exists() else [])
    if existing:
        raise SystemExit(f'Refusing to overwrite: {existing}')
    if sha256(WORKBOOK) != WORKBOOK_SHA:
        raise SystemExit('Workbook hash differs from the recorded download; refusing to proceed')

    modules, sets = needed_genes()
    zf = zipfile.ZipFile(WORKBOOK)
    result, arrays, meta = {}, {}, {}
    for species, sheet in SHEETS.items():
        print(f'scanning {species}', flush=True)
        wanted = set(sets[species])
        r = scan(zf, sheet, wanted)
        cols = sorted(r['order'])
        names = [r['order'][i] for i in cols]
        if species == 'mouse':
            library_names = names
        elif names != library_names:
            raise SystemExit('The two sheets do not carry the same libraries in the same order')
        genes = sorted(r['keep'])
        mat = np.zeros((len(genes), len(cols)), dtype=np.int64)
        index = {c: j for j, c in enumerate(cols)}
        for gi, g in enumerate(genes):
            for c, v in r['keep'][g].items():
                mat[gi, index[c]] = v
        arrays[f'{species}_counts'] = mat
        arrays[f'{species}_genes'] = np.array(genes, dtype=object)
        arrays[f'{species}_totals'] = np.array([r['totals'][c] for c in cols], dtype=np.int64)
        meta[species] = {'rows_scanned': r['rows'], 'genes_wanted': len(wanted), 'genes_found': len(genes),
                         'libraries': len(cols)}
        print(f'   {species}: {r["rows"]} rows, {len(genes)} of {len(wanted)} wanted genes found', flush=True)
    arrays['libraries'] = np.array(library_names, dtype=object)
    np.savez_compressed(NPZ, **arrays)

    lines = ['library\tmouse_total_counts\thuman_total_counts']
    for j, lib in enumerate(library_names):
        lines.append(f'{lib}\t{arrays["mouse_totals"][j]}\t{arrays["human_totals"][j]}')
    (OUT / 'library_totals.tsv').write_text('\n'.join(lines) + '\n', encoding='utf-8')

    record = {
        'stage': '3a', 'purpose': 'stream the count workbook and cache the declared gene set only',
        'completed_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
        'script': rel(Path(__file__).resolve()), 'script_sha256': sha256(Path(__file__).resolve()),
        'config_sha256': sha256(CONFIG), 'workbook_sha256': WORKBOOK_SHA,
        'workbook_bytes': WORKBOOK.stat().st_size,
        'tracked_inputs': {rel(p): sha256(p) for p in [CONTRACT / 'frozen_modules.json',
                                                       CONTRACT / 'module_membership.tsv', ES1, ORTHOLOGS]},
        'species_sheets': {'mouse': 'sheet1', 'human': 'sheet2',
                           'note': 'the workbook names its sheets by species, which resolves the stage 1 hold on '
                                   'which assigned genome is which'},
        'per_species': meta, 'libraries': len(library_names),
        'module_sizes': {k: len(v) for k, v in modules.items()},
        'cache': {'path': rel(NPZ), 'sha256': sha256(NPZ), 'bytes': NPZ.stat().st_size},
        'python': sys.version.split()[0],
        'outputs': {n: sha256(OUT / n) for n in OUTPUTS if n != 'extract_run.json'},
    }
    (OUT / 'extract_run.json').write_text(json.dumps(record, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({'per_species': meta, 'libraries': len(library_names)}, indent=2))


if __name__ == '__main__':
    main()
