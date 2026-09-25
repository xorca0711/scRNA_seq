"""Stage 2: freeze the shared, development-specific and lesion-specific modules.

Applies the partition rule pre-registered in config/shared_component.json to three
independently defined source lists. Reads frozen upstream outputs; never reruns or
edits the scripts that produced them. Refuses to overwrite its own outputs.
"""
from __future__ import annotations

import datetime as dt
import hashlib
import itertools
import json
import platform
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
CONFIG = HERE / 'config/shared_component.json'
OUT = HERE / 'tables'
OUTPUTS = ['frozen_modules.json', 'module_membership.tsv', 'module_overlap.tsv', 'freeze_run.json']

GUO_XLSX = HERE / 'sources/guo_2019_supplementary_data_2.xlsx'
ES1_MODULES = ROOT / 'Research Article/epithelial_state_specificity/modules.json'
U6_SPEC = ROOT / 'Research Article/gate2_C3_yu_lee_choi_min_2026/trials/u6_specificity/module_specification.json'
CONTROLS = ['HALLMARK_P53_PATHWAY', 'HALLMARK_HYPOXIA', 'HALLMARK_INFLAMMATORY_RESPONSE']
IDENTITY = ['AT1_published_400', 'AT2_published_400']


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT)).replace('\\', '/')


def load_guo(expected_sha: str) -> list[str]:
    import pandas as pd
    got = sha256(GUO_XLSX)
    if got != expected_sha:
        raise SystemExit(f'Guo source hash mismatch: {got} != {expected_sha}; refusing to freeze')
    table = pd.read_excel(GUO_XLSX, sheet_name='Drop-seq signature', header=3).iloc[:, :6]
    table = table.dropna(subset=['Gene', 'Group'])
    genes = table.loc[table.Group == 'AT1/AT2', 'Gene'].tolist()
    bad = [g for g in genes if not isinstance(g, str) or not g.strip()]
    if bad:
        raise SystemExit(f'Non-text gene entries in the Guo list, possibly spreadsheet date conversion: {bad}')
    if len(genes) != 100 or len(set(genes)) != 100:
        raise SystemExit(f'Guo AT1/AT2 list has {len(genes)} entries, {len(set(genes))} unique; expected 100 unique')
    return genes


def module_genes(spec: dict, name: str) -> list[str]:
    matches = [m for m in spec['modules'] if m['name'] == name]
    if len(matches) != 1:
        raise SystemExit(f'Expected exactly one module named {name}, found {len(matches)}')
    return list(matches[0]['genes'])


def main() -> None:
    OUT.mkdir(exist_ok=True)
    existing = [n for n in OUTPUTS if (OUT / n).exists()]
    if existing:
        raise SystemExit(f'Refusing to overwrite frozen outputs: {existing}')

    config = json.loads(CONFIG.read_text(encoding='utf-8'))
    exclusions = set(config['operational_exclusions']['genes'])
    es1 = json.loads(ES1_MODULES.read_text(encoding='utf-8'))
    u6 = json.loads(U6_SPEC.read_text(encoding='utf-8'))

    raw = {
        'D': load_guo(config['source_lists']['D_developmental_transitional']['member_sha256']),
        'I': module_genes(es1, 'ADI_published_400'),
        'L': module_genes(u6, 'HPCS_author_top100'),
    }
    for key, expected in [('I', 400), ('L', 100)]:
        if len(raw[key]) != expected or len(set(raw[key])) != expected:
            raise SystemExit(f'Source list {key} has {len(raw[key])} entries; expected {expected} unique')

    removed = {k: sorted(set(v) & exclusions) for k, v in raw.items()}
    lists = {k: set(v) - exclusions for k, v in raw.items()}

    universe = set().union(*lists.values())
    count = {g: sum(g in lists[k] for k in 'DIL') for g in universe}
    modules = {
        'shared_remodelling': sorted(g for g in universe if count[g] >= 2),
        'development_specific': sorted(lists['D'] - lists['I'] - lists['L']),
        'lesion_specific': sorted(lists['L'] - lists['D'] - lists['I']),
        'injury_residual': sorted(lists['I'] - lists['D'] - lists['L']),
    }
    core = sorted(lists['D'] & lists['I'] & lists['L'])

    for a, b in itertools.combinations(modules, 2):
        if set(modules[a]) & set(modules[b]):
            raise SystemExit(f'Partition is not disjoint: {a} and {b} overlap')
    if set().union(*map(set, modules.values())) != universe:
        raise SystemExit('Partition does not cover the union of the three source lists')
    if not set(core) <= set(modules['shared_remodelling']):
        raise SystemExit('Core genes must sit inside the shared module')

    identity = {n: set(module_genes(es1, n)) for n in IDENTITY}
    controls = {n: set(module_genes(es1, n)) for n in CONTROLS}

    rows = []
    for name, genes in modules.items():
        for g in genes:
            rows.append({
                'gene': g, 'module': name,
                'in_D': g in lists['D'], 'in_I': g in lists['I'], 'in_L': g in lists['L'],
                'core_all_three': g in core,
                **{f'in_{n}': g in identity[n] for n in IDENTITY},
                **{f'in_{n}': g in controls[n] for n in CONTROLS},
            })
    header = list(rows[0].keys())
    lines = ['\t'.join(header)] + ['\t'.join(str(r[h]) for h in header) for r in rows]
    (OUT / 'module_membership.tsv').write_text('\n'.join(lines) + '\n', encoding='utf-8')

    sets = {**{k: set(v) for k, v in modules.items()}, **identity, **controls,
            'source_D': lists['D'], 'source_I': lists['I'], 'source_L': lists['L']}
    overlap = ['a\tb\tn_a\tn_b\tn_overlap\tjaccard']
    for a, b in itertools.combinations(sets, 2):
        inter = len(sets[a] & sets[b])
        union = len(sets[a] | sets[b])
        overlap.append(f'{a}\t{b}\t{len(sets[a])}\t{len(sets[b])}\t{inter}\t{inter / union if union else 0:.4f}')
    (OUT / 'module_overlap.tsv').write_text('\n'.join(overlap) + '\n', encoding='utf-8')

    now = dt.datetime.now(dt.timezone.utc).isoformat()
    frozen = {
        'version': 'A5_A11-shared-component-2026-09-25',
        'frozen_utc': now,
        'partition_rule': config['partition_rule'],
        'operational_exclusions': sorted(exclusions),
        'removed_by_exclusion': removed,
        'source_sizes_after_exclusion': {k: len(v) for k, v in lists.items()},
        'sources': {
            'D': {'path': rel(GUO_XLSX), 'sha256': sha256(GUO_XLSX), 'doi': '10.1038/s41467-018-07770-1',
                  'sheet': 'Drop-seq signature', 'group': 'AT1/AT2'},
            'I': {'path': rel(ES1_MODULES), 'sha256': sha256(ES1_MODULES), 'module': 'ADI_published_400',
                  'doi': '10.1038/s41467-020-17358-3'},
            'L': {'path': rel(U6_SPEC), 'sha256': sha256(U6_SPEC), 'module': 'HPCS_author_top100'},
        },
        'core_all_three': core,
        'modules': {k: {'n': len(v), 'genes': v} for k, v in modules.items()},
        'contract_modules': ['shared_remodelling', 'development_specific', 'lesion_specific'],
        'reference_modules': ['injury_residual'],
    }
    (OUT / 'frozen_modules.json').write_text(json.dumps(frozen, indent=2) + '\n', encoding='utf-8')

    import pandas as pd
    record = {
        'stage': 2, 'completed_utc': now,
        'script': rel(Path(__file__).resolve()), 'script_sha256': sha256(Path(__file__).resolve()),
        'config_sha256': sha256(CONFIG),
        'python': sys.version.split()[0], 'platform': platform.platform(), 'pandas': pd.__version__,
        'inputs': {rel(p): sha256(p) for p in [GUO_XLSX, ES1_MODULES, U6_SPEC, CONFIG]},
        'outputs': {n: sha256(OUT / n) for n in OUTPUTS if n != 'freeze_run.json'},
        'expression_data_read': False,
    }
    (OUT / 'freeze_run.json').write_text(json.dumps(record, indent=2) + '\n', encoding='utf-8')

    print(json.dumps({'removed_by_exclusion': removed,
                      'source_sizes_after_exclusion': frozen['source_sizes_after_exclusion'],
                      'module_sizes': {k: len(v) for k, v in modules.items()},
                      'core_all_three': core}, indent=2))


if __name__ == '__main__':
    main()
