"""Freeze the external A5 test definitions before reading expression counts."""
import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    import pandas as pd
    ap = argparse.ArgumentParser()
    ap.add_argument('--source-root', type=Path, required=True)
    ap.add_argument('--data-root', type=Path, required=True)
    args = ap.parse_args()
    config = HERE / 'config/strunz_test_contract.json'
    out = HERE / 'tables/external_test_modules.json'
    record = HERE / 'tables/external_freeze_run.json'
    if any(p.exists() for p in (config, out, record)):
        raise SystemExit('Refusing to overwrite frozen A5 definitions')
    shared = ROOT / 'RQ_Specified/A5_A11_shared_component_contract'
    contract = json.loads((shared / 'config/shared_component.json').read_text())
    guo = args.source_root / 'RQ_Specified/A5_A11_shared_component_contract/sources/guo_2019_supplementary_data_2.xlsx'
    if digest(guo) != contract['source_lists']['D_developmental_transitional']['member_sha256']:
        raise SystemExit('Guo workbook hash differs from frozen source')
    t = pd.read_excel(guo, sheet_name='Drop-seq signature', header=3).iloc[:, :6].dropna(subset=['Gene', 'Group'])
    sets = {k: set(t.loc[t.Group == k, 'Gene']) for k in ['AT1/AT2', 'AT1', 'AT2']}
    assert all(len(g) == 100 for g in sets.values())
    excluded = set(contract['operational_exclusions']['genes'])
    full = sets['AT1/AT2'] - excluded
    identity_free = full - sets['AT1'] - sets['AT2']
    assert len(full) == 99 and len(identity_free) == 57
    gmt = args.data_root / 'raw_data/msigdb/mh.all.v2024.1.Mm.symbols.gmt'
    control_names = ['HALLMARK_P53_PATHWAY', 'HALLMARK_HYPOXIA', 'HALLMARK_INFLAMMATORY_RESPONSE',
                     'HALLMARK_UNFOLDED_PROTEIN_RESPONSE', 'HALLMARK_E2F_TARGETS', 'HALLMARK_G2M_CHECKPOINT']
    controls = {}
    for line in gmt.read_text().splitlines():
        name, _, *genes = line.split('\t')
        if name in control_names:
            controls[name] = set(genes)
    assert set(controls) == set(control_names)
    modules = {'Guo_AT1_AT2_external': full, 'Guo_minus_identity': identity_free,
               'Guo_minus_identity_and_controls': identity_free - set.union(*controls.values()),
               'Guo_AT1': sets['AT1'] - excluded, 'Guo_AT2': sets['AT2'] - excluded, **controls}
    cfg = {
        'schema': 'a5-strunz-external-test/v1', 'date': '2026-09-25',
        'authorization': 'Owner requested proceeding with the review revisions and biological rationale; commit before counts.',
        'cohort': 'GSE141259 high-resolution epithelial experiment',
        'unit': 'mouse', 'day_min': 2, 'day_max': 21, 'exclude_sample_prefix': 'NC',
        'case_label': 'Krt8+ ADI', 'reference_label': 'AT2 activated', 'secondary_reference_label': 'AT2',
        'depth_umi': 500, 'cell_floor': 30, 'unit_floor': 3, 'coverage_floor': 0.7,
        'primary_module': 'Guo_AT1_AT2_external',
        'measurement': 'Expected detection after 500-UMI sampling without replacement; equal-mouse paired difference in percentage points',
        'estimand': 'Equal-mouse mean paired difference, conditional on eligible mice and the sampled day 2-21 schedule',
        'primary_test': 'two-sided one-sample t test of paired differences, alpha 0.05, 95% t interval',
        'meaningful_effect_margin': None,
        'secondary_family': 'Holm across three tests; ineligible tests enter correction as p=1',
        'secondary_tests': [
            {'module': 'Guo_minus_identity', 'reference': 'AT2 activated'},
            {'module': 'Guo_minus_identity_and_controls', 'reference': 'AT2 activated'},
            {'module': 'Guo_AT1_AT2_external', 'reference': 'AT2'}],
        'controls': control_names,
        'descriptive': ['original shared contract modules and Strunz-filtered 51-gene module',
                        'Guo identity axes and six Hallmark axes', 'equal-day mean', 'leave-one-mouse-out means'],
        'interpretation': 'External-source signature recruitment; no untouched-cohort, lineage, mechanism or successful-repair claim',
        'on_primary_failure': 'Stop, do not lower depth/cell floors or switch reference',
    }
    config.parent.mkdir(exist_ok=True)
    config.write_text(json.dumps(cfg, indent=2) + '\n')
    out.write_text(json.dumps({'modules': {k: sorted(v) for k, v in modules.items()},
                               'operational_exclusions': sorted(excluded),
                               'identity_source': 'Guo AT1 and AT2 lists; no Strunz gene ranking used'}, indent=2) + '\n')
    record.write_text(json.dumps({'completed_utc': datetime.now(timezone.utc).isoformat(),
                                 'counts_read': False, 'script_sha256': digest(Path(__file__)),
                                 'sources': {str(p): digest(p) for p in [guo, gmt]},
                                 'outputs': {str(p.relative_to(HERE)): digest(p) for p in [config, out]}}, indent=2) + '\n')
    print(json.dumps({k: len(v) for k, v in modules.items()}))


if __name__ == '__main__':
    main()
