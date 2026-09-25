"""Power of the pre-registered Kim 2020 test, computed before any Kim score exists.

Uses only the discovery cohort's tracked per-patient differences to set the effect
size. Simulates the exact two-sided Wilcoxon signed-rank test and computes the
paired t test's power analytically, at the discovery effect and at half of it.
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
OUTPUTS = ['power.tsv', 'power_run.json']
DISCOVERY = ROOT / 'Research Article/gate2_C3_yu_lee_choi_min_2026/trials/u6_human_specificity/paired_program_values.csv'
N_PATIENTS = 8
ALPHA = 0.05
SIMULATIONS = 20000
SEED = 20260925
# The discovery run's own primary is the broad type 2 compartment label, as its
# validation record ("primary_broad_contrasts") and report show. The first version
# of this script used the narrower type 2 label; that output is kept under
# tables/superseded_power_label_AT2/ and the narrower label stays as a sensitivity.
LABELS = [('__broad__', 'primary: broad type 2 compartment, the discovery run primary'),
          ('AT2', 'sensitivity: narrower type 2 label')]


def main() -> None:
    import numpy as np
    import pandas as pd
    from scipy import stats

    OUT.mkdir(exist_ok=True)
    existing = [n for n in OUTPUTS if (OUT / n).exists()]
    if existing:
        raise SystemExit(f'Refusing to overwrite: {existing}')
    values = pd.read_csv(DISCOVERY)
    rng = np.random.default_rng(SEED)
    rows, discovery = [], {}
    for discovery_label, role in LABELS:
        d = values[(values.module == 'HPCS_without_ADI_or_operational_markers') & (values.config == 'unc20_pooled')
                   & (values.label == discovery_label) & (values.cell_floor == 50) & (values.case == 'LUAD')
                   & (values.reference == 'normal')].difference.to_numpy()
        if len(d) != 23:
            raise SystemExit(f'Expected 23 discovery patients for {discovery_label}, found {len(d)}')
        dz = d.mean() / d.std(ddof=1)
        discovery[discovery_label] = {'role': role, 'n': int(len(d)), 'mean': round(float(d.mean()), 4),
                                      'sd': round(float(d.std(ddof=1)), 4),
                                      'standardized_effect': round(float(dz), 4), 'positive': int((d > 0).sum())}
        for scenario, effect in [('discovery effect', dz), ('half the discovery effect', dz / 2)]:
            hits = 0
            for _ in range(SIMULATIONS):
                x = rng.normal(effect, 1.0, N_PATIENTS)
                if stats.wilcoxon(x, alternative='two-sided', method='exact').pvalue < ALPHA:
                    hits += 1
            df = N_PATIENTS - 1
            crit = stats.t.ppf(1 - ALPHA / 2, df)
            ncp = effect * np.sqrt(N_PATIENTS)
            t_power = (1 - stats.nct.cdf(crit, df, ncp)) + stats.nct.cdf(-crit, df, ncp)
            rows.append({'discovery_label': discovery_label, 'role': role, 'scenario': scenario,
                         'standardized_effect': round(float(effect), 4), 'n_patients': N_PATIENTS,
                         'wilcoxon_exact_power': round(hits / SIMULATIONS, 4),
                         'paired_t_power': round(float(t_power), 4)})
    pd.DataFrame(rows).to_csv(OUT / 'power.tsv', sep='\t', index=False)

    record = {
        'purpose': 'planning calculation for the pre-registered Kim 2020 test; no Kim data read',
        'completed_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
        'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'discovery_input': str(DISCOVERY.relative_to(ROOT)).replace('\\', '/'),
        'discovery_input_sha256': hashlib.sha256(DISCOVERY.read_bytes()).hexdigest(),
        'discovery': discovery,
        'supersedes': 'tables/superseded_power_label_AT2/, which used only the narrower type 2 label',
        'n_patients': N_PATIENTS, 'alpha': ALPHA, 'simulations': SIMULATIONS, 'seed': SEED,
        'python': sys.version.split()[0], 'scipy': stats.__name__ and __import__('scipy').__version__,
        'outputs': {'power.tsv': hashlib.sha256((OUT / 'power.tsv').read_bytes()).hexdigest()},
    }
    (OUT / 'power_run.json').write_text(json.dumps(record, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({'discovery': record['discovery'], 'power': rows}, indent=2))


if __name__ == '__main__':
    main()
