"""Expose correlation-model dependence without selecting a new primary analysis."""
from pathlib import Path
import json
import sys

PAPER = Path(__file__).resolve().parents[1]
ROOT = PAPER.parents[1]
sys.path.insert(0, str(ROOT))
from analysis.lib.provenance import sha256_file, code_identity, write_json_atomic


def main():
    import pandas as pd
    import numpy as np
    out = PAPER / 'trials/u6_completion'
    jobs = [(c, PAPER / f'trials/u5_ipf_pathways/{c}_camera.csv', 'FDR_global') for c in ['GSE136831', 'GSE135893']]
    jobs += [('GSE300288', PAPER / 'trials/u4_mouse_niche/camera.csv', 'FDR_global_early'),
             ('GSE308103', PAPER / 'trials/u5_human_niche/camera.csv', 'FDR_global')]
    rows = []
    for cohort, path, qcol in jobs:
        s = pd.read_csv(path); s = s[s.cell_floor.eq(50)]
        if 'config' in s: s = s[s.config.eq('unc20_pooled')]
        a = s[s.correlation_setting.eq('estimated')]
        b = s[s.correlation_setting.eq('fixed001')]
        keys = [k for k in ['compartment', 'comp', 'label', 'case', 'reference', 'set'] if k in s]
        matched = a.merge(b, on=keys, suffixes=('_estimated', '_fixed'), validate='one_to_one')
        assert len(a) == len(b) == len(matched)
        assert np.array_equal(matched.NGenes_estimated, matched.NGenes_fixed)
        rows.append(dict(cohort=cohort, tests=len(a), primary_estimated_q05=int(a[qcol].lt(.05).sum()),
                         fixed001_q05=int(b[qcol].lt(.05).sum()), median_estimated_correlation=float(a.Correlation.median()),
                         estimated_correlation_q25=float(a.Correlation.quantile(.25)), estimated_correlation_q75=float(a.Correlation.quantile(.75)),
                         identical_test_and_gene_count_sets=True, source_sha256=sha256_file(path)))
    table = pd.DataFrame(rows); table.to_csv(out / 'pathway_correlation_sensitivity.csv', index=False)
    write_json_atomic(out / 'pathway_correlation_validation.json', dict(status='passed', code=code_identity(ROOT, __file__),
        cohorts=len(table), identical_test_membership_checked=True,
        interpretation='Primary and fixed-correlation sensitivity results remain separate; a change in model assumptions is not a rescue by assay coverage or proof of absence. No criterion or model was changed.'))
    print(table[['cohort', 'tests', 'primary_estimated_q05', 'fixed001_q05', 'median_estimated_correlation']].to_string(index=False))


if __name__ == '__main__':
    main()
