"""Check whether the planned 0.5/0.7/0.8 assay-coverage gates change any set."""
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
    out.mkdir(exist_ok=True)
    human_sets = PAPER / 'trials/u5_ipf_spec/pathway_genes.tsv'
    jobs = [(c, ROOT / f'analysis/corrections/statistics/cache/{c}_subtypes_counts.csv.gz', human_sets) for c in ['GSE136831', 'GSE135893']]
    jobs.append(('GSE300288', PAPER / 'cache/u3_lineage_annotation/subtype_counts.csv.gz', PAPER / 'trials/u4_resources/pathway_genes.tsv'))
    jobs.append(('GSE308103', PAPER / 'cache/u5_human_full/GSM9237901/unc20_pseudobulks.npz', human_sets))
    # Resolve the first deposited library rather than relying on accession ordering.
    jobs[-1] = (jobs[-1][0], PAPER / 'cache/u5_human_full' / json.loads((PAPER / 'trials/u5_human_full/acquisition_specification.json').read_text())['jobs'][0]['gsm'] / 'unc20_pseudobulks.npz', human_sets)
    rows = []
    for cohort, assay, source in jobs:
        genes = set(np.load(assay)['genes']) if assay.suffix == '.npz' else set(pd.read_csv(assay, usecols=[0]).iloc[:, 0])
        members = pd.read_csv(source, sep='\t')
        for name, g in members.groupby('set'):
            source_genes = set(g.gene)
            fraction = len(source_genes & genes) / len(source_genes)
            rows.append(dict(cohort=cohort, set=name, source_genes=len(source_genes), assayed_genes=len(source_genes & genes),
                             assayed_fraction=fraction, pass_05=fraction >= .5, pass_07=fraction >= .7, pass_08=fraction >= .8))
    result = pd.DataFrame(rows)
    result.to_csv(out / 'pathway_assay_fraction_sensitivity.csv', index=False)
    changed = result[(result.pass_05 != result.pass_07) | (result.pass_08 != result.pass_07)]
    state = dict(status='passed_identical_eligibility' if changed.empty else 'additional_coverage_sensitivity_required',
                 code=code_identity(ROOT, __file__), assessed_set_cohort_rows=len(result), changed_rows=len(changed),
                 minimum_assayed_fraction=float(result.assayed_fraction.min()),
                 conclusion='All declared coverage gates select identical pathway sets; identical inputs/designs give the existing CAMERA results and BH families without rerunning models.' if changed.empty else 'Review the changed set/cohort rows before claiming coverage sensitivity complete.',
                 source_memberships={str(p.relative_to(PAPER)): sha256_file(p) for p in {j[2] for j in jobs}})
    write_json_atomic(out / 'pathway_coverage_validation.json', state)
    print(json.dumps({k: v for k, v in state.items() if k not in ['code', 'source_memberships']}))


if __name__ == '__main__':
    main()
