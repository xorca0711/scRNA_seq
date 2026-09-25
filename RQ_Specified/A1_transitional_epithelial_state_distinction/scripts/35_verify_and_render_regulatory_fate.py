"""Independent source/table checks and scientific figures; preserve earlier runs."""
import csv
import hashlib
import itertools
import json
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
import openpyxl
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parents[1]
ROOT = BASE.parents[1]
TABLE = BASE / 'tables/regulatory_fate'
FIG = BASE / 'figures/regulatory_fate'


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    report = BASE / 'reports/regulatory_fate_verification.json'
    assert not report.exists() and not FIG.exists(), 'Refusing to overwrite verified output'
    checks = {}
    for name, keys in [('regulatory_fate_input_audit.json', ['inputs', 'outputs']),
                       ('regulatory_fate_analysis_run.json', ['source_hashes', 'output_hashes'])]:
        record = json.loads((BASE / 'reports' / name).read_text())
        for key in keys:
            for p, h in record[key].items():
                assert digest(BASE / p) == h, p
        checks[name] = sum(len(record[key]) for key in keys)
    baseline = json.loads((BASE / 'tables/evidence_closure/prior_numerical_sha256.json').read_text())
    baseline.pop((BASE / 'figures/README.md').relative_to(ROOT).as_posix())
    for p, h in baseline.items():
        assert digest(ROOT / p) == h, p
    checks['preserved_preclosure_numerical_artifacts'] = len(baseline)
    for name in ['closure_identity_run.json', 'cd44_closure_run.json']:
        record = json.loads((BASE / 'reports' / name).read_text())
        for p, h in record['output_sha256'].items():
            assert digest(BASE / p) == h, p
        checks['preserved_' + name] = len(record['output_sha256'])

    # Recompute animal-region ratios from original workbook, without analysis helpers.
    w = openpyxl.load_workbook(BASE / 'cache/regulatory_fate/ap1_DC5.xlsx', data_only=True, read_only=True)
    mice = pd.read_csv(TABLE / 'ap1_mice.tsv', sep='\t')
    original = []
    for name in ['wt SeV in situ', 'wt SeV de novo', 'mut SeV in situ', 'mut SeV de novo']:
        frame = pd.DataFrame(list(w[name].values)[2:11], columns=['mouse', 'field', 'total', 'positive', 'pct'])
        for animal, d in frame.groupby('mouse'):
            genotype = name.split()[0]
            region = 'in_situ' if 'in situ' in name else 'de_novo'
            original.append(dict(mouse=f'{genotype}_{animal}', genotype=genotype, region=region,
                                 percent=100 * d.positive.sum()/d.total.sum()))
    source = pd.DataFrame(original)
    match = source.merge(mice, on=['mouse', 'genotype', 'region'], validate='one_to_one')
    assert len(match) == 12 and np.allclose(match.percent, match.count_weighted_percent)
    pivot = source.pivot(index=['genotype', 'mouse'], columns='region', values='percent')
    diff = pivot.de_novo - pivot.in_situ
    observed = diff.loc['mut'].mean() - diff.loc['wt'].mean()
    effects = pd.read_csv(TABLE / 'ap1_region_interaction.tsv', sep='\t')
    assert np.isclose(observed, effects.iloc[0].effect_pp)
    null = []
    a = diff.to_numpy()
    for subset in itertools.combinations(range(6), 3):
        mask = np.isin(np.arange(6), subset)
        null.append(a[mask].mean() - a[~mask].mean())
    assert sum(abs(v) >= abs(observed)-1e-10 for v in null)/20 == effects.iloc[0].p_two_sided == .1
    checks['AP1_source_mouse_region_fractions'] = len(match)
    checks['AP1_exhaustive_allocations'] = len(null)

    # Independent signed-set classification using numeric arrays.
    tw = openpyxl.load_workbook(BASE / 'cache/regulatory_fate/tp53_table1.xlsx', data_only=True, read_only=True)
    shared = pd.DataFrame(list(tw['MDM2-KO Shared'].values)[1:])
    a, b = shared[2].to_numpy(float), shared[8].to_numpy(float)
    signed = pd.read_csv(TABLE / 'tp53_signed_set_summary.tsv', sep='\t')
    expected = dict(both_up=int(((a>0)&(b>0)).sum()), both_down=int(((a<0)&(b<0)).sum()), opposite=int((a*b<0).sum()))
    actual = signed[signed.supplied_set == 'MDM2-KO Shared'].set_index('sign_class').genes.to_dict()
    assert expected == actual == dict(both_up=493, both_down=266, opposite=109)
    checks['TP53_shared_gene_signs'] = len(shared)
    values = pd.read_csv(TABLE / 'tsutsui_source_values.tsv', sep='\t')
    summary = pd.read_csv(TABLE / 'tsutsui_endpoint_summary.tsv', sep='\t')
    recalculated = values.groupby(['panel', 'marker', 'condition']).value.agg(['mean', 'std', 'size']).reset_index()
    joined = summary.merge(recalculated, on=['panel', 'marker', 'condition'], validate='one_to_one', suffixes=('', '_check'))
    assert len(joined) == 65 and np.allclose(joined['mean'], joined.mean_check) and np.allclose(joined.sd, joined['std'])
    assert (joined.n == joined['size']).all()
    checks['Tsutsui_endpoint_recomputations'] = len(joined)

    FIG.mkdir(parents=True)
    plt.rcParams.update({'font.size': 10, 'svg.fonttype': 'none', 'axes.spines.top': False, 'axes.spines.right': False})
    colors = {'wt': '#23779c', 'mut': '#ba5944'}
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.7), gridspec_kw={'width_ratios': [1.3, 1]})
    for genotype in ['wt', 'mut']:
        data = mice[mice.genotype == genotype]
        for animal, d in data.groupby('mouse'):
            vals = d.set_index('region').loc[['in_situ', 'de_novo'], 'count_weighted_percent'].to_numpy()
            ax[0].plot([0, 1], vals, '-o', color=colors[genotype], alpha=.65, linewidth=1.3, markersize=5)
        ax[0].plot([], [], '-o', color=colors[genotype], label='Wild type' if genotype=='wt' else 'AP-1 mutant')
        x = 0 if genotype == 'wt' else 1
        vals = diff.loc[genotype].to_numpy()
        ax[1].scatter(x+np.linspace(-.07,.07,3), vals, color=colors[genotype], s=45)
        ax[1].hlines(vals.mean(), x-.16, x+.16, colors='black', linewidth=2)
    ax[0].set(xticks=[0,1], xticklabels=['In situ\ninjured regions', 'De novo\nintact regions'],
              ylabel='HOPX-positive / GFP-labelled cells (%)', title='Each line connects regions from one mouse')
    ax[0].legend(frameon=False, loc='upper center')
    ax[1].axhline(0, color='.6', linewidth=.8)
    ax[1].set(xticks=[0,1], xticklabels=['Wild type', 'AP-1 mutant'],
              ylabel='De novo − in situ (percentage points)', title=f'Genotype × region contrast: {observed:.2f} pp')
    fig.suptitle('AP-1 loss changes HOPX acquisition in opposite directions by region', fontsize=13)
    fig.text(.05, .01, '3 mice/genotype; fields aggregated within mouse. Exact two-sided permutation p = 0.10 (20 allocations).\nHOPX acquisition is a differentiation readout; mature AT1 function and chromatin mediation were not measured here.', fontsize=9)
    fig.tight_layout(rect=[0,.11,1,.94])
    for ext in ['png','svg']:
        fig.savefig(FIG / ('ap1_region_outcome.'+ext), dpi=180)
    plt.close(fig)

    # Plot culture observations, keeping medium-switch and knockdown experiments separate.
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.8), gridspec_kw={'width_ratios':[1,1.4]})
    reporter = values[(values.panel=='6i') & values.condition.str.startswith('iATCs')]
    for x, condition in enumerate(['iATCs→iAT2s', 'iATCs→iAT1s']):
        nums = reporter[reporter.condition==condition].value.to_numpy()
        ax[0].scatter(x+np.linspace(-.06,.06,len(nums)), nums, color=['#23779c','#b76b24'][x], s=35)
        ax[0].hlines(nums.mean(), x-.13, x+.13, colors='black')
    ax[0].set(yscale='log', xticks=[0,1], xticklabels=['AT2 induction','AT1 induction'],
              ylabel='AGER HiBiT luminescence\n(absolute unit unspecified)', title='Medium-switch capacity: iATC origin')
    genes = ['TP63','MMP7','ITGB6','KRT17','COL1A1']
    for offset, condition, color, label in [(-.13,'siATF3','#ba5944','ATF3 knockdown (n=3)'),(.13,'siHNF1b','#23779c','HNF1B knockdown (n=4)')]:
        for i, gene in enumerate(genes):
            nums = values[(values.panel=='9g') & (values.marker==gene) & (values.condition==condition)].value.to_numpy()
            ax[1].scatter(i+offset+np.linspace(-.025,.025,len(nums)), nums, color=color, s=22, label=label if i==0 else None)
            ax[1].hlines(nums.mean(), i+offset-.08, i+offset+.08, color='black')
    ax[1].axhline(1, color='.6', linestyle='--', linewidth=1)
    ax[1].set(xticks=range(5), xticklabels=genes, ylabel='Expression relative to siCont',
              title='Separate experiment: transitional markers')
    ax[1].legend(frameon=False, fontsize=9, loc='upper center', bbox_to_anchor=(.5,-.12), ncol=2)
    fig.suptitle('Differentiation capacity and regulatory response are distinct pieces of evidence', fontsize=13)
    fig.text(.05,.01,'Dots are author-declared independent culture experiments from one parental iPSC line; black bars are means.\nThese panels do not show that knockdown or p300/CBP inhibition restores the lineage fate of the same cells.',fontsize=9)
    fig.tight_layout(rect=[0,.12,1,.94])
    for ext in ['png','svg']:
        fig.savefig(FIG / ('tsutsui_capacity_and_response.'+ext), dpi=180)
    plt.close(fig)
    report.write_text(json.dumps(dict(utc=datetime.now(timezone.utc).isoformat(), checks=checks,
        code_sha256=digest(Path(__file__)), figures={p.relative_to(BASE).as_posix(): digest(p) for p in FIG.iterdir()},
        visual_review='pending; record inspection separately'), indent=2)+'\n')
    print(json.dumps(checks))


if __name__ == '__main__':
    main()
