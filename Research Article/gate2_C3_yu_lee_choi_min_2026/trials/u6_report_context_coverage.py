"""Show verified designs and primary biological-unit coverage without pooling contexts."""
from pathlib import Path
import json

PAPER = Path(__file__).resolve().parents[1]
ROOT = PAPER.parents[1]


def main():
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    out = PAPER / 'trials/u6_completion'
    out.mkdir(exist_ok=True)
    rows = []
    for cohort in ['GSE136831', 'GSE135893']:
        e = pd.read_csv(PAPER / f'trials/u5_ipf_pathways/{cohort}_eligibility.csv')
        for row in e[e.cell_floor.eq(50) & e.label.eq('__broad__')].itertuples():
            rows.append(dict(cohort=cohort, comp=row.compartment, case='IPF', reference='control', n_case=row.n_IPF, n_reference=row.n_control, eligible=row.eligible))
    e = pd.read_csv(PAPER / 'trials/u4_mouse_niche/eligibility.csv')
    rename = dict(alveolar='AT2', fibroblast='fibroblasts', myeloid='macrophages')
    for row in e[e.cell_floor.eq(50) & e.label.eq('__broad__')].itertuples():
        rows.append(dict(cohort='GSE300288', comp=rename[row.compartment], case='anti-IL1B', reference='IgG', n_case=row.n_antiIL1B, n_reference=row.n_control, eligible=row.eligible))
    independent = pd.DataFrame(rows)
    independent.to_csv(out / 'primary_broad_independent_arm_coverage.csv', index=False)
    all_pairs = pd.read_csv(PAPER / 'trials/u5_human_niche/paired_eligibility.csv')
    paired = all_pairs[all_pairs.config.eq('unc20_pooled') & all_pairs.cell_floor.eq(50) & all_pairs.label.eq('__broad__')]
    paired.to_csv(out / 'primary_human_broad_pair_coverage.csv', index=False)
    assert not paired.duplicated(['comp', 'case', 'reference']).any()
    design = pd.DataFrame([
        ['Direct blockade', 'GSE300288', '7 early mice (4 IgG / 3 anti-IL1B)', 'Intervention; KAC identity unavailable'],
        ['Fibrotic niche', 'GSE136831 / GSE135893', '60 / 22 donors in full source scans', 'Observational; subtype coverage differs'],
        ['Human lesions', 'GSE308103', '75 libraries / 23 patient labels', 'Paired histologies; reference-compatible cells'],
        ['Human spatial companion', 'GSE307534', '56 sections / 25 patient labels', 'Shares 23 RNA patient labels; no ROI labels'],
        ['Pathological post-viral repair', 'GSE267226 / GSE267228', '5 human donors / 4 mice', 'Whole-sample RNA; no deposited coordinates'],
        ['Repair / HPCS specificity', 'GSE247130 / 310539 / 277777', 'Pooled wells / source identifiers', 'Descriptive; not independent cell replicates'],
    ], columns=['Context', 'Deposit', 'Observed sampling', 'Interpretation boundary'])
    design.to_csv(out / 'context_design_map.csv', index=False)
    pal = json.loads((ROOT / 'analysis/config/palette.json').read_text())
    cmap = LinearSegmentedColormap.from_list('coverage', [pal['surface'], pal['sequential_ramp'][2]])
    fig = plt.figure(figsize=(16, 10.5), facecolor=pal['surface'])
    gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1], width_ratios=[1, 1.7])
    ax = fig.add_subplot(gs[0, :]); ax.axis('off')
    table = ax.table(cellText=design.values, colLabels=design.columns, cellLoc='left', colLoc='left',
                     colWidths=[.19, .23, .26, .32], bbox=[0, 0, 1, 1])
    table.auto_set_font_size(False); table.set_fontsize(9)
    for (i, j), cell in table.get_celld().items():
        cell.set_facecolor(pal['deemph'] if i == 0 else pal['surface'])
        cell.set_edgecolor(pal['grid']); cell.set_text_props(color=pal['ink'])
    comps = ['AT2', 'fibroblasts', 'macrophages']
    ax = fig.add_subplot(gs[1, 0])
    cohorts = ['GSE136831', 'GSE135893', 'GSE300288']
    vals = np.zeros((3, 3)); labels = {}
    for i, cohort in enumerate(cohorts):
        for j, comp in enumerate(comps):
            row = independent[independent.cohort.eq(cohort) & independent.comp.eq(comp)].iloc[0]
            vals[i, j] = min(row.n_case, row.n_reference)
            labels[(i, j)] = f'{row.n_case} / {row.n_reference}'
    ax.imshow(vals, cmap=cmap, vmin=0, vmax=max(10, vals.max()), aspect='auto')
    for (i, j), label in labels.items():
        ax.text(j, i, label, ha='center', va='center', color=pal['ink'], fontsize=12)
    ax.set_xticks(range(3), ['Alveolar / AT2', 'Fibroblast', 'Myeloid / Mph'], fontsize=9)
    ax.set_yticks(range(3), ['IPF 136831\nIPF / control', 'IPF 135893\nIPF / control', 'Mouse 300288\nanti-IL1B / IgG'], fontsize=10)
    ax.set_title('Independent-unit contrasts: case / reference counts', loc='left', fontsize=11)
    ax.tick_params(length=0)
    for spine in ax.spines.values(): spine.set_visible(False)
    ax = fig.add_subplot(gs[1, 1])
    contrasts = [(x, 'normal') for x in ['AAH', 'AIS', 'MIA', 'LUAD']] + [('LUAD', x) for x in ['AAH', 'AIS', 'MIA']]
    vals = np.zeros((3, 7))
    for i, comp in enumerate(comps):
        for j, (case, reference) in enumerate(contrasts):
            row = paired[paired.comp.eq(comp) & paired.case.eq(case) & paired.reference.eq(reference)]
            vals[i, j] = int(row.complete_patients.iloc[0]) if len(row) else 0
    ax.imshow(vals, cmap=cmap, vmin=0, vmax=max(10, vals.max()), aspect='auto')
    for i in range(3):
        for j in range(7): ax.text(j, i, f'{int(vals[i,j])}', ha='center', va='center', color=pal['ink'], fontsize=12)
    ax.set_xticks(range(7), [f'{a}\n− {b}' for a, b in contrasts], fontsize=9)
    ax.set_yticks(range(3), ['AT2-like', 'Fibroblasts', 'Macrophages'], fontsize=10)
    ax.set_title('Human lesions: complete paired patients', loc='left', fontsize=11)
    ax.tick_params(length=0)
    for spine in ax.spines.values(): spine.set_visible(False)
    fig.suptitle('Context and coverage determine which questions can be evaluated', x=.02, ha='left', fontsize=15)
    fig.text(.02, .02, 'Primary niche floor: ≥50 cells per compartment and ≥3 biological units per arm or complete patients. Counts vary by question.\nMouse myeloid is CSF1R+ and includes non-macrophages; alveolar cells include mixed AT2/AT1 states. Human lesion labels are reference-compatible.\nThe mouse epithelial-fraction endpoint instead requires ≥100 alveolar cells, leaving 3 controls / 2 treated mice, and additionally lacks a KAC classifier.\nPatient counts overlap across histologies and assays. No displayed counts may be summed into independent replication.', fontsize=9, color=pal['ink_2'])
    fig.subplots_adjust(left=.12, right=.98, bottom=.19, top=.92, wspace=.35, hspace=.30)
    for ext in ['png', 'svg']:
        fig.savefig(PAPER / 'figures' / f'context_design_and_coverage.{ext}', dpi=220, facecolor=pal['surface'])
    plt.close(fig)
    print(json.dumps(dict(independent_arm_rows=len(independent), paired_rows=len(paired), context_rows=len(design))))


if __name__ == '__main__':
    main()
