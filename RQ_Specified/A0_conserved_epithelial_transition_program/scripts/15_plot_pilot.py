"""Render cohort coverage, the frozen discovery gates and qualifying gene effects."""
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[1]
OUT = BASE / 'tables/pilot_v1'
FIG = BASE / 'figures/pilot_v1'
BLUE = '#28699C'
GRAY = '#BBC3CA'


def main():
    if FIG.exists():
        raise SystemExit('Refusing to overwrite pilot figures')
    FIG.mkdir(parents=True)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'svg.hashsalt': 'A0-pilot-v1'})
    cfg = json.loads((BASE / 'config/pilot_v1.json').read_text())
    frozen = json.loads((OUT / 'frozen_programme.json').read_text())
    coverage = pd.read_csv(OUT / 'eligibility.tsv', sep='\t')
    genes = pd.read_csv(OUT / 'discovery_gene_effects.tsv.gz', sep='\t')
    panels = []

    def save(fig, name):
        for suffix in ['png', 'svg']:
            kwargs = {'metadata': {'Date': None}} if suffix == 'svg' else {}
            fig.savefig(FIG / f'{name}.{suffix}', dpi=180, bbox_inches='tight', **kwargs)
        plt.close(fig)
        panels.append(name)

    fig, axes = plt.subplots(1, 3, figsize=(12.2, 5.2), layout='constrained')
    for ax, role, title in zip(axes, ['D1', 'D2', 'V1'],
                                ['Repair: Strunz mice', 'Development: Sountoulidis donors', 'Intestine: Haber mice']):
        d = coverage[coverage.role == role].copy()
        d['minimum'] = d[['start', 'intermediate', 'destination']].min(axis=1)
        d = d.sort_values(['minimum', 'unit'])
        ax.barh(d.unit, d.minimum, color=[BLUE if p else GRAY for p in d.eligible])
        ax.axvline(cfg['cell_floor'], color='#B14F32', linestyle='--', linewidth=1)
        ax.set_title(f'{role}  {title}\n{int(d.eligible.sum())} eligible biological units', fontsize=10)
        ax.set_xlabel('Cells in the least represented state')
        ax.tick_params(axis='y', labelsize=8)
        ax.grid(axis='x', alpha=.15)
    fig.suptitle('Eligibility: all three states must have at least 30 cells', fontsize=13)
    save(fig, '01_biological_unit_coverage')

    allowed = ~(genes.label_excluded | genes.operational_excluded)
    detection = np.ones(len(genes), dtype=bool)
    for role in ['D1', 'D2']:
        detection &= genes[f'{role}_intermediate_detection_fraction'] + 1e-12 >= cfg['discovery']['minimum_detected_unit_fraction']
    masks = [np.ones(len(genes), dtype=bool), allowed.to_numpy(), allowed.to_numpy() & detection]
    labels = ['Common one-to-one orthologs', 'After fixed marker / symbol exclusions', 'Intermediate detection in both contexts']
    for role, endpoint, label in [('D1', 'start', 'Repair: intermediate > AT2'),
                                  ('D1', 'destination', 'Repair: intermediate > AT1'),
                                  ('D2', 'start', 'Development: intermediate > distal'),
                                  ('D2', 'destination', 'Development: intermediate > proximal')]:
        masks.append(masks[-1] & genes[f'{role}_vs_{endpoint}_pass'].to_numpy())
        labels.append(label)
    counts = [int(m.sum()) for m in masks]
    pd.DataFrame({'display_order': range(1, len(counts)+1), 'cumulative_gate': labels,
                  'genes_remaining': counts}).to_csv(OUT / 'discovery_gate_counts.tsv', sep='\t', index=False)
    fig, ax = plt.subplots(figsize=(10, 4.7), layout='constrained')
    yy = np.arange(len(labels))
    ax.barh(yy, counts, color=[GRAY]*(len(labels)-1)+[BLUE])
    ax.set_yticks(yy, labels)
    ax.invert_yaxis()
    ax.set_xscale('symlog', linthresh=1)
    ax.set_xlim(0, max(counts)*2.8)
    for y, count in zip(yy, counts):
        ax.text(max(count, .05), y, f'  {count:,}', va='center')
    ax.axvline(cfg['discovery']['minimum_genes'], color='#B14F32', ls='--', lw=1)
    ax.set_xlabel('Genes remaining (compressed scale; dashed line = 20-gene minimum)')
    ax.set_title('Discovery requires every gate; no endpoint can substitute for another')
    save(fig, '02_discovery_gates')

    d = genes[genes.eligible].sort_values(['minimum_median_effect', 'human'], ascending=[False, True])
    if len(d) > 50:
        d = d[d.selected]
    fig, axes = plt.subplots(1, 2, figsize=(9, max(3.5, .27*len(d)+1.8)), layout='constrained', sharey=True)
    if len(d):
        columns = [f'{role}_vs_{end}' for role in ['D1', 'D2'] for end in ['start', 'destination']]
        values = d[[c+'_median' for c in columns]].to_numpy()
        positive = d[[c+'_positive_fraction' for c in columns]].to_numpy()
        for ax, matrix, title, vmax in zip(axes, [values, positive],
                                          ['Median log2(CPM + 1) difference', 'Fraction of biological units > 0'],
                                          [max(.2, values.max()), 1]):
            im = ax.imshow(matrix, cmap='Blues', vmin=0, vmax=vmax, aspect='auto')
            ax.set_xticks(range(4), ['Repair\nvs AT2', 'Repair\nvs AT1', 'Development\nvs distal', 'Development\nvs proximal'], fontsize=8)
            ax.set_yticks(range(len(d)), d.human)
            ax.set_title(title, fontsize=10)
            for i in range(len(d)):
                for j in range(4):
                    ax.text(j, i, f'{matrix[i,j]:.2f}', ha='center', va='center', fontsize=8,
                            color='white' if matrix[i,j] > .65*vmax else '#17232D')
            fig.colorbar(im, ax=ax, fraction=.04, pad=.025)
    else:
        for ax in axes:
            ax.set_axis_off()
        axes[0].text(.1, .5, 'No genes satisfy every fixed discovery gate.', transform=axes[0].transAxes)
    fig.suptitle(f"{frozen['eligible_genes']} qualifying genes; {frozen['selected_genes']} selected\n"
                 'Selected discovery effects are descriptive, not independent validation', fontsize=12)
    save(fig, '03_qualifying_gene_effects')
    transfer_path = OUT / 'transfer_unit_differences.tsv'
    if transfer_path.exists():
        transfer = pd.read_csv(transfer_path, sep='\t')
        summary = pd.read_csv(OUT / 'transfer_summary.tsv', sep='\t')
        fig, axes = plt.subplots(1, 3, figsize=(11.8, 4.7), layout='constrained')
        for ax, role, title, endpoints in zip(axes, ['D1','D2','V1'],
                ['Repair: selected discovery', 'Development: selected discovery', 'Intestine: frozen transfer'],
                [['vs AT2','vs AT1'],['vs distal','vs proximal'],['vs stem','vs mature']]):
            d = transfer[transfer.role == role]
            units = sorted(d.unit.unique())
            color = BLUE if role != 'V1' else '#B14F32'
            for offset, unit in zip(np.linspace(-.13,.13,len(units)), units):
                y = [100*float(d.loc[(d.unit==unit)&(d.endpoint==e),'difference'].iloc[0]) for e in ['start','destination']]
                ax.plot(np.array([0,1])+offset,y,'o-',color=color,alpha=.6,lw=.65,ms=4)
            for i, endpoint in enumerate(['start','destination']):
                median = float(summary.loc[(summary.role==role)&(summary.endpoint==endpoint),'median_difference'].iloc[0])
                ax.plot([i-.22,i+.22],[100*median,100*median],color='#15232F',lw=2.5)
            ax.axhline(0,color='#697680',lw=.8,ls='--')
            ax.set_xticks([0,1],endpoints)
            ax.set_xlim(-.4,1.4)
            ax.set_title(f'{title}\nn = {len(units)} biological units',fontsize=10)
            ax.set_ylabel('Intermediate minus endpoint\n(mean percentile-rank score points)')
            ax.grid(axis='y',alpha=.15)
        fig.suptitle('The same frozen programme must exceed both endpoints in the new tissue',fontsize=12)
        save(fig, '04_frozen_transfer')
    record = {'generator_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'input_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                              for p in [OUT/'eligibility.tsv', OUT/'discovery_gene_effects.tsv.gz', OUT/'frozen_programme.json']},
              'figures': {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(FIG.iterdir())},
              'visual_review': 'pending'}
    if transfer_path.exists():
        for path in [transfer_path,OUT/'transfer_summary.tsv']:
            record['input_sha256'][path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
    (OUT / 'figure_render.json').write_text(json.dumps(record, indent=2)+'\n')
    print(json.dumps({'figures': panels, 'gate_counts': counts}))


if __name__ == '__main__':
    main()
