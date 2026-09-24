"""Render research-question figures from saved diagnostics and completed results."""
from pathlib import Path
import json
import sys
from datetime import datetime, timezone

PAPER = Path(__file__).resolve().parents[1]
ROOT = PAPER.parents[1]
sys.path.insert(0, str(ROOT))
from analysis.lib.provenance import sha256_file, code_identity, write_json_atomic


def main():
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, Normalize
    from matplotlib.lines import Line2D
    from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
    from matplotlib.ticker import PercentFormatter

    out = PAPER / 'trials/u7_proposal_figures'
    figdir = PAPER / 'figures'
    assert json.loads((out / 'preparation_run_record.json').read_text())['status'] == 'completed'
    pal = json.loads((ROOT / 'analysis/config/palette.json').read_text())
    cat = pal['categorical']
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'axes.titlesize': 12,
                         'axes.labelsize': 10, 'text.color': pal['ink'], 'axes.labelcolor': pal['ink'],
                         'xtick.color': pal['ink_2'], 'ytick.color': pal['ink_2'], 'svg.fonttype': 'none',
                         'figure.facecolor': pal['surface'], 'axes.facecolor': pal['surface'],
                         'savefig.facecolor': pal['surface']})
    seq = LinearSegmentedColormap.from_list('repository_magnitude', pal['sequential_ramp'])
    div = LinearSegmentedColormap.from_list('repository_difference', [cat['2'], pal['surface'], cat['1']])
    inputs = []
    outputs = []
    validations = []

    def read(rel):
        path = PAPER / rel
        inputs.append(dict(path=rel, sha256=sha256_file(path)))
        return pd.read_csv(path)

    def write_table(frame, name):
        assert not frame.empty, name
        frame.to_csv(out / name, index=False)

    def clean(ax, grid=None):
        for side in ['top', 'right']:
            ax.spines[side].set_visible(False)
        for side in ['left', 'bottom']:
            ax.spines[side].set_color(pal['axis'])
        ax.tick_params(length=3)
        if grid:
            ax.grid(axis=grid, color=pal['grid'], lw=.6)
            ax.set_axisbelow(True)

    def export(fig, stem):
        for ext in ['png', 'svg']:
            path = figdir / f'{stem}.{ext}'
            fig.savefig(path, dpi=200, bbox_inches='tight', pad_inches=.15)
            outputs.append(dict(path=path.relative_to(PAPER).as_posix(), sha256=sha256_file(path)))
        plt.close(fig)
        print('Rendered', stem, flush=True)

    def note(fig, text):
        fig.text(.035, .025, text, fontsize=9, color=pal['ink_2'], va='bottom', linespacing=1.5)

    def paired_violin(ax, pair, ylabel, title, percent=False):
        assert pair.index.is_unique and len(pair) >= 10
        values = [pair['normal'].to_numpy(), pair['LUAD'].to_numpy()]
        parts = ax.violinplot(values, positions=[0, 1], widths=.7, showextrema=False, showmedians=False)
        for body, color in zip(parts['bodies'], [pal['muted'], cat['1']]):
            body.set_facecolor(color)
            body.set_edgecolor(color)
            body.set_alpha(.16)
        jitter = np.linspace(-.10, .10, len(pair))
        for j, (_, row) in enumerate(pair.iterrows()):
            ax.plot([jitter[j], 1+jitter[j]], [row['normal'], row['LUAD']], color=pal['axis'], lw=.7, zorder=2)
        for x, label, color in [(0, 'normal', pal['muted']), (1, 'LUAD', cat['1'])]:
            ax.scatter(x+jitter, pair[label], color=color, s=24, edgecolor=pal['surface'], lw=.4, zorder=3)
            ax.plot([x-.13, x+.13], [pair[label].median()]*2, color=pal['ink'], lw=2, zorder=4)
        ax.set_xticks([0, 1], ['Normal', 'LUAD'])
        ax.set_ylabel(ylabel)
        ax.set_title(f'{title}\n{len(pair)} paired patients', loc='left')
        if percent:
            ax.set_ylim(-.02, 1.03)
            ax.set_yticks(np.linspace(0, 1, 6))
            ax.yaxis.set_major_formatter(PercentFormatter(1))
        clean(ax, 'y')

    # D0: geometry and source-allocation uncertainty; no new cell labels.
    cells = read('trials/u7_proposal_figures/umap_display_values.csv')
    assert cells.cell_id.is_unique and np.isfinite(cells[['UMAP1', 'UMAP2']]).all().all()
    sources = read('trials/u5_human_sources/IL1B_source_fractions.csv')
    src = sources[(sources.uncertainty == .2) & (sources.broad == 'Unassigned')]
    pair = src.pivot(index='patient', columns='histology', values='fraction_of_observed_IL1B_counts')[['normal', 'LUAD']].dropna()
    write_table(pair.reset_index(), 'D0_paired_unassigned_IL1B.csv')
    populations = ['AT2-like', 'Fibroblasts', 'Macrophages', 'Other assigned', 'Unassigned']
    colors = dict(zip(populations, [cat['1'], cat['2'], cat['3'], cat['4'], pal['muted']]))
    cells['display_population'] = cells.broad.where(cells.broad.isin(populations), 'Other assigned')
    fig, axs = plt.subplots(2, 2, figsize=(12.5, 10))
    fig.subplots_adjust(left=.07, right=.94, bottom=.16, top=.88, hspace=.45, wspace=.36)
    fig.suptitle('D0 | Where does source uncertainty enter the niche analysis?', x=.035, ha='left', fontsize=16, y=.97)
    fig.text(.035, .922, f'GSE308103: {len(cells):,} display cells from 555,480 QC nuclei; 23 patients', color=pal['ink_2'])
    plot = cells.sample(frac=1, random_state=20260924)
    axs[0, 0].scatter(plot.UMAP1, plot.UMAP2, c=plot.display_population.map(colors), s=2, alpha=.65, rasterized=True, linewidths=0)
    axs[0, 0].set_title('A  Existing reference-compatible populations', loc='left')
    axs[0, 0].legend(handles=[Line2D([], [], marker='o', ls='', color=colors[x], label=x, markersize=5) for x in populations],
                     fontsize=8, frameon=False, ncol=3, loc='upper left', bbox_to_anchor=(0, -.10))
    sc = axs[0, 1].scatter(plot.UMAP1, plot.UMAP2, c=plot.ann_finest_level_uncertainty, cmap=seq, vmin=0, vmax=1,
                         s=2, rasterized=True, linewidths=0)
    axs[0, 1].set_title('B  Fine-label uncertainty', loc='left')
    fig.colorbar(sc, ax=axs[0, 1], fraction=.045, pad=.025, label='Uncertainty (primary cutoff: 0.2)')
    axs[1, 0].scatter(plot.UMAP1, plot.UMAP2, c=pal['deemph'], s=2, rasterized=True, linewidths=0)
    positive = plot[plot.IL1B_count > 0].sort_values('IL1B_log1p_normalized_10000')
    sc = axs[1, 0].scatter(positive.UMAP1, positive.UMAP2, c=positive.IL1B_log1p_normalized_10000,
                         cmap=seq, vmin=0, vmax=float(positive.IL1B_log1p_normalized_10000.max()),
                         s=3, rasterized=True, linewidths=0)
    axs[1, 0].set_title('C  Measured IL1B RNA', loc='left')
    fig.colorbar(sc, ax=axs[1, 0], fraction=.045, pad=.025, label='log1p(counts per 10,000)')
    for ax in [axs[0, 0], axs[0, 1], axs[1, 0]]:
        ax.set_xlabel('UMAP 1')
        ax.set_ylabel('UMAP 2')
        ax.set_xticks([])
        ax.set_yticks([])
        clean(ax)
    paired_violin(axs[1, 1], pair, 'Fraction of recovered IL1B counts', 'D  RNA in unassigned cells', percent=True)
    note(fig, 'A-C: one fixed UMAP of saved reference-conditioned latent vectors; sampling balances patient / histology / source label.\n'
              'Point density is not abundance. D uses all QC cells; lines pair patients and black bars show medians.\n'
              'Label uncertainty and RNA allocation do not identify a new state or measure cytokine secretion. No labels were changed.')
    export(fig, 'derived_D0_source_annotation_context')
    validations.append(dict(check='D0_display_and_source_units', passed=bool(len(pair) == 23 and cells.patient.nunique() == 23)))

    # D1: donor PCA and existing program differences.
    pca = read('trials/u7_proposal_figures/AT2_pseudobulk_PCA_values.csv')
    pc_info = json.loads((out / 'PCA_summary.json').read_text())
    programs = read('trials/u6_human_specificity/patient_histology_scores.csv')
    programs = programs[(programs.config == 'unc20_pooled') & (programs.label == '__broad__') & (programs.cell_floor == 50) & programs.eligible]
    pairs = read('trials/u6_human_specificity/paired_program_values.csv')
    pairs = pairs[(pairs.config == 'unc20_pooled') & (pairs.label == '__broad__') & (pairs.cell_floor == 50)]
    hpcs = 'HPCS_without_ADI_or_operational_markers'
    hp = programs[programs.module == hpcs].pivot(index='patient', columns='histology', values='mean_logCPM')[['normal', 'LUAD']].dropna()
    write_table(hp.reset_index(), 'D1_HPCS_paired_scores.csv')
    hist_colors = {'normal': pal['muted'], 'AAH': cat['3'], 'AIS': cat['4'], 'MIA': cat['2'], 'LUAD': cat['1']}
    hist_markers = {'normal': 'o', 'AAH': '^', 'AIS': 's', 'MIA': 'D', 'LUAD': 'P'}
    fig, axs = plt.subplots(2, 2, figsize=(13.3, 10.8))
    fig.subplots_adjust(left=.10, right=.93, bottom=.19, top=.89, hspace=.47, wspace=.48)
    fig.suptitle('D1 | Which epithelial programs accompany lesion context?', x=.035, ha='left', y=.97, fontsize=16)
    fig.text(.035, .927, 'Reference-compatible AT2-like populations; histology does not establish malignant-cell identity', color=pal['ink_2'])
    for _, group in pca[pca.histology.isin(['normal', 'LUAD'])].groupby('patient'):
        if len(group) == 2:
            axs[0, 0].plot(group.PC1, group.PC2, color=pal['grid'], lw=.7, zorder=0)
    for hist, color in hist_colors.items():
        g = pca[pca.histology == hist]
        axs[0, 0].scatter(g.PC1, g.PC2, s=42, color=color, marker=hist_markers[hist], label=f'{hist} (n={len(g)})', edgecolor=pal['surface'], lw=.5)
    axs[0, 0].set_title('A  Patient-histology pseudobulk PCA', loc='left')
    axs[0, 0].set_xlabel(f'PC1 ({pc_info["explained_variance_ratio"][0]:.1%} variance)')
    axs[0, 0].set_ylabel(f'PC2 ({pc_info["explained_variance_ratio"][1]:.1%} variance)')
    axs[0, 0].legend(frameon=False, fontsize=8, ncol=2)
    clean(axs[0, 0])
    paired_violin(axs[0, 1], hp, 'Reduced-HPCS mean TMM log2 CPM', 'B  Shared-plasticity program')
    comparisons = [('AAH', 'normal'), ('AIS', 'normal'), ('MIA', 'normal'), ('LUAD', 'normal'), ('LUAD', 'AAH'), ('LUAD', 'AIS'), ('LUAD', 'MIA')]
    hpairs = pairs[pairs.module == hpcs]
    labels = []
    for i, (case, reference) in enumerate(comparisons):
        g = hpairs[(hpairs.case == case) & (hpairs.reference == reference)]
        axs[1, 0].scatter(g.difference, i+np.linspace(-.12, .12, len(g)), color=cat['1'], s=19)
        axs[1, 0].scatter(g.difference.mean(), i, marker='D', s=38, color=pal['ink'], zorder=3)
        labels.append(f'{case} - {reference}  (n={len(g)})')
    axs[1, 0].set_yticks(range(7), labels, fontsize=9)
    axs[1, 0].invert_yaxis()
    axs[1, 0].axvline(0, color=pal['axis'], lw=1)
    axs[1, 0].set_xlabel('Within-patient difference (log2 CPM)')
    axs[1, 0].set_title('C  Reduced-HPCS contrasts', loc='left')
    clean(axs[1, 0], 'x')
    write_table(hpairs, 'D1_HPCS_all_contrasts.csv')
    modules = {hpcs: 'HPCS reduced', 'Han_ISR_without_HPCS_ADI_or_operational_markers': 'ISR reduced',
               'AT2_published_holdout': 'AT2 program', 'AT1_published_holdout': 'AT1 program',
               'HALLMARK_HYPOXIA': 'Hypoxia', 'HALLMARK_INFLAMMATORY_RESPONSE': 'Inflammatory'}
    heat = pairs[(pairs.case == 'LUAD') & (pairs.reference == 'normal') & pairs.module.isin(modules)]
    mat = heat.pivot(index='patient', columns='module', values='difference').loc[:, list(modules)].sort_index()
    assert mat.shape == (23, 6) and mat.notna().all().all()
    limit = float(np.abs(mat.to_numpy()).max())
    im = axs[1, 1].imshow(mat, cmap=div, vmin=-limit, vmax=limit, aspect='auto', interpolation='none')
    axs[1, 1].set_xticks(range(6), list(modules.values()), rotation=35, ha='right', fontsize=8)
    axs[1, 1].set_yticks(range(23), [f'P{x}' for x in mat.index], fontsize=7)
    axs[1, 1].set_title('D  LUAD - normal, fixed programs', loc='left')
    fig.colorbar(im, ax=axs[1, 1], fraction=.045, pad=.025, label='Paired difference (log2 CPM)')
    write_table(heat, 'D1_program_heatmap_values.csv')
    note(fig, 'A: centered log2(CPM+1), 2,000 variable expressed genes; no histology-directed feature selection or batch correction.\n'
              'B-D: existing primary TMM scores, uncertainty <=0.2 and >=50 cells; bars in B are medians and diamonds in C are means.\n'
              'The HPCS score also increases in repair/IPF comparisons. Different panels reuse patients; these are not progression trajectories.')
    export(fig, 'derived_D1_epithelial_programs')
    reference = hpairs[(hpairs.case == 'LUAD') & (hpairs.reference == 'normal')].set_index('patient').difference.sort_index()
    validations.append(dict(check='D1_paired_scores_match_released_differences', passed=bool(np.allclose((hp.LUAD-hp.normal).sort_index(), reference))))

    # D2a: components and measured edge differences; all independent patient values retained.
    context = read('trials/u5_human_sources/broad_IL1_context.csv')
    genes = ['IL1B', 'IL1A', 'IL1R1', 'IL1RAP', 'IL1RN', 'IL1R2', 'SIGIRR']
    context = context[(context.uncertainty == .2) & context.histology.isin(['normal', 'LUAD']) &
                      context.broad.isin(['Macrophages', 'Fibroblasts', 'AT2-like']) & context.gene.isin(genes) & (context.cells >= 50)]
    dots = context.groupby(['broad', 'histology', 'gene']).agg(n_patients=('patient', 'nunique'),
            donor_mean_expression=('mean_normalized_10000', 'mean'), donor_mean_detection=('detection_fraction', 'mean')).reset_index()
    write_table(dots, 'D2_component_dot_values.csv')
    write_table(context, 'D2_component_donor_values.csv')
    norm = read('cache/u5_human_niche/normalized_components.csv.gz')
    norm = norm[(norm.config == 'unc20_pooled') & (norm.label == '__broad__') & (norm.cell_floor == 50) &
                (norm.prior_count == 1) & norm.histology.isin(['normal', 'LUAD'])]
    choices = [('macrophages', 'IL1B'), ('fibroblasts', 'IL1R1'), ('fibroblasts', 'IL1RAP'),
               ('AT2', 'IL1R1'), ('AT2', 'IL1RAP')]
    component_pairs = []
    for comp, gene in choices:
        g = norm[(norm.comp == comp) & (norm.gene == gene)].pivot(index='patient', columns='histology', values='logCPM').dropna()
        for patient, row in g.iterrows():
            component_pairs.append(dict(comp=comp, gene=gene, patient=patient, difference=row.LUAD-row.normal))
    cp = pd.DataFrame(component_pairs)
    write_table(cp, 'D2_paired_components.csv')
    edges = read('trials/u5_human_niche/primary_compatibility_patient_values.csv')
    edges = edges[(edges['view'] == 'broad') & (edges.resource == 'consensus') & (edges.case == 'LUAD') & (edges.reference == 'normal')]
    selections = [('IL1B', 'IL1R1_IL1RAP', 'macrophages', 'fibroblasts', 'IL1B: macrophage to fibroblast'),
                  ('IL1B', 'IL1R1_IL1RAP', 'macrophages', 'AT2', 'IL1B: macrophage to AT2-like'),
                  ('TGFB1', 'TGFBR1_TGFBR2', 'macrophages', 'fibroblasts', 'TGFB1: macrophage to fibroblast'),
                  ('CCL2', 'CCR2', 'fibroblasts', 'macrophages', 'CCL2: fibroblast to macrophage')]
    fig = plt.figure(figsize=(13, 9.8))
    grid = fig.add_gridspec(2, 2, left=.09, right=.94, bottom=.16, top=.88, height_ratios=[1.0, 1], wspace=.6, hspace=.65)
    ax = fig.add_subplot(grid[0, :])
    fig.suptitle('D2a | Separate ligand, recipient and regulatory RNA components', x=.035, ha='left', y=.97, fontsize=16)
    fig.text(.035, .924, 'GSE308103 | observed RNA profiles and LUAD-versus-normal paired contrasts', color=pal['ink_2'])
    row_order = [(b, h) for b in ['Macrophages', 'Fibroblasts', 'AT2-like'] for h in ['normal', 'LUAD']]
    vmax = float(np.log1p(dots.donor_mean_expression).max())
    for i, (broad, hist) in enumerate(row_order):
        g = dots[(dots.broad == broad) & (dots.histology == hist)].set_index('gene').loc[genes]
        sc = ax.scatter(range(len(genes)), [i]*len(genes), s=20+180*g.donor_mean_detection,
                        c=np.log1p(g.donor_mean_expression), cmap=seq, vmin=0, vmax=vmax, edgecolor=pal['ink_2'], lw=.4)
    ax.set_xticks(range(len(genes)), genes)
    row_labels = []
    for broad, hist in row_order:
        ns = dots[(dots.broad == broad) & (dots.histology == hist)].n_patients
        row_labels.append(f'{broad} | {hist} (n={ns.min()}-{ns.max()})' if ns.min() != ns.max() else f'{broad} | {hist} (n={ns.iloc[0]})')
    ax.set_yticks(range(len(row_order)), row_labels, fontsize=9)
    ax.invert_yaxis()
    ax.set_ylim(len(row_order)-.5, -.5)
    ax.set_xlim(-.5, len(genes)-.5)
    ax.axvline(1.5, color=pal['axis'], lw=.8)
    ax.axvline(3.5, color=pal['axis'], lw=.8)
    ax.set_title('A  Ligands | activating receptor components | regulatory context', loc='left', pad=12)
    fig.colorbar(sc, ax=ax, fraction=.025, pad=.025, label='log1p(donor-mean counts/10,000)')
    ax.legend(handles=[ax.scatter([], [], s=20+180*f, facecolor='none', edgecolor=pal['ink_2'], label=f'{f:.0%}') for f in [.1, .5, 1]],
              title='Mean detection', loc='upper center', bbox_to_anchor=(.5, -.13), ncol=3, frameon=False, fontsize=8, title_fontsize=8)
    clean(ax)
    ax = fig.add_subplot(grid[1, 0])
    labs = []
    for i, (comp, gene) in enumerate(choices):
        g = cp[(cp.comp == comp) & (cp.gene == gene)]
        ax.scatter(g.difference, i+np.linspace(-.10, .10, len(g)), s=21, color=cat['1'])
        ax.scatter(g.difference.mean(), i, s=36, marker='D', color=pal['ink'])
        labs.append(f'{comp}: {gene}\n(n={len(g)})')
    ax.set_yticks(range(len(labs)), labs, fontsize=8)
    ax.invert_yaxis()
    ax.axvline(0, color=pal['axis'], lw=1)
    ax.set_xlabel('Paired component difference (log2 CPM)')
    ax.set_title('B  Which component changes?', loc='left')
    clean(ax, 'x')
    ax = fig.add_subplot(grid[1, 1])
    plotted_edges = []
    labs = []
    for i, (lig, rec, source, target, label) in enumerate(selections):
        g = edges[(edges.ligand == lig) & (edges.receptor == rec) & (edges.source == source) & (edges.target == target)]
        assert len(g) >= 3 and g.patient.is_unique
        ax.scatter(g.paired_difference, i+np.linspace(-.10, .10, len(g)), s=21, color=cat['1'])
        ax.scatter(g.paired_difference.mean(), i, s=36, marker='D', color=pal['ink'])
        labs.append(label.replace(': ', ':\n')+f' (n={len(g)})')
        plotted_edges.append(g)
    ax.set_yticks(range(len(labs)), labs, fontsize=8)
    ax.invert_yaxis()
    ax.axvline(0, color=pal['axis'], lw=1)
    ax.set_xlabel('Paired RNA-compatibility difference')
    ax.set_title('C  Fixed directional candidates', loc='left')
    clean(ax, 'x')
    write_table(pd.concat(plotted_edges), 'D2_paired_edge_values.csv')
    note(fig, 'A: >=50 cells per donor / compartment; equal donor weighting within each histology. Components use all-assay normalization.\n'
              'B-C: each point is one complete patient; diamonds are means, with no cell-level P values.\n'
              'IL1RN, IL1R2 and SIGIRR are regulatory context. RNA compatibility does not measure cytokine release or causal communication.')
    export(fig, 'derived_D2_source_recipient_components')
    validations.append(dict(check='D2_exact_component_and_edge_pairs', passed=bool(not cp.duplicated(['comp', 'gene', 'patient']).any())))

    # D2b: all broad-receiver pathway rows for the chosen contrast, including null results.
    camera = read('trials/u5_human_niche/camera.csv')
    family = camera[(camera.config == 'unc20_pooled') & (camera.cell_floor == 50)]
    pw = family[(family.label == '__broad__') & (family.case == 'LUAD') & (family.reference == 'normal')].copy()
    pathways = sorted(pw['set'].unique())
    col_order = [(comp, setting) for comp in ['AT2', 'fibroblasts', 'macrophages'] for setting in ['estimated', 'fixed001']]
    names = {'AT2': 'AT2-like', 'fibroblasts': 'Fibroblast', 'macrophages': 'Macrophage'}
    fig, ax = plt.subplots(figsize=(12.7, 8.1))
    fig.subplots_adjust(left=.39, right=.94, bottom=.25, top=.83)
    fig.suptitle('D2b | Recipient enrichment depends on the correlation model', x=.035, ha='left', y=.97, fontsize=16)
    fig.text(.035, .911, 'LUAD - normal | all 21 declared, eligible broad-receiver tests in this contrast', color=pal['ink_2'])
    for i, pathway in enumerate(pathways):
        for j, (comp, setting) in enumerate(col_order):
            g = pw[(pw['set'] == pathway) & (pw.comp == comp) & (pw.correlation_setting == setting)]
            if g.empty:
                ax.scatter(j, i, marker='x', s=25, color=pal['axis'], linewidths=.6)
                continue
            assert len(g) == 1
            row = g.iloc[0]
            color = cat['1'] if row.Direction == 'Up' else cat['2']
            size = 22+30*min(4, -np.log10(max(row.FDR_global, 1e-300)))
            ax.scatter(j, i, s=size, facecolor=color if row.FDR_global < .05 else pal['surface'], edgecolor=color, lw=1.2,
                       marker='^' if row.Direction == 'Up' else 'v')
    ax.set_yticks(range(len(pathways)), [s.replace('HALLMARK_', '').replace('REACTOME_', '').replace('_', ' ').title() for s in pathways], fontsize=9)
    labels = []
    for comp, setting in col_order:
        g = pw[(pw.comp == comp) & (pw.correlation_setting == setting)]
        labels.append(f'{names[comp]}\n'+('Primary' if setting == 'estimated' else 'Sensitivity')+f'\nn={int(g.n_patients.iloc[0])}')
    ax.set_xticks(range(6), labels, fontsize=9)
    ax.set_xlim(-.6, 5.6)
    ax.set_ylim(len(pathways)-.5, -.5)
    ax.axvline(1.5, color=pal['axis'], lw=.8)
    ax.axvline(3.5, color=pal['axis'], lw=.8)
    ax.grid(color=pal['grid'], lw=.4)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0)
    counts = family.groupby('correlation_setting').FDR_global.agg(total='size', q05=lambda s: int((s < .05).sum()))
    note(fig, 'Up triangles: up; down triangles: down. Filled: global q < 0.05; open: tested without meeting that threshold.\n'
              'Size = -log10(q), capped at 4; x = pathway not assigned to this receiver in the declared test family.\n'
              f'Across all 279 primary-family tests: estimated correlation {counts.loc["estimated", "q05"]}/279; fixed 0.01 sensitivity {counts.loc["fixed001", "q05"]}/279.\n'
              'Both columns use the existing patient-blocked CAMERA analysis; competitive direction and q values are not GSEA NES or effect sizes.')
    write_table(pw, 'D2_pathway_values.csv')
    export(fig, 'derived_D2_recipient_enrichment')
    validations.append(dict(check='D2_pathway_pairs_and_full_family_counts', passed=bool(len(pw) == 42 and counts.loc['estimated', 'q05'] == 0 and counts.loc['fixed001', 'q05'] == 148)))

    # D2c: show every expression-supported ligand, not just IL1B or top ranks.
    fits = read('trials/u5_human_ligand_targets/primary_candidate_rankings.csv')
    stability = read('trials/u5_human_ligand_targets/candidate_ranking_stability.csv')
    take = (fits.label == '__broad__') & (fits.case == 'LUAD') & (fits.reference == 'normal') & (fits.direction == 'up') & (fits.scope == 'focused_triad')
    selected = fits[take].merge(stability, on=['comp', 'label', 'case', 'reference', 'direction', 'scope', 'ligand'], validate='one_to_one', suffixes=('', '_stability'))
    write_table(selected, 'D2_ligand_target_values.csv')
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 7.7))
    fig.subplots_adjust(left=.08, right=.97, bottom=.23, top=.78, wspace=.48)
    fig.suptitle('D2c | Is IL1B a leading candidate for recipient target changes?', x=.035, ha='left', y=.97, fontsize=16)
    fig.text(.035, .918, 'LUAD - normal upregulated targets | all eligible focused-triad candidates, per receiver', color=pal['ink_2'])
    for ax, comp in zip(axes, ['AT2', 'fibroblasts', 'macrophages']):
        g = selected[selected.comp == comp].sort_values('rank')
        for i, row in enumerate(g.itertuples()):
            color = cat['1'] if row.ligand == 'IL1B' else pal['muted']
            if np.isfinite(row.minimum_pearson) and np.isfinite(row.maximum_pearson):
                ax.plot([row.minimum_pearson, row.maximum_pearson], [i, i], color=pal['axis'], lw=2, zorder=1)
            ax.scatter(row.pearson, i, s=55 if row.ligand == 'IL1B' else 30, color=color,
                       marker='D' if row.ligand == 'IL1B' else 'o', zorder=2)
        ax.set_yticks(range(len(g)), [f'{int(r.rank)}. {r.ligand}' for r in g.itertuples()], fontsize=9)
        ax.set_ylim(len(g)-.5, -.5)
        ax.set_title(f'{names[comp]}\n{len(g)} eligible candidates', loc='left')
        ax.set_xlabel('Prior-target Pearson correlation')
        ax.axvline(0, color=pal['axis'], lw=.8)
        clean(ax, 'x')
    note(fig, 'Points: full-data fit; horizontal lines: eligible leave-one-patient-out ranges, not confidence intervals. Blue diamonds: IL1B.\n'
              'IL1B ranks 12/14 in AT2-like cells and 4/10 in fibroblasts; it is not expression-eligible in the macrophage panel.\n'
              'Different target/background sets prevent causal strength comparisons across receivers. Prior weights are unsigned; fit is not activation.\n'
              'IL1B omission coverage: AT2 23/23; fibroblasts 19/22. All candidate-specific omission counts are included in the plotted table.')
    export(fig, 'derived_D2_ligand_target_candidates')
    validations.append(dict(check='D2_all_ligand_candidates_retained', passed=bool(len(selected) == 32)))

    # D4: an experimental design, explicitly distinct from measured results.
    fig, ax = plt.subplots(figsize=(12.6, 6.8))
    fig.subplots_adjust(left=.02, right=.98, top=.82, bottom=.18)
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 5)
    ax.axis('off')
    fig.suptitle('D4 | Test recovery after IL-1beta withdrawal', x=.035, ha='left', y=.95, fontsize=17)
    fig.text(.035, .888, 'PROPOSED EXPERIMENT | design schematic only; no outcome data or expected-response curves', color=pal['ink_2'])
    boxes = [(.1, 2.7, 3.3, 1.8, 'Exposure schedule', 'Vehicle / transient / sustained\nTime-matched sampling\nMeasured cytokine exposure'),
             (4.3, 2.7, 3.3, 1.8, 'Recipient context', 'Epithelial-only / fibroblast co-culture\nRecipient-specific IL1R1 perturbation\nIndependent biological replicates'),
             (8.5, 2.7, 3.3, 1.8, 'Withdrawal and follow-up', 'Remove the external stimulus\nVerify target engagement\nFollow recovery over time'),
             (2.2, .1, 7.6, 1.7, 'Discriminating measurements', 'Epithelial mature-state recovery and viable yield\nFibroblast response, protein measurements and functional readouts\nLineage measurements if persistence or fate is the endpoint')]
    for x, y, w, h, title, body in boxes:
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.04', facecolor=pal['surface'], edgecolor=cat['1'], lw=1.4))
        ax.text(x+w/2, y+h-.35, title, ha='center', va='center', weight='bold', fontsize=11)
        ax.text(x+w/2, y+.65, body, ha='center', va='center', fontsize=9, linespacing=1.65)
    for start, end in [((3.45, 3.6), (4.22, 3.6)), ((7.65, 3.6), (8.42, 3.6)), ((10.15, 2.6), (8.4, 1.87)), ((5.95, 2.6), (5.95, 1.87))]:
        ax.add_patch(FancyArrowPatch(start, end, arrowstyle='-|>', mutation_scale=14, color=pal['ink_2'], lw=1.1))
    note(fig, 'RQ4: does recovery after withdrawal depend on exposure history and fibroblast IL-1 reception?\n'
              'Resolve dose, sampling times, genotype controls and replication in a dedicated experimental protocol.\n'
              'The completed cross-sectional RNA analysis cannot establish reversibility, persistent fate or malignant transformation.')
    export(fig, 'derived_D4_withdrawal_experiment_proposal')

    assert all(row['passed'] for row in validations), validations
    pd.DataFrame(inputs).drop_duplicates('path').to_csv(out / 'render_inputs.csv', index=False)
    write_json_atomic(out / 'render_validation.json', dict(status='passed', checks=validations, figure_count=6,
                      measured_figure_count=5, proposed_experiment_schematics=1, no_new_inferential_tests=True,
                      visual_review_status='pending', generated_utc=datetime.now(timezone.utc).isoformat()))
    write_json_atomic(out / 'render_run_record.json', dict(status='completed_pending_visual_review', code=code_identity(ROOT, __file__),
                      completed_utc=datetime.now(timezone.utc).isoformat(), outputs=outputs,
                      matplotlib_version=matplotlib.__version__, proposal_vs_results='Five data figures and one explicitly labeled experimental design'))


if __name__ == '__main__':
    main()
