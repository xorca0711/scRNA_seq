"""All-QC IL1B source attribution and processing/antagonist RNA context."""
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

    out = PAPER / 'trials/u5_human_sources'
    assert json.loads((out / 'run_record.json').read_text())['status'] == 'completed'
    sources = pd.read_csv(out / 'IL1B_source_fractions.csv')
    context = pd.read_csv(out / 'IL1_context_profiles.csv')
    parity = pd.read_csv(out / 'all_cell_count_parity.csv')
    assert len(parity) == 150 and parity.exact_count_parity.all()
    sums = sources.groupby(['patient', 'histology', 'uncertainty']).fraction_of_observed_IL1B_counts.sum(min_count=1)
    assert np.allclose(sums.dropna(), 1)
    context = context.groupby(['patient', 'histology', 'uncertainty', 'broad', 'gene']).agg(
        cells=('cells', 'sum'), counts=('count_sum', 'sum'), detected_cells=('detected_cells', 'sum'),
        sum_normalized_10000=('sum_normalized_10000', 'sum')).reset_index()
    context['detection_fraction'] = context.detected_cells / context.cells
    context['mean_normalized_10000'] = context.sum_normalized_10000 / context.cells
    context.to_csv(out / 'broad_IL1_context.csv', index=False)
    primary = sources[sources.uncertainty.eq(.2)]
    summary = primary.groupby(['histology', 'broad']).agg(
        n_patients=('patient', 'nunique'), median_count_fraction=('fraction_of_observed_IL1B_counts', 'median'),
        min_count_fraction=('fraction_of_observed_IL1B_counts', 'min'), max_count_fraction=('fraction_of_observed_IL1B_counts', 'max'),
        median_detection_fraction=('detection_fraction', 'median')).reset_index()
    summary.to_csv(out / 'IL1B_source_summary.csv', index=False)
    histologies = ['normal', 'AAH', 'AIS', 'MIA', 'LUAD']
    broad = ['Macrophages', 'Other myeloid', 'Fibroblasts', 'AT2-like', 'Other epithelial', 'Other assigned', 'Unassigned']
    pal = json.loads((ROOT / 'analysis/config/palette.json').read_text())
    fig, axes = plt.subplots(1, 5, figsize=(16, 6.5), sharey=True)
    fig.set_facecolor(pal['surface'])
    for ax, hist in zip(axes, histologies):
        for i, label in enumerate(broad):
            g = primary[primary.histology.eq(hist) & primary.broad.eq(label)]
            ax.scatter(g.fraction_of_observed_IL1B_counts, i + np.linspace(-.15, .15, len(g)), s=16, color=pal['categorical']['1'])
            if len(g):
                ax.scatter(g.fraction_of_observed_IL1B_counts.median(), i, marker='D', s=40, facecolor=pal['surface'], edgecolor=pal['ink'], zorder=3)
        n = primary[primary.histology.eq(hist)].patient.nunique()
        ax.set_title(f'{hist} | {n} patients', loc='left', fontsize=11)
        ax.set_xlim(-.02, 1.02)
        ax.set_xticks([0, .5, 1])
        ax.set_yticks(range(len(broad)), broad)
        ax.set_xlabel('Fraction of observed IL1B counts', fontsize=9)
        ax.grid(axis='x', color=pal['grid'], lw=.5)
        ax.set_facecolor(pal['surface'])
        ax.tick_params(length=0, labelsize=9)
        for spine in ax.spines.values():
            spine.set_visible(False)
    axes[0].invert_yaxis()
    fig.suptitle('Human IL1B RNA sources, including cells without confident reference labels', x=.02, ha='left', fontsize=14)
    fig.text(.02, .035, 'Dots = patient/histology aggregates; diamonds = within-histology medians. All QC cells enter exactly one group.\nCount fractions depend on recovered cell mixtures and library sizes; they are not tissue abundance or secretion fractions.\nUnassigned cells remain visible. Processing, receptor, decoy and antagonist RNA values are available separately; none measures mature IL-1β.', fontsize=9, color=pal['ink_2'])
    fig.tight_layout(rect=[0, .17, 1, .91])
    for ext in ['png', 'svg']:
        fig.savefig(PAPER / 'figures' / f'human_IL1B_sources.{ext}', dpi=220, facecolor=pal['surface'])
    plt.close(fig)
    unassigned = primary[primary.broad.eq('Unassigned')].groupby('histology').agg(
        patients=('patient', 'nunique'), median_IL1B_count_fraction=('fraction_of_observed_IL1B_counts', 'median')).reset_index()
    report = ['# Human IL-1 source and context profiles', '',
              'All QC cells from 75 libraries are represented, including cells whose finest reference-label uncertainty exceeds the primary 0.2 cutoff. Repeated libraries are pooled within patient and histology. Raw counts agree exactly with the retained source panel in all 150 library/cutoff checks. Both 0.2 and 0.3 confidence views are preserved.', '',
              'The figure shows the fraction of observed IL1B counts assigned to each broad population. This is a recovered-cell and library-dependent RNA allocation, not a tissue-composition correction or a secretion measurement. Patient-level detection fractions and full-library-normalized mean expression accompany raw count fractions. Missing or unassayed panel genes are recorded rather than filled with biological zeros.', '',
              'A substantial fraction of IL1B RNA belongs to cells without confident finest-level labels. The medians below are patient-level count fractions, not pooled cohort fractions. This limits attribution to specific cell types and makes the confidently mapped niche comparisons conditional on label retention; it does not invalidate the measured RNA in the unassigned population.', '',
              '| Histology | Patients | Median IL1B count fraction in unassigned cells |', '|---|---:|---:|']
    for row in unassigned.itertuples():
        report.append(f'| {row.histology} | {row.patients} | {row.median_IL1B_count_fraction:.1%} |')
    report += ['', 'IL1RN, IL1R2 and SIGIRR describe antagonist/decoy context. NLRP3, PYCARD, CASP1 and GSDMD are RNA measurements of processing-related components; they cannot establish inflammasome assembly, cleavage, release or mature IL-1β protein. Receptor and antagonist expression should be read alongside ligand RNA before proposing an activation model.', '',
              '![IL1B source allocation](../../figures/human_IL1B_sources.png)', '',
              '[Patient source fractions](IL1B_source_fractions.csv), [broad processing/receptor/antagonist context](broad_IL1_context.csv), [source summary](IL1B_source_summary.csv), [assay coverage](source_panel_coverage.csv), [count checks](all_cell_count_parity.csv).']
    (out / 'REPORT.md').write_text('\n'.join(report) + '\n', encoding='utf-8')
    qa = dict(status='passed', library_cutoff_parity_checks=len(parity), source_fraction_sums_checked=len(sums), all_QC_cells_accounted_for=True)
    (out / 'validation.json').write_text(json.dumps(qa, indent=2) + '\n')
    print(json.dumps(qa))


if __name__ == '__main__':
    main()
