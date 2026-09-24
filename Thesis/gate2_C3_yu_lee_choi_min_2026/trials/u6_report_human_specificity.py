"""Report paired histology contrasts of frozen source programs in candidate AT2 cells."""
from pathlib import Path
import json

PAPER = Path(__file__).resolve().parents[1]
ROOT = PAPER.parents[1]
CONTRASTS = [(x, 'normal') for x in ['AAH', 'AIS', 'MIA', 'LUAD']] + [('LUAD', x) for x in ['AAH', 'AIS', 'MIA']]


def main():
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    out = PAPER / 'trials/u6_human_specificity'
    scores = pd.read_csv(out / 'patient_histology_scores.csv')
    keys = ['config', 'label', 'cell_floor', 'module']
    assert not scores.duplicated(keys + ['patient', 'histology']).any()
    pairs, summaries = [], []
    for key, group in scores.groupby(keys):
        for case, reference in CONTRASTS:
            joined = group[group.histology.eq(case)].merge(
                group[group.histology.eq(reference)], on='patient', suffixes=('_case', '_reference'), validate='one_to_one')
            covered = bool(group.eligible.all())
            values = joined.mean_logCPM_case - joined.mean_logCPM_reference
            n = len(joined)
            base = dict(zip(keys, key), case=case, reference=reference)
            loo = (values.sum() - values) / (n - 1) if n > 1 else np.array([np.nan])
            summaries.append(dict(base, n_patients=n, eligible=covered and n >= 3,
                                  source_coverage=float(group.assayed_source_fraction.min()),
                                  reason='assayed_source_gene_fraction_below_0.7' if not covered else ('fewer_than_3_complete_patients' if n < 3 else 'eligible_descriptive_contrast'),
                                  mean_difference=float(values.mean()) if covered and n else np.nan,
                                  positive_patients=int((values > 0).sum()) if covered and n else np.nan,
                                  loo_min=float(np.min(loo)) if covered and n > 1 else np.nan,
                                  loo_max=float(np.max(loo)) if covered and n > 1 else np.nan))
            if covered:
                for patient, difference in zip(joined.patient, values):
                    pairs.append(dict(base, patient=patient, difference=float(difference)))
    summary, values = pd.DataFrame(summaries), pd.DataFrame(pairs)
    summary.to_csv(out / 'paired_program_summary.csv', index=False)
    values.to_csv(out / 'paired_program_values.csv', index=False)
    primary = summary[summary.config.eq('unc20_pooled') & summary.cell_floor.eq(50) & summary.label.eq('__broad__')]
    dots = values[values.config.eq('unc20_pooled') & values.cell_floor.eq(50) & values.label.eq('__broad__')]
    for row in primary[primary.mean_difference.notna()].itertuples():
        x = dots[dots.module.eq(row.module) & dots.case.eq(row.case) & dots.reference.eq(row.reference)]
        assert x.patient.is_unique and len(x) == row.n_patients
        assert np.isclose(x.difference.mean(), row.mean_difference)
    sensitivity = summary.merge(primary[['case', 'reference', 'module', 'mean_difference', 'eligible']].rename(
        columns={'mean_difference': 'primary_difference', 'eligible': 'primary_eligible'}), on=['case', 'reference', 'module'], validate='many_to_one')
    sensitivity['same_direction_when_both_eligible'] = np.where(
        sensitivity.eligible & sensitivity.primary_eligible, np.sign(sensitivity.mean_difference).eq(np.sign(sensitivity.primary_difference)), np.nan)
    sensitivity.to_csv(out / 'program_sensitivity_comparison.csv', index=False)

    modules = ['HPCS_author_top100', 'HPCS_without_ADI_or_operational_markers', 'ADI_published_holdout',
               'Han_ISR_source_symbols', 'Han_ISR_without_HPCS_ADI_or_operational_markers']
    names = ['HPCS source', 'HPCS minus ADI/labels', 'ADI holdout', 'ISR source', 'ISR minus HPCS/ADI/labels']
    pal = json.loads((ROOT / 'analysis/config/palette.json').read_text())
    fig, axes = plt.subplots(2, 4, figsize=(16, 8), sharey=True)
    fig.set_facecolor(pal['surface'])
    for ax, (case, reference) in zip(axes.flat, CONTRASTS):
        for i, module in enumerate(modules):
            x = dots[dots.case.eq(case) & dots.reference.eq(reference) & dots.module.eq(module)]
            ax.scatter(x.difference, i + np.linspace(-.15, .15, len(x)), s=18, color=pal['categorical']['1'])
            gate = primary[primary.case.eq(case) & primary.reference.eq(reference) & primary.module.eq(module)]
            if len(gate) and gate.source_coverage.iloc[0] < .7:
                ax.text(.5, i, f'Assay coverage {gate.source_coverage.iloc[0]:.0%} <70%',
                        transform=ax.get_yaxis_transform(), ha='center', va='center', fontsize=8, color=pal['ink_2'])
            if len(x) >= 3:
                ax.scatter(x.difference.mean(), i, marker='D', s=45, facecolor=pal['surface'], edgecolor=pal['ink'], zorder=3)
        n = dots[dots.case.eq(case) & dots.reference.eq(reference) & dots.module.eq(modules[0])].patient.nunique()
        ax.set_title(f'{case} − {reference} | n={n} patients', fontsize=10, loc='left')
        ax.set_yticks(range(len(modules)), names)
        ax.axvline(0, color=pal['axis'], lw=.8)
        ax.grid(axis='x', color=pal['grid'], lw=.5)
        ax.set_facecolor(pal['surface'])
        ax.tick_params(labelsize=8, length=0)
        ax.set_xlabel('Paired mean log2 CPM difference', fontsize=9)
        for spine in ax.spines.values():
            spine.set_visible(False)
        if n == 0:
            ax.text(.5, .5, 'No complete patient pairs', transform=ax.transAxes, ha='center')
    axes[0, 0].invert_yaxis()
    axes.flat[-1].axis('off')
    fig.suptitle('Human lesion context: source programs in reference-compatible AT2-like cells', x=.02, ha='left', fontsize=14)
    fig.text(.02, .025, 'Each dot is one patient with both histologies; diamonds mark means with at least three paired patients. No cell-level P values.\nPrimary: confidence uncertainty ≤0.2, ≥50 cells per patient/histology, pooled repeated libraries, all-gene TMM.\nHealthy-reference compatibility does not establish nonmalignant identity; HPCS/ISR scores do not assign KAC, HPCS or malignant states.', fontsize=9, color=pal['ink_2'])
    fig.tight_layout(rect=[0, .15, 1, .92])
    for ext in ['png', 'svg']:
        fig.savefig(PAPER / 'figures' / f'human_epithelial_program_specificity.{ext}', dpi=220, facecolor=pal['surface'])
    plt.close(fig)

    report = ['# Human epithelial program contrasts', '',
              'Frozen source-defined HPCS, ADI, DATP, alveolar, stress and Han ISR programs were scored in healthy-reference-compatible AT2-like pseudobulks. These candidate labels can include lesional cells; they are not an independent malignant-cell classifier.', '',
              'Each difference retains the same patient on both sides. Repeated libraries are pooled before scoring. All-assay TMM log2 CPM uses prior count 1; primary coverage is 50 cells and uncertainty ≤0.2. The 30/100-cell, uncertainty ≤0.3 and largest-library analyses are separate sensitivities. Source coverage uses the original mouse source-list denominator after strict one-to-one ortholog mapping, with a 0.7 gate. Missing source genes are not zeros.', '',
              '## Reduced HPCS program', '']
    for row in primary[primary.module.eq('HPCS_without_ADI_or_operational_markers')].itertuples():
        report.append(f'- {row.case} minus {row.reference}: {row.n_patients} paired patients, mean {row.mean_difference:+.3f} log2 CPM, {int(row.positive_patients) if pd.notna(row.positive_patients) else 0}/{row.n_patients} positive; {row.reason}.')
    report += ['', 'The ADI holdout has 64.1% source coverage and the six-gene DATP/PATS holdout has 66.7%, so both remain ineligible under the 70% rule. Their empty score rows are missing coverage, not zero biological effects. The complete source HPCS list retains 78.0%, and the overlap-reduced HPCS list retains 79.1%; the measured assay subset is disclosed rather than treated as the complete original signature.', '',
               'These descriptive program contrasts cannot replace the prespecified source-defined KAC/NF-κB association. The public KAC classifier remains unavailable. Similar program scores across lesion, repair or fibrosis contexts establish neither common ancestry nor transformation. Overlapping patients across comparisons are not independent replications.', '',
               '![Human epithelial programs](../../figures/human_epithelial_program_specificity.png)', '',
               '[All paired summaries](paired_program_summary.csv), [patient values](paired_program_values.csv), [sensitivities](program_sensitivity_comparison.csv), [source modules](module_specification.json).']
    (out / 'REPORT.md').write_text('\n'.join(report) + '\n', encoding='utf-8')
    qa = dict(status='passed', primary_broad_contrasts=len(primary), eligible_primary_broad_contrasts=int(primary.eligible.sum()),
              patient_unique=True, primary_effects_recomputed=True, no_cell_level_inference=True)
    (out / 'validation.json').write_text(json.dumps(qa, indent=2) + '\n')
    print(json.dumps(qa))


if __name__ == '__main__':
    main()
