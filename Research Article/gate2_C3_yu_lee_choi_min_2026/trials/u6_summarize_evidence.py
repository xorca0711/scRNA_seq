"""Assemble a reviewable evidence register after every feasible package finishes."""
from pathlib import Path
import json

PAPER = Path(__file__).resolve().parents[1]
ROOT = PAPER.parents[1]


def main():
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    def read(relative):
        return json.loads((PAPER / relative).read_text())

    human = read('trials/u5_human_niche/validation.json')
    targets = read('trials/u5_human_ligand_targets/validation.json')
    source = read('trials/u5_human_sources/validation.json')
    programs = read('trials/u6_human_specificity/validation.json')
    for record in [human, targets, source, programs]:
        assert record['status'] == 'passed'
    coverage = read('trials/u6_completion/pathway_coverage_validation.json')
    assert coverage['status'] == 'passed_identical_eligibility'
    out = PAPER / 'trials/u6_completion'
    matrix = pd.DataFrame([
        ['Direct mouse blockade', 'RNA niches measured', '0 / 46 primary q<0.05', 'KAC endpoint unidentifiable', 'No KAC classifier; 3 controls / 2 treated at phenotype floor', 'trials/u4_mouse_niche/REPORT.md'],
        ['IPF: two cohorts', 'RNA niches / targets measured', 'No primary q<0.05', 'IL1B compatibility differs by cohort', 'Observational; complete triads sparse', 'trials/u5_ligand_targets/REPORT.md'],
        ['Human lesion cohort', 'Paired RNA niches / targets measured', f'{human["primary_pathway_q05"]} / {human["primary_pathway_tests"]} primary q<0.05', 'Candidate atlas compartments', 'Much IL1B remains unassigned; no KAC classifier', 'trials/u5_human_niche/REPORT.md'],
        ['Human spatial companion', '56 section maps / paired means', 'Descriptive whole-section values', 'Pathological regions unidentifiable', 'No independent ROI labels; shared RNA patients', 'trials/u5_spatial_context/REPORT.md'],
        ['Post-viral repair context', '9 whole-sample RNA profiles', 'Descriptive; small biological n', 'Neighborhoods unidentifiable', 'No deposited coordinates; anti-CD8 in mice', 'trials/u5_spatial_context/REPORT.md'],
        ['Epithelial specificity', 'Source and reduced programs measured', 'Reduced HPCS also rises in repair/IPF', 'Increased score not neoplasia-specific', 'Pooled wells, sparse pairs, coupled source labels', 'trials/u6_specificity/REPORT.md'],
        ['Han ISR source', '129-symbol signature evaluated', '24 / 24 technical-seed directions agree', 'Original expression arm conditional', 'No reusable original processed matrices found', 'trials/u6_isr_extension/REPORT.md'],
    ], columns=['Context', 'Completed measurement', 'Evidence at declared setting', 'Interpretation', 'Material limit', 'Report'])
    matrix.to_csv(out / 'evidence_matrix.csv', index=False)
    pal = json.loads((ROOT / 'analysis/config/palette.json').read_text())
    fig, ax = plt.subplots(figsize=(18, 8.2), facecolor=pal['surface'])
    ax.axis('off')
    display = matrix.drop(columns='Report').copy()
    # Explicit wraps keep complete qualifications legible in the standalone export.
    import textwrap
    limits = [25, 32, 29, 31, 40]
    values = [[textwrap.fill(str(value), width=limits[j]) for j, value in enumerate(row)] for row in display.values]
    tab = ax.table(cellText=values, colLabels=display.columns, cellLoc='left', colLoc='left',
                   colWidths=[.15, .20, .20, .20, .25], bbox=[0, 0, 1, 1])
    tab.auto_set_font_size(False); tab.set_fontsize(10)
    for (i, j), cell in tab.get_celld().items():
        cell.set_facecolor(pal['deemph'] if i == 0 else pal['surface'])
        cell.set_edgecolor(pal['grid']); cell.set_text_props(color=pal['ink'])
    fig.suptitle('Evidence review: measured, inconclusive and unidentifiable are different outcomes', x=.02, ha='left', fontsize=15)
    fig.text(.02, .025, 'Primary pathway q-values use the declared estimated-correlation CAMERA families. Fixed-correlation sensitivities do not replace them.\nRNA compatibility, unsigned ligand-target ranks and shared program scores do not establish secretion, communication, lineage or malignant conversion.\nThis register is for interpretation review; it does not promote new causal claims into the repository claim contract.', fontsize=10, color=pal['ink_2'])
    fig.subplots_adjust(left=.02, right=.98, top=.91, bottom=.15)
    for ext in ['png', 'svg']:
        fig.savefig(PAPER / 'figures' / f'cross_context_evidence_review.{ext}', dpi=220, facecolor=pal['surface'])
    plt.close(fig)
    p = pd.read_csv(PAPER / 'trials/u5_human_niche/camera.csv')
    primary = p[p.config.eq('unc20_pooled') & p.cell_floor.eq(50) & p.correlation_setting.eq('estimated')]
    significant = primary[primary.FDR_global.lt(.05)].copy()
    significant.to_csv(out / 'human_primary_pathways_q05.csv', index=False)
    hpcs = pd.read_csv(PAPER / 'trials/u6_human_specificity/paired_program_summary.csv')
    hpcs = hpcs[hpcs.config.eq('unc20_pooled') & hpcs.label.eq('__broad__') & hpcs.cell_floor.eq(50) & hpcs.module.eq('HPCS_without_ADI_or_operational_markers')]
    text = ['# Completed analysis: evidence for joint review', '',
            'The feasible public-data analyses and their declared sensitivities are complete. No threshold was lowered after viewing results. Unsupported identities, absent spatial annotations and unavailable processed data remain explicit limits, rather than being replaced by proxy claims. The initial batch remains preserved in INITIAL_RUN_REPORT.md; this document summarizes the completed continuation.', '',
            '## Findings that remain available for interpretation', '',
            '1. **RNA niches can be measured without establishing an IL-1-specific causal circuit.** Both IPF cohorts support donor-level source, receptor, pathway and compatibility analyses. The estimated canonical macrophage-to-fibroblast IL1B compatibility contrast is near zero in GSE136831 (+0.005) and positive in GSE135893 (+1.285 log2-CPM-based units); this is not an equivalence test. No primary IPF pathway passes global q<0.05. The two cohorts do not jointly establish one uniform IL-1 response.', '',
            f'2. **Human paired contrasts are directly evaluated.** All 75 GSE308103 libraries from 23 patient labels were processed. The primary pathway family has {human["primary_pathway_q05"]}/{human["primary_pathway_tests"]} q<0.05 tests; {human["primary_RNA_contrasts"]} primary RNA-compatibility contrasts retain paired patient values. Conditional ligand-target fits yield {targets["candidate_rows"]} expression-supported candidate rows across overlapping views/scopes, not independent validation. Unassigned cells carry median IL1B count fractions of 52–72% across histologies; confident-subtype results cannot establish the dominant source across all recovered cells. Detailed findings and patient coverage are in the human reports.', '',
            '3. **An increased HPCS score is not specific to neoplasia in these comparisons.** The overlap-reduced source signature is higher in all seven evaluable pooled repair/developmental libraries and all three eligible KRT5-/KRT17+ versus AT2 IPF donor pairs. This is descriptive counterevidence to score specificity; it does not establish that repair or fibrotic cells are HPCS, KAC or malignant.', '',
            '4. **Spatial measurements are available at a narrower scale than the proposed niche hypothesis.** All 56 human sections and nine post-viral matrices were processed. Deposited human coordinates support measured maps and patient-level whole-section summaries. Missing independent pathology regions prevent regional enrichment claims; absent post-viral coordinates prevent neighborhood tests. Whole-section means cannot reproduce or refute a region-specific source result.', '',
            '## Human reduced-HPCS contrasts', '', '| Contrast | Paired patients | Mean log2 CPM difference | Positive patients | Gate |', '|---|---:|---:|---:|---|']
    for row in hpcs.itertuples():
        text.append(f'| {row.case} − {row.reference} | {row.n_patients} | {row.mean_difference:+.3f} | {int(row.positive_patients) if pd.notna(row.positive_patients) else 0} | {row.reason} |')
    text += ['', 'These atlas-compatible AT2-like contrasts remain separate from the unavailable source-defined KAC endpoint. Different histology contrasts reuse patients and are not independent confirmations.', '',
             '## Questions the public data still cannot identify', '',
             '- The primary KAC response to anti-IL-1β, and the KAC/NF-κB joint association: the public author classifier is missing. At the 100-alveolar-cell phenotype floor only three controls and two treated animals remain; lowering the floor does not recover KAC identity.',
             '- Later mouse blockade effects tied to treatment strategy: the required treatment-history crosswalk remains unresolved.',
             '- Mature IL-1β release, neutralization efficacy, mechanical stiffness, functional fibroblast suppression, irreversible fate, lineage ancestry and chromatin memory: RNA measurements do not directly assay these outcomes.',
             '- Independent confirmation from shared specimens, pooled wells or sorting gates: these units cannot be relabelled as additional animals or cohorts. GSE222901/GSE300293 received contextual source/design audits, not a new independent expression-validation run. The GSE277777 source-label pilot is a consistency/specificity analysis with unresolved biological independence.',
             '- Original Han developmental expression contrasts: the public source audit recovered the ISR signature and code, but not the reusable processed matrices needed for that conditional arm.', '',
             '## What a change in criteria would and would not do', '',
             'The predeclared cell floors, confidence, prior-count, sampling, LR resource/detection and correlation sensitivities are retained beside the primary results. All 52 cohort/pathway assay-coverage checks reach at least 92%, so the 50%, 70% and 80% pathway-coverage cutoffs select the same sets. A lower cutoff cannot fix absent labels, sample identity, regions or replication. Tested results without q<0.05 are inconclusive at the declared family threshold; they are not evidence of equivalence or absence.', '',
             '## Next interpretation decisions', '',
             'Review the evidence by biological claim rather than by the number of nominally positive panels: reproducible niche association; IL-1 specificity versus generic inflammatory/stress programs; and specificity of dysplastic/neoplastic states versus repair. A focused next experiment would vary IL-1β duration and withdrawal while measuring recipient-specific response, viable mature fate, protein release and lineage outcomes in independent biological replicates. Missing author KAC annotations and independent spatial region labels are targeted data requests; no outreach has been sent.', '',
             '![Context evidence](figures/cross_context_evidence_review.png)', '',
             '[Full completion register](WORK_PACKAGES.md) · [Figure gallery](README.md#figure-gallery) · [Evaluation rules](EVALUATION_RULES.md) · [Evidence matrix](trials/u6_completion/evidence_matrix.csv) · [Human primary pathway discoveries](trials/u6_completion/human_primary_pathways_q05.csv).']
    correlation = pd.read_csv(out / 'pathway_correlation_sensitivity.csv')
    detail = ['The pathway findings depend materially on the correlation assumption. The matched comparisons below use identical test membership and gene counts; the primary analysis remains unchanged.', '',
              '| Cohort | Tests | Estimated-correlation q<0.05 | Fixed-0.01 q<0.05 | Median estimated correlation |',
              '|---|---:|---:|---:|---:|']
    for row in correlation.itertuples():
        detail.append(f'| {row.cohort} | {row.tests} | {row.primary_estimated_q05} | {row.fixed001_q05} | {row.median_estimated_correlation:.3f} |')
    detail += ['', 'These are alternative modeling assumptions, not different biological replications. The choice of estimated correlation is a project decision rather than a universal field standard. The discrepancy should be discussed before making a pathway-level claim; it cannot be resolved merely by choosing the setting with more discoveries.', '',
               '[Correlation sensitivity and checks](trials/u6_completion/pathway_correlation_validation.json).', '']
    at = text.index('## Next interpretation decisions')
    text[at:at] = detail
    (PAPER / 'EVIDENCE_REVIEW.md').write_text('\n'.join(text) + '\n', encoding='utf-8')
    print(json.dumps(dict(contexts=len(matrix), human_primary_discoveries=len(significant), status='evidence_review_ready')))


if __name__ == '__main__':
    main()
