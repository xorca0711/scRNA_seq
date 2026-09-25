"""Report paired human target eligibility, expression-supported fits and omissions."""
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
    from matplotlib.colors import LinearSegmentedColormap

    out = PAPER / 'trials/u5_human_ligand_targets'
    record = json.loads((out / 'filter_run_record.json').read_text())
    assert record['status'] in ['completed', 'completed_no_eligible_target_sets']
    e = pd.read_csv(out / 'target_eligibility.csv')
    full_e = e[e.omitted_patient.eq('__none__')]
    parity = pd.read_csv(out / 'full_fit_DE_parity.csv')
    assert parity.maximum_q_difference.max() < 1e-8
    warning_audit = pd.read_csv(out / 'primary_prior_warning_audit.csv')
    assert len(warning_audit) == int(full_e.eligible.sum())
    assert warning_audit.nonfinite_exactly_constant.all()
    assert warning_audit.maximum_finite_score_difference.max() < 1e-12
    ranks = pd.read_csv(out / 'primary_candidate_rankings.csv') if record['primary_candidate_rows'] else pd.DataFrame()
    stable = pd.read_csv(out / 'candidate_ranking_stability.csv') if len(ranks) else pd.DataFrame()
    keys = ['comp', 'label', 'case', 'reference', 'direction', 'scope']
    if len(ranks):
        assert not ranks.duplicated(keys + ['ligand']).any()
        for _, g in ranks.groupby(keys):
            assert np.array_equal(g['rank'], g.pearson.rank(method='min', ascending=False))
            assert g.ranked_ligands.eq(len(g)).all()
        assert (stable.eligible_omissions <= stable.planned_omissions).all()
    comp_order = ['AT2', 'fibroblasts', 'macrophages']
    rows = []
    for comp in comp_order:
        for case, reference in CONTRASTS:
            for direction in ['up', 'down']:
                gate = full_e[full_e.comp.eq(comp) & full_e.label.eq('__broad__') & full_e.case.eq(case) & full_e.reference.eq(reference) & full_e.direction.eq(direction)]
                z = ranks[ranks.comp.eq(comp) & ranks.label.eq('__broad__') & ranks.case.eq(case) & ranks.reference.eq(reference) & ranks.direction.eq(direction) & ranks.scope.eq('focused_triad') & ranks.ligand.eq('IL1B')] if len(ranks) else ranks
                status = 'fewer_than_three_paired_patients' if gate.empty else ('target_gate_not_met' if not gate.eligible.iloc[0] else ('source_receiver_expression_gate_not_met' if z.empty else 'eligible'))
                rows.append(dict(comp=comp, case=case, reference=reference, direction=direction, status=status,
                                 n_patients=int(gate.n_patients.iloc[0]) if len(gate) else np.nan,
                                 mapped_targets=int(gate.mapped_targets.iloc[0]) if len(gate) else np.nan,
                                 pearson=float(z.pearson.iloc[0]) if len(z) else np.nan,
                                 rank=float(z['rank'].iloc[0]) if len(z) else np.nan,
                                 ranked_ligands=int(z.ranked_ligands.iloc[0]) if len(z) else 0))
    plotted = pd.DataFrame(rows)
    plotted.to_csv(out / 'broad_IL1B_figure_values.csv', index=False)
    pal = json.loads((ROOT / 'analysis/config/palette.json').read_text())
    cmap = LinearSegmentedColormap.from_list('prior_fit', [pal['categorical']['2'], pal['surface'], pal['categorical']['1']])
    cmap.set_bad(pal['deemph'])
    maximum = max(.1, float(plotted.pearson.abs().max())) if plotted.pearson.notna().any() else .1
    fig, axes = plt.subplots(3, 1, figsize=(13, 10.5))
    fig.set_facecolor(pal['surface'])
    for panel, ax in enumerate(axes):
        data = np.full((3, 7), np.nan)
        texts = {}
        for i, comp in enumerate(comp_order):
            for j, (case, reference) in enumerate(CONTRASTS):
                x = plotted[plotted.comp.eq(comp) & plotted.case.eq(case) & plotted.reference.eq(reference)]
                if panel == 0:
                    up, down = x[x.direction.eq('up')].iloc[0], x[x.direction.eq('down')].iloc[0]
                    texts[(i, j)] = '<3 patients' if pd.isna(up.mapped_targets) else f'{int(up.mapped_targets)} up\n{int(down.mapped_targets)} down'
                    data[i, j] = 0 if pd.notna(up.mapped_targets) else np.nan
                else:
                    row = x[x.direction.eq('up' if panel == 1 else 'down')].iloc[0]
                    data[i, j] = row.pearson
                    texts[(i, j)] = f'{row.pearson:.3f}\nrank {int(row["rank"])}/{row.ranked_ligands}' if row.status == 'eligible' else {
                        'fewer_than_three_paired_patients': '<3 patients', 'target_gate_not_met': 'Targets <10',
                        'source_receiver_expression_gate_not_met': 'Expression\ngate'}[row.status]
        im = ax.imshow(data, cmap=cmap, vmin=-maximum, vmax=maximum, aspect='auto')
        for (i, j), label in texts.items():
            ax.text(j, i, label, ha='center', va='center', fontsize=9, color=pal['ink'])
        ax.set_xticks(range(7), [f'{a} − {b}' for a, b in CONTRASTS], fontsize=9)
        ax.set_yticks(range(3), ['AT2-like', 'Fibroblasts', 'Macrophages'])
        ax.tick_params(length=0)
        ax.set_title(['A  Mapped DE target counts (each direction needs ≥10)', 'B  IL1B prior fit to upregulated targets', 'C  IL1B prior fit to downregulated targets'][panel], loc='left', fontsize=11)
        for spine in ax.spines.values():
            spine.set_visible(False)
    fig.suptitle('Human paired ligand–target prioritization: eligibility precedes ranking', x=.02, ha='left', fontsize=14)
    fig.text(.02, .025, 'Broad receiver views; focused fixed-source triad candidates. Patient-blocked DE/BH, q<0.05 targets; same-patient expression support.\nScores are unsigned-prior Pearson fits, not activation or inhibition. Downregulated-target fit does not establish ligand repression.\nRanks use all expression-eligible planned ligands. Full patient-omission refits and conditional support are reported in the stability table.\nGrey cells are unevaluable, not zero effects; healthy-reference-compatible labels do not establish nonmalignant identity.', fontsize=9, color=pal['ink_2'])
    fig.subplots_adjust(left=.12, right=.98, top=.92, bottom=.21, hspace=.45)
    cbax = fig.add_axes([.68, .125, .28, .015])
    colorbar = fig.colorbar(im, cax=cbax, orientation='horizontal')
    colorbar.set_label('Unsigned-prior Pearson fit', fontsize=8, labelpad=2)
    colorbar.ax.tick_params(labelsize=8)
    for ext in ['png', 'svg']:
        fig.savefig(PAPER / 'figures' / f'human_ligand_target_eligibility_and_fit.{ext}', dpi=220, facecolor=pal['surface'])
    plt.close(fig)
    relevant = stable[stable.label.eq('__broad__') & stable.scope.eq('focused_triad') & stable.ligand.isin(['IL1A', 'IL1B'])] if len(stable) else stable
    report = ['# Paired human ligand–target prioritization', '',
              f'Completed {len(parity)} full-data receiver/contrast DE parity checks. Of {len(full_e)} fitted primary receiver/contrast/direction target sets, {int(full_e.eligible.sum())} pass the ten-mapped-target gate. Expression filtering leaves {len(ranks)} candidate rows across broad and subtype views and the two planned source scopes; these overlapping views are not independent confirmations.', '',
              'The statistic is Pearson correlation between the checksum-verified NicheNet v2 unsigned prior and a binary receiver target vector. Each receiver is fit with a patient-blocked design; up- and downregulated q<0.05 target sets remain separate. Every eligible omitted-patient analysis refits filtering, TMM, voom, DE/BH and target selection. Source/receiver expression must retain at least three same-case-histology patients with the same source label, ≥50 cells per group and ≥10% detection of all required subunits. Omission eligibility is checked again.', '',
              'Canonical IL-1 candidates require IL1R1 and IL1RAP; IL1R2 and SIGIRR do not qualify as activating receptors. All assayed genes in the exact pseudobulks are available for expression gates. Candidate ranks are restricted to planned ligand families; raw prior scores without expression support are not findings. Neither a high rank nor complete omission coverage constitutes causal or independent validation.', '',
              f'The original R log contains an aggregate warning notice. A follow-up reconstruction of all {len(warning_audit)} eligible primary target sets reproduced only zero-standard-deviation warnings from constant prior columns. Every undefined primary correlation exactly matches a constant column, and every finite score agrees within 1e-12. Undefined scores are excluded from candidate ranks, not replaced by zero. The audit does not claim that the original run saved each individual warning message.', '',
              '## Broad-receiver IL-1 stability', '',
              '| Receiver | Contrast | Targets | Ligand | Rank | Eligible/planned omissions | Rank range |',
              '|---|---|---|---|---:|---:|---|']
    for row in relevant.itertuples():
        rank_range = f'{row.min_rank:g}–{row.max_rank:g}' if pd.notna(row.min_rank) else 'not estimable'
        report.append(f'| {row.comp} | {row.case} − {row.reference} | {row.direction} | {row.ligand} | {row.primary_rank:g} | {row.eligible_omissions}/{row.planned_omissions} | {rank_range} |')
    if relevant.empty:
        report.append('| No eligible broad-receiver IL-1 ranks | | | | | | |')
    report += ['', '![Human ligand targets](../../figures/human_ligand_target_eligibility_and_fit.png)', '',
               '[Target eligibility](target_eligibility.csv), [plotted gate/fit values](broad_IL1B_figure_values.csv), [full-data DE parity](full_fit_DE_parity.csv), [constant-prior score audit](primary_prior_warning_audit.csv).']
    if len(ranks):
        report += ['', '[Candidate ranks](primary_candidate_rankings.csv), [omission stability](candidate_ranking_stability.csv), [fixed source support](fixed_source_receiver_support.csv), [prior edge assay coverage](prior_edge_assay_coverage.csv).']
    (out / 'REPORT.md').write_text('\n'.join(report) + '\n', encoding='utf-8')
    qa = dict(status='passed', full_fit_DE_parity_checks=len(parity), primary_target_sets=len(full_e),
              eligible_primary_target_sets=int(full_e.eligible.sum()), candidate_rows=len(ranks), rank_recomputation_passed=True,
              omission_denominators_checked=True, broad_IL1B_figure_cells=len(plotted),
              primary_prior_score_reconstructions=len(warning_audit), all_undefined_primary_scores_have_constant_prior=True)
    (out / 'validation.json').write_text(json.dumps(qa, indent=2) + '\n')
    print(json.dumps(qa))


if __name__ == '__main__':
    main()
