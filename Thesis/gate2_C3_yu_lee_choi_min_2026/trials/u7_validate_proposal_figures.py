"""Check source parity, plotted units, provenance and inspected figure hashes."""
from pathlib import Path
import json
import sys
from datetime import datetime, timezone

PAPER = Path(__file__).resolve().parents[1]
ROOT = PAPER.parents[1]
sys.path.insert(0, str(ROOT))
from analysis.lib.provenance import sha256_file, write_json_atomic, code_identity


def main():
    import numpy as np
    import pandas as pd
    out = PAPER / 'trials/u7_proposal_figures'
    prep = json.loads((out / 'preparation_run_record.json').read_text())
    render = json.loads((out / 'render_run_record.json').read_text())
    visual = json.loads((out / 'visual_review.json').read_text())
    assert prep['status'] == 'completed'
    assert json.loads((out / 'render_validation.json').read_text())['status'] == 'passed'
    for record in [prep, render]:
        path = ROOT / record['code']['entrypoint']
        assert sha256_file(path) == record['code']['entrypoint_sha256'], str(path)
    hashes = pd.concat([pd.read_csv(out / x) for x in ['preparation_inputs.csv', 'render_inputs.csv']]).drop_duplicates('path')
    for row in hashes.itertuples():
        assert sha256_file(PAPER / row.path) == row.sha256, row.path
    cells = pd.read_csv(out / 'umap_display_values.csv')
    coverage = pd.read_csv(out / 'embedding_sampling_coverage.csv')
    sources = pd.read_csv(PAPER / 'trials/u5_human_sources/IL1B_source_fractions.csv')
    source = sources[sources.uncertainty == .2]
    comparison = coverage.merge(source, on=['patient', 'histology', 'broad'], how='outer', validate='one_to_one')
    assert comparison.all_QC_cells.notna().all() and comparison.cells.notna().all()
    assert np.array_equal(comparison.all_QC_cells, comparison.cells)
    actual = cells.groupby(['patient', 'histology', 'broad']).size().rename('actual').reset_index()
    check = coverage.merge(actual, on=['patient', 'histology', 'broad'], validate='one_to_one')
    assert np.array_equal(check.display_cells, check.actual)
    assert len(cells) == 34178 and cells.cell_id.is_unique and cells.patient.nunique() == 23
    assert int(coverage.all_QC_cells.sum()) == 555480 and coverage.display_cells.between(1, 75).all()
    assert np.isfinite(cells[['UMAP1', 'UMAP2', 'IL1B_log1p_normalized_10000']]).all().all()
    pca = pd.read_csv(out / 'AT2_pseudobulk_PCA_values.csv')
    assert len(pca) == 70 and not pca.duplicated(['patient', 'histology']).any()
    assert pca.cells.ge(50).all() and pca.comp.eq('AT2').all() and pca.label.eq('__broad__').all()
    assert np.allclose(pca[[f'PC{i}' for i in range(1, 6)]].mean(), 0, atol=1e-8)
    programs = pd.read_csv(out / 'D1_program_heatmap_values.csv')
    assert not programs.duplicated(['patient', 'module']).any() and len(programs) == 23*6
    paired = pd.read_csv(out / 'D1_HPCS_paired_scores.csv').set_index('patient')
    differences = pd.read_csv(out / 'D1_HPCS_all_contrasts.csv')
    delta = differences[(differences.case == 'LUAD') & (differences.reference == 'normal')].set_index('patient').difference
    assert np.allclose((paired.LUAD-paired.normal).sort_index(), delta.sort_index())
    pathways = pd.read_csv(out / 'D2_pathway_values.csv')
    assert len(pathways) == 42 and len(pathways[pathways.correlation_setting == 'estimated']) == 21
    assert not pathways.duplicated(['comp', 'set', 'correlation_setting']).any()
    targets = pd.read_csv(out / 'D2_ligand_target_values.csv')
    assert targets.groupby('comp').size().to_dict() == {'AT2': 14, 'fibroblasts': 10, 'macrophages': 8}
    assert not targets.duplicated(['comp', 'ligand']).any()
    for row in render['outputs']:
        assert sha256_file(PAPER / row['path']) == row['sha256'], row['path']
    pngs = [row for row in render['outputs'] if row['path'].endswith('.png')]
    assert len(pngs) == 6 and len(visual['figures']) == 6
    for row in pngs:
        reviewed = visual['figures'][Path(row['path']).name]
        assert reviewed['status'] == 'visually_reviewed' and reviewed['png_sha256'] == row['sha256']
    record = dict(status='passed', checked_utc=datetime.now(timezone.utc).isoformat(), code=code_identity(ROOT, __file__),
                  input_hashes_checked=len(hashes), source_label_cell_count_parity=True, full_QC_cells=555480,
                  display_cells=len(cells), patients=23, PCA_units=70, program_pairs_verified=True,
                  pathway_rows=42, target_candidates=32, visually_reviewed_figures=6,
                  measured_figures=5, proposal_schematics=1, no_new_inferential_tests=True)
    write_json_atomic(out / 'validation.json', record)
    print(json.dumps({k: v for k, v in record.items() if k != 'code'}))


if __name__ == '__main__':
    main()
