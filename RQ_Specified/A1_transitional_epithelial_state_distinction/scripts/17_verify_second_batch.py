"""Verify saved A1 evidence without fitting models or replacing run records.

Run with the scientific environment. Writes a fresh report plus complete
denominator/window diagnostics; refuses existing output directories.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import gzip
import hashlib
import importlib.util
import json
import sys
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[1]
S = BASE / 'tables/ire1_stability_2026-09-25'
D = BASE / 'tables/direct_marks_2026-09-25'


def sha(path):
    with path.open('rb') as h:
        return hashlib.file_digest(h, 'sha256').hexdigest()


def read(path):
    return pd.read_csv(path, sep='\t')


def close(actual, expected):
    np.testing.assert_allclose(actual, expected, rtol=1e-10, atol=1e-12, equal_nan=True)


def bh(p):
    p = np.asarray(p)
    order = np.argsort(p)
    adjusted = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p)+1))[::-1])[::-1]
    result = np.empty(len(p))
    result[order] = np.minimum(1, adjusted)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-id', default='second_batch_verification')
    args = parser.parse_args()
    if not args.run_id.replace('_', '').replace('-', '').isalnum():
        parser.error('Invalid run ID')
    out = BASE / 'tables' / args.run_id
    rp = BASE / 'reports' / f'{args.run_id}.json'
    if out.exists() or rp.exists():
        raise SystemExit('Refusing to overwrite verification evidence')
    hashes = {}
    for name in ['ire1_run', 'ire1_stability_run', 'direct_mark_run', 'second_batch_summary_run']:
        path = BASE / 'reports' / f'{name}.json'
        record = json.loads(path.read_text())
        assert record['status'] == 'completed'
        hashes[path.relative_to(BASE).as_posix()] = sha(path)
        for family in ['input_sha256', 'output_sha256']:
            for relative, digest in record.get(family, {}).items():
                if relative not in hashes:
                    hashes[relative] = sha(BASE / relative)
                assert hashes[relative] == digest, relative
    archive = BASE / 'reports/execution_sources/second_batch_2026-09-25'
    for name, run in [('11_run_ire1_stability.py', 'ire1_stability_run'),
                      ('12_fit_ire1_stability.R', 'ire1_stability_run'),
                      ('15_quantify_direct_marks.py', 'direct_mark_run'),
                      ('16_summarize_second_batch.py', 'second_batch_summary_run')]:
        record = json.loads((BASE/'reports'/f'{run}.json').read_text())
        expected = record['code_sha256']
        if isinstance(expected, dict):
            expected = expected[name]
        assert sha(archive/name) == expected
    contract = read(BASE/'tables/ire1/input_contract.tsv').set_index('key').value
    for key, name in [('counts', 'processed/ire1/counts.tsv'), ('samples', 'tables/ire1/sample_manifest.tsv')]:
        with (BASE/name).open('rb') as h:
            assert hashlib.file_digest(h, 'md5').hexdigest() == contract[key+'_md5']
    meta = read(BASE/'tables/ire1/sample_manifest.tsv')
    primary = read(BASE/'tables/ire1/gene_effects.tsv').set_index('feature_id')
    full = read(S/'all_gene_sensitivity.tsv')
    focus = read(S/'focus_effects.tsv')
    diag = read(S/'fit_diagnostics.tsv')
    assert meta.mouse.nunique() == 10 and len(primary) == 14811
    assert set(diag.run) == {f'omit_{m}' for m in meta.mouse} | {'S061'}
    expected_focus = set(read(BASE/'tables/ire1/marker_mapping.tsv').feature_id) | set(primary.index[primary.FDR < .05])
    pd.testing.assert_frame_equal(focus.reset_index(drop=True),
        full[full.feature_id.isin(expected_focus)].reset_index(drop=True))
    for row in diag.itertuples():
        m = meta[meta.batch == 'S061'] if row.run == 'S061' else meta[meta.mouse != int(row.omitted_mouse)]
        norm = read(S/f'{row.run}_normalization.tsv')
        assert list(m.sample_id) == list(norm.sample_id)
        design = [np.ones(len(m)), (m.sex == 'M').astype(int), (m.group == 'KIRA8').astype(int)]
        if row.run != 'S061':
            design.append((m.batch == 'S135').astype(int))
        rank = np.linalg.matrix_rank(np.column_stack(design))
        assert rank == row.rank == len(design)
        assert len(m)-rank == row.residual_df and len(m) == row.n
        assert (m.group == 'Vehicle').sum() == row.Vehicle >= 3
        assert (m.group == 'KIRA8').sum() == row.KIRA8 >= 3
        f = full[full.run == row.run].set_index('feature_id')
        assert list(f.index) == list(primary.index)
        close(f.FDR, bh(f.PValue))
        assert (f.FDR < .05).sum() == row.genes_FDR05
        close(f.logFC.corr(primary.logFC), row.Pearson_to_primary)
        close(f.logFC.corr(primary.logFC, method='spearman'), row.Spearman_to_primary)
    sets = read(BASE/'tables/ire1/gene_set_mapping_tested.tsv')
    coverage = read(BASE/'tables/ire1/gene_set_coverage.tsv')
    pathways = read(S/'pathway_sensitivity.tsv')
    eligible = coverage[coverage.eligible].set_index('gene_set')
    for name, f in pathways.groupby('run'):
        assert set(f.gene_set) == set(eligible.index) == {'TGFbeta', 'Gene Ontology 0006986'}
        close(f.FDR, bh(f.PValue))
        for row in f.itertuples():
            members = set(sets[(sets.gene_set == row.gene_set) & sets.retained].feature_id)
            assert members <= set(primary.index) and len(members) == row.NGenes == eligible.loc[row.gene_set, 'unique_tested_genes']
    summary = read(BASE/'tables/second_batch_summary/ire1_focus_stability.tsv')
    markers = read(BASE/'tables/ire1/marker_mapping.tsv')
    small = read(S/'S135_descriptive.tsv').set_index('feature_id')
    assert len(summary) == 12 and (summary.selection == 'frozen_marker').sum() == 8
    assert set(summary[summary.selection == 'post_selection_primary_FDR_hit'].feature_id) == set(primary.index[primary.FDR < .05])
    close(small.log2_ratio_offset_0_5 if 'log2_ratio_offset_0_5' in small else small['log2_ratio_offset_0.5'],
          np.log2((small.KIRA8_mean_CPM+.5)/(small.Vehicle_mean_CPM+.5)))
    assert not any(c in small for c in ['PValue', 'FDR'])
    for row in summary.itertuples():
        p = primary.loc[row.feature_id]
        loo = full[(full.feature_id == row.feature_id) & (full.kind == 'leave_one_mouse_out')]
        within = full[(full.feature_id == row.feature_id) & (full.run == 'S061')].iloc[0]
        close([row.primary_logFC, row.primary_FDR, row.LOMO_min_logFC, row.LOMO_max_logFC,
               row.LOMO_min_FDR, row.LOMO_max_FDR, row.LOMO_max_abs_change, row.S061_logFC,
               row.S061_FDR, row.S135_descriptive_log2_ratio],
              [p.logFC, p.FDR, loo.logFC.min(), loo.logFC.max(), loo.FDR.min(), loo.FDR.max(),
               abs(loo.logFC-p.logFC).max(), within.logFC, within.FDR, small.loc[row.feature_id, 'log2_ratio_offset_0.5']])
        assert row.LOMO_n == len(loo) == 10
        assert row.LOMO_same_direction == (np.sign(loo.logFC) == np.sign(p.logFC)).sum()
        assert row.LOMO_FDR05 == (loo.FDR < .05).sum()
        if row.selection == 'frozen_marker':
            assert row.label == markers.set_index('feature_id').loc[row.feature_id, 'symbol']
    print('J1: provenance, 11 fits, gene universe, BH, pathways and summaries verified', flush=True)
    sys.path.insert(0, str(BASE/'cache/pybigtools'))
    import pybigtools
    spec = importlib.util.spec_from_file_location('direct', BASE/'scripts/15_quantify_direct_marks.py')
    direct = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(direct)
    cfg = json.loads((BASE/'config/direct_mark_loci.json').read_text())
    inventory = json.loads((BASE/'reports/direct_mark_input_inventory.json').read_text())
    loci = read(D/'loci.tsv')
    for ref in inventory['files']:
        if ref.get('accession') != 'reference' or ref.get('assembly') == 'mm10':
            continue
        actual = pd.DataFrame(direct.annotations(BASE/ref['path'], ref['assembly'], cfg['gene_groups']))
        saved = loci[loci.assembly == ref['assembly']].reset_index(drop=True)
        pd.testing.assert_frame_equal(actual[saved.columns], saved, check_dtype=False)
    sig = read(D/'histone_signals.tsv')
    manifest = read(D/'sample_manifest.tsv')
    assert len(manifest) == 24 and not manifest.duplicated(['state', 'replicate', 'mark']).any()
    assert len(sig) == 1104 and sig.gene.nunique() == 23
    assert set(manifest.assembly) == set(sig.assembly) == {'T2T-CHM13v2.0'}
    expected_samples = {(s, r, m) for s in ['iAT1','iAT2','iATCs'] for r in ['CUT1','CUT2'] for m in ['H3','H3K27ac','H3K4me3','H3K27me3']}
    assert set(zip(manifest.state, manifest.replicate, manifest.mark)) == expected_samples
    inv_by_gsm = {r.get('gsm'): r for r in inventory['files'] if r.get('gsm')}
    chroms = None
    checks = []
    for sample in manifest.itertuples():
        source = inv_by_gsm[sample.gsm]
        assert source['title'] == f'{sample.state}_{sample.replicate}_{sample.mark}'
        assert source['path'] == sample.source_path and source['sha256'] == sample.sha256
        bw = pybigtools.open(str(BASE/sample.source_path))
        if chroms is None:
            chroms = bw.chroms()
        assert bw.chroms() == chroms
        # Base-resolution API independently checks saved integral/coverage for
        # every track at the three loci frozen for the figure (72 windows).
        for row in sig[(sig.gsm == sample.gsm) & (sig.window == 'promoter') & sig.gene.isin(cfg['locus_panels'])].itertuples():
            values = bw.values(row.chrom, int(row.start), int(row.end), missing=float('nan'))
            arr = np.asarray(values)
            mean = np.nan_to_num(arr, nan=0).mean()
            covered = np.isfinite(arr).mean()
            close([mean, covered], [row.mean_CPM, row.covered_fraction])
            checks.append(dict(gsm=sample.gsm, gene=row.gene, mean_CPM=mean, covered_fraction=covered))
        bw.close()
    # Native CHM13 identifiers/sizes, additionally checked against a compact
    # assembly index when available; no liftover or cross-assembly joins.
    fai = BASE/'cache/followup_sources/chm13v2.0.fa.gz.fai'
    assembly_index_verified = False
    if fai.exists():
        ref_sizes = {f[0]: int(f[1]) for f in (line.split() for line in fai.read_text().splitlines())}
        assert all(ref_sizes[c] == size for c, size in chroms.items())
        hashes[fai.relative_to(BASE).as_posix()] = sha(fai)
        assembly_index_verified = True
    assert chroms['chr1'] == 248387328 and chroms['chr2'] == 242696752
    li = loci[loci.assembly == 'T2T-CHM13v2.0'].set_index('gene')
    h3 = sig[sig.mark == 'H3'].set_index(['gene','state','replicate','window'])
    for row in sig.itertuples():
        locus = li.loc[row.gene]
        width = cfg['promoter_half_width_bp' if row.window == 'promoter' else 'broad_promoter_half_width_bp']
        assert row.chrom == locus.chrom and row.tss == locus.tss and row.strand == locus.strand
        assert row.start == max(0, row.tss-width) and row.end == min(chroms[row.chrom], row.tss+width)
        c = h3.loc[(row.gene,row.state,row.replicate,row.window)]
        close([row.H3_mean_CPM,row.H3_covered_fraction], [c.mean_CPM,c.covered_fraction])
        ok = row.covered_fraction >= .8 and c.covered_fraction >= .8 and row.mean_CPM > 0 and c.mean_CPM > 0
        assert row.ratio_eligible == ok
        close(row.log2_mark_over_H3, np.log2(row.mean_CPM/c.mean_CPM) if ok else np.nan)
    contrasts = read(D/'histone_contrasts.tsv')
    idx = sig.set_index(['gene','state','replicate','window','mark'])
    for row in contrasts.itertuples():
        a = idx.loc[(row.gene,'iATCs',row.replicate,row.window,row.mark)]
        b = idx.loc[(row.gene,row.contrast.removeprefix('iATCs_minus_'),row.replicate,row.window,row.mark)]
        assert row.eligible == bool(a.ratio_eligible and b.ratio_eligible)
        close(row.log2_ratio_difference, a.log2_mark_over_H3-b.log2_mark_over_H3)
        close(row.raw_mark_log2_ratio, np.log2(a.mean_CPM/b.mean_CPM) if a.mean_CPM > 0 and b.mean_CPM > 0 else np.nan)
        close(row.H3_log2_ratio, np.log2(a.H3_mean_CPM/b.H3_mean_CPM) if a.H3_mean_CPM > 0 and b.H3_mean_CPM > 0 else np.nan)
    keys = ['gene','gene_group','mark','replicate','contrast']
    windows = contrasts[contrasts.window == 'promoter'].merge(
        contrasts[contrasts.window == 'broad_promoter'], on=keys, suffixes=('_1kb','_5kb'), validate='one_to_one')
    windows['both_windows_eligible'] = windows.eligible_1kb & windows.eligible_5kb
    windows['same_direction_across_windows'] = np.where(windows.both_windows_eligible,
        np.sign(windows.log2_ratio_difference_1kb) == np.sign(windows.log2_ratio_difference_5kb), np.nan)
    for width in ['1kb','5kb']:
        windows[f'H3_changes_direction_{width}'] = np.where(windows[f'eligible_{width}'],
            np.sign(windows[f'raw_mark_log2_ratio_{width}']) != np.sign(windows[f'log2_ratio_difference_{width}']), np.nan)
    concordance = read(BASE/'tables/second_batch_summary/histone_replicate_concordance.tsv')
    ci = contrasts.set_index(['gene','mark','window','contrast','replicate'])
    for row in concordance.itertuples():
        a = ci.loc[(row.gene,row.mark,row.window,row.contrast,'CUT1')]
        b = ci.loc[(row.gene,row.mark,row.window,row.contrast,'CUT2')]
        assert row.both_eligible == bool(a.eligible and b.eligible)
        close([row.CUT1_log2_ratio_difference,row.CUT2_log2_ratio_difference], [a.log2_ratio_difference,b.log2_ratio_difference])
        if row.both_eligible:
            assert row.direction_concordant == (np.sign(a.log2_ratio_difference) == np.sign(b.log2_ratio_difference))
    print('J2: 24 headers, 72 base-resolution windows, annotations and all controls/contrasts verified', flush=True)
    domains = read(D/'methylation_domain_overlap.tsv')
    assert len(domains) == 414 and domains.gene.nunique() == 23
    domain_inputs = {p: read(BASE/p) for p in domains.source_path.unique()}
    for frame in domain_inputs.values():
        if 'type' not in frame:
            frame['type'] = 'PMD'
    for row in domains.itertuples():
        source = domain_inputs[row.source_path]
        sub = source[(source.chr == row.chrom) & (source.type == row.domain)]
        shift = int(row.coordinate_assumption == 'one_based_inclusive')
        sub = sub[(sub.end > row.start) & (sub.start-shift < row.end)]
        # Independent per-base boolean union rather than the interval helper.
        mask = np.zeros(row.end-row.start, dtype=bool)
        for segment in sub.itertuples():
            left = max(row.start,int(segment.start)-shift)-row.start
            right = min(row.end,int(segment.end))-row.start
            mask[left:right] = True
        close(mask.mean(), row.domain_overlap_fraction)
    dp = domains.pivot(index=['gene','day','domain'],columns='coordinate_assumption',values='domain_overlap_fraction')
    coordinate_delta = abs(dp.one_based_inclusive-dp.zero_based_half_open_sensitivity).max()
    assert coordinate_delta <= .0005 + 1e-12
    out.mkdir(parents=True)
    windows.to_csv(out/'histone_window_and_H3_sensitivity.tsv',sep='\t',index=False)
    pd.DataFrame(checks).to_csv(out/'base_resolution_checks.tsv',sep='\t',index=False)
    pd.DataFrame(sorted(chroms.items()), columns=['chrom','size']).to_csv(out/'track_chromosome_sizes.tsv',sep='\t',index=False)
    stat = []
    for mark, f in windows.groupby('mark'):
        stat.append(dict(mark=mark, possible_preparation_contrasts=len(f), eligible_1kb=int(f.eligible_1kb.sum()),
            eligible_5kb=int(f.eligible_5kb.sum()), both_windows_eligible=int(f.both_windows_eligible.sum()),
            same_direction_across_windows=int(f.same_direction_across_windows.sum()),
            H3_changes_direction_1kb=int(f.H3_changes_direction_1kb.sum()),
            H3_changes_direction_5kb=int(f.H3_changes_direction_5kb.sum())))
    pd.DataFrame(stat).to_csv(out/'histone_sensitivity_counts.tsv',sep='\t',index=False)
    record = dict(status='passed',utc=datetime.now(timezone.utc).isoformat(),code_sha256=sha(Path(__file__)),
        verified_sha256=hashes,archived_executed_sources_verified=4,primary_MD5_verified=True,
        fits_verified=11,retained_genes=14811,whole_family_BH_verified=True,
        histone_headers_verified=24,base_resolution_windows_verified=len(checks),
        assembly_index_verified=assembly_index_verified,histone_signals_verified=len(sig),
        domain_intersections_verified=len(domains),maximum_coordinate_overlap_change=float(coordinate_delta),
        output_sha256={p.relative_to(BASE).as_posix():sha(p) for p in sorted(out.iterdir())})
    rp.write_text(json.dumps(record,indent=2)+'\n')
    print('J3: 414 domain intersections independently verified; maximum coordinate change', coordinate_delta)
    print(pd.DataFrame(stat).to_string(index=False))


if __name__ == '__main__':
    main()
