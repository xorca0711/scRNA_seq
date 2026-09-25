"""Prepare bounded display embeddings and donor PCA without refitting inference."""
from pathlib import Path
import os
for key in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMBA_NUM_THREADS']:
    os.environ.setdefault(key, '2')
import sys
import json
import time
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[2]
PAPER = ROOT / 'Research Article/gate2_C3_yu_lee_choi_min_2026'
sys.path.insert(0, str(ROOT))
from analysis.lib.provenance import sha256_file, code_identity, write_json_atomic
sys.path.insert(0, str(PAPER / 'trials'))
from u5_human_aggregate import TRIAD


def main():
    import numpy as np
    import pandas as pd
    import h5py
    from scipy.sparse import csr_matrix
    from sklearn.decomposition import PCA

    out = ROOT / 'analysis/figures/rq/il1b_context'
    cache = out / 'cache'
    out.mkdir(parents=True, exist_ok=True)
    cache.mkdir(exist_ok=True)
    record = out / 'preparation_run_record.json'
    assert not record.exists(), 'Existing preparation is immutable; reuse its saved coordinates.'
    specification = dict(
        frozen_utc=datetime.now(timezone.utc).isoformat(), scope='Post-analysis display diagnostics; no new DE, enrichment or label fitting',
        cohort='GSE308103', seed=20260924,
        umap=dict(input='Saved 30-dimensional frozen-reference query latent vectors',
                  display_sampling='At most 75 cells per patient, histology and existing source-profile broad label; repeated libraries pooled before sampling',
                  n_neighbors=30, min_dist=0.35, metric='euclidean', n_components=2, n_epochs=300, n_jobs=1,
                  interpretation='Reference-conditioned display geometry; not abundance, lineage or independent annotation validation'),
        pca=dict(unit='patient-histology broad AT2-like pseudobulk', primary_uncertainty=0.2, cell_floor=50,
                 normalization='log2(CPM + 1), full assay library totals; diagnostic normalization, not the TMM inferential model',
                 genes='CPM >= 1 in at least 5 aggregates; top 2000 by variance without using histology labels',
                 scaling='Gene centering, no variance scaling, no batch correction', n_components=5),
        source_counts='Existing all-QC per-patient source fractions; no display sampling used in quantitative summaries',
        program_contrast='Existing primary uncertainty 0.2 / 50-cell / broad AT2 paired values; no new hypothesis tests',
        pathway_view='All eligible broad-receiver LUAD-normal rows; primary and fixed001 sensitivity side by side',
        ligand_target_view='All eligible focused_triad, broad AT2/fibroblast, LUAD-normal, up-target candidates; no top-hit selection')
    write_json_atomic(out / 'specification.json', specification)
    start = time.monotonic()
    state = dict(status='running', started_utc=specification['frozen_utc'], code=code_identity(ROOT, __file__))
    write_json_atomic(record, state)
    inputs = []

    def input_path(rel):
        path = PAPER / rel
        inputs.append(dict(path=rel, bytes=path.stat().st_size, sha256=sha256_file(path)))
        return path

    try:
        jobs = json.loads(input_path('trials/u5_human_full/acquisition_specification.json').read_text())['jobs']
        frames = []
        for job in jobs:
            gsm = job['gsm']
            frame = pd.read_csv(input_path(f'cache/u5_human_full/{gsm}/cell_annotations.csv.gz'))
            frame['library_row'] = np.arange(len(frame))
            assert frame.gsm.eq(gsm).all() and frame.patient.astype(str).eq(str(job['patient'])).all()
            frame['broad'] = 'Other assigned'
            frame.loc[frame.ann_level_1.eq('Epithelial'), 'broad'] = 'Other epithelial'
            frame.loc[frame.ann_level_2.eq('Myeloid'), 'broad'] = 'Other myeloid'
            for comp, labels in TRIAD.items():
                frame.loc[frame.ann_finest_level.isin(labels), 'broad'] = {'AT2': 'AT2-like', 'fibroblasts': 'Fibroblasts', 'macrophages': 'Macrophages'}[comp]
            frame.loc[frame.ann_finest_level_uncertainty.gt(.2), 'broad'] = 'Unassigned'
            frames.append(frame)
        full = pd.concat(frames, ignore_index=True)
        assert full.cell_id.is_unique and len(full) == 555480 and full.patient.nunique() == 23
        rng = np.random.default_rng(specification['seed'])
        selected = np.sort(np.concatenate([rng.choice(ix, min(75, len(ix)), replace=False)
                                          for ix in full.groupby(['patient', 'histology', 'broad'], sort=True).groups.values()]))
        sample = full.loc[selected].copy().reset_index(drop=True)
        coverage = full.groupby(['patient', 'histology', 'broad']).size().rename('all_QC_cells').to_frame()
        coverage['display_cells'] = sample.groupby(['patient', 'histology', 'broad']).size()
        coverage.reset_index().to_csv(out / 'embedding_sampling_coverage.csv', index=False)
        assert coverage.display_cells.le(75).all() and coverage.display_cells.gt(0).all()
        # Full-population QC summaries, not estimates from the balanced display sample.
        full.groupby(['patient', 'histology', 'broad']).agg(
            cells=('cell_id', 'size'), median_library_size=('full_library_size', 'median'),
            median_detected_genes=('n_genes_by_counts', 'median'), median_mito_percent=('pct_counts_mt', 'median'),
            median_uncertainty=('ann_finest_level_uncertainty', 'median')).reset_index().to_csv(out / 'donor_QC_summary.csv', index=False)
        latent = np.empty((len(sample), 30), dtype=np.float32)
        sample['IL1B_count'] = 0
        for gsm, ix in sample.groupby('gsm', sort=True).groups.items():
            rows = sample.loc[ix, 'library_row'].to_numpy()
            z = np.load(input_path(f'cache/u5_human_full/{gsm}/query_latent.npz'))['latent']
            assert z.shape == (int(full.gsm.eq(gsm).sum()), 30)
            latent[ix] = z[rows]
            with h5py.File(input_path(f'cache/u5_human_full/{gsm}/annotated_panel.h5ad'), 'r') as h:
                # Check cell order explicitly because saved latent rows have positional IDs.
                cell_ids = h['obs'][h['obs'].attrs['_index']].asstr()[:]
                assert np.array_equal(cell_ids[rows], sample.loc[ix, 'cell_id'])
                genes = h['var'][h['var'].attrs['_index']].asstr()[:]
                gene_index = list(genes).index('IL1B')
                g = h['X']
                x = csr_matrix((g['data'][:], g['indices'][:], g['indptr'][:]), shape=tuple(g.attrs['shape']))
                sample.loc[ix, 'IL1B_count'] = x[rows, gene_index].toarray().ravel()
        sample['IL1B_log1p_normalized_10000'] = np.log1p(sample.IL1B_count * 10000 / sample.full_library_size)
        np.savez_compressed(cache / 'display_latent.npz', latent=latent)
        sample.to_csv(cache / 'display_cells_before_embedding.csv.gz', index=False, compression='gzip')
        state.update(stage='UMAP', display_cells=len(sample))
        write_json_atomic(record, state)
        print(f'Prepared {len(sample):,} display cells from {len(full):,} QC nuclei; fitting saved-latent UMAP.', flush=True)
        # Import only after preparing the diagnostic inputs; no training is performed.
        import umap
        params = specification['umap']
        model = umap.UMAP(**{key: params[key] for key in ['n_neighbors', 'min_dist', 'metric', 'n_components', 'n_epochs', 'n_jobs']},
                          random_state=specification['seed'], transform_seed=specification['seed'])
        coords = model.fit_transform(latent)
        assert coords.shape == (len(sample), 2) and np.isfinite(coords).all()
        sample['UMAP1'], sample['UMAP2'] = coords.T
        # Only plotted fields; full source matrices remain in ignored caches.
        sample[['cell_id', 'gsm', 'patient', 'histology', 'broad', 'ann_finest_level_uncertainty',
                'IL1B_count', 'IL1B_log1p_normalized_10000', 'UMAP1', 'UMAP2']].to_csv(out / 'umap_display_values.csv', index=False)
        print('UMAP coordinates saved; preparing patient-level PCA.', flush=True)

        meta = pd.read_csv(input_path('trials/u5_human_niche/unc20_pooled_triad_units.csv'))
        meta = meta[(meta.comp == 'AT2') & (meta.label == '__broad__') & (meta.cells >= 50)].copy()
        counts_path = input_path('cache/u5_human_niche/unc20_pooled_triad_counts.csv.gz')
        gene_col = pd.read_csv(counts_path, nrows=0).columns[0]
        counts = pd.read_csv(counts_path, usecols=[gene_col] + meta.unit_id.tolist(), index_col=0)
        counts = counts.loc[:, meta.unit_id]
        assert np.array_equal(counts.sum(axis=0).to_numpy(), meta.full_library_sum.to_numpy())
        cpm = counts.to_numpy(dtype=float) * 1e6 / meta.full_library_sum.to_numpy()[None, :]
        keep = (cpm >= 1).sum(axis=1) >= 5
        log = np.log2(cpm[keep] + 1)
        genes = counts.index.to_numpy()[keep]
        variance = log.var(axis=1, ddof=1)
        order = np.argsort(-variance, kind='stable')[:2000]
        pca = PCA(n_components=5, svd_solver='full')
        scores = pca.fit_transform(log[order].T)
        for i in range(5):
            meta[f'PC{i+1}'] = scores[:, i]
        meta.to_csv(out / 'AT2_pseudobulk_PCA_values.csv', index=False)
        pd.DataFrame(dict(gene=genes[order], input_variance=variance[order],
                          PC1_loading=pca.components_[0], PC2_loading=pca.components_[1])).to_csv(out / 'AT2_PCA_features.csv', index=False)
        write_json_atomic(out / 'PCA_summary.json', dict(units=len(meta), patients=int(meta.patient.nunique()),
                          expressed_genes=int(keep.sum()), selected_genes=len(order), explained_variance_ratio=pca.explained_variance_ratio_.tolist()))
        pd.DataFrame(inputs).drop_duplicates('path').to_csv(out / 'preparation_inputs.csv', index=False)
        import importlib.metadata as im
        versions = {pkg: im.version(pkg) for pkg in ['numpy', 'pandas', 'scipy', 'scikit-learn', 'umap-learn', 'h5py']}
        state.update(status='completed', stage='diagnostic_inputs_ready', all_QC_cells=len(full), display_cells=len(sample),
                     sampled_strata=len(coverage), PCA_units=len(meta), patients=23, versions=versions,
                     no_new_statistical_tests=True, no_annotation_changes=True)
    except Exception as exc:
        state.update(status='failed', error=repr(exc))
        raise
    finally:
        state.update(elapsed_seconds=round(time.monotonic()-start, 2), updated_utc=datetime.now(timezone.utc).isoformat())
        write_json_atomic(record, state)


if __name__ == '__main__':
    main()
