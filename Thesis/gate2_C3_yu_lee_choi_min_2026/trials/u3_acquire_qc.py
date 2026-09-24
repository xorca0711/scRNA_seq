"""Checkpointed acquisition and QC; provisional coverage, no treatment inference."""
from __future__ import annotations
import argparse
import csv
from datetime import datetime, timedelta, timezone
import gzip
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import sys
import time
import urllib.request

for key in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMBA_NUM_THREADS']:
    os.environ.setdefault(key, '2')
PAPER = Path(__file__).resolve().parents[1]
ROOT = PAPER.parents[1]
sys.path.insert(0, str(ROOT))
from analysis.lib.provenance import archive_existing_record, code_identity, sha256_file, write_json_atomic


def read_lines(path):
    with gzip.open(path, 'rt', encoding='utf-8') as handle:
        return [line.rstrip('\n\r').split('\t') for line in handle if line.strip()]


def read_counts(files):
    import numpy as np
    from scipy.io import mmread
    features = read_lines(files['features'])
    barcodes = [x[0] for x in read_lines(files['barcodes'])]
    ids = [x[0] for x in features]
    if len(set(ids)) != len(ids) or len(set(barcodes)) != len(barcodes):
        raise ValueError('Nonunique feature IDs or within-sample barcodes')
    if any(len(x) > 2 and x[2] != 'Gene Expression' for x in features):
        raise ValueError('Mixed feature modalities require separate handling')
    with gzip.open(files['matrix'], 'rb') as handle:
        coo = mmread(handle).tocoo()
    if coo.shape != (len(features), len(barcodes)):
        raise ValueError('Matrix dimension mismatch')
    if not np.all(np.isfinite(coo.data)) or np.any(coo.data < 0) or np.any(coo.data != np.floor(coo.data)):
        raise ValueError('Matrix is not finite nonnegative integer counts')
    if coo.data.size and coo.data.max() > np.iinfo(np.int64).max:
        raise ValueError('Count exceeds int64 capacity')
    coo.data = coo.data.astype(np.int64)
    x = coo.T.tocsr()
    duplicates = int(coo.nnz-x.nnz)
    x.eliminate_zeros()
    return x, features, barcodes, duplicates


def download(item, cache, spec, ledger):
    path = cache / item['filename']
    if path.name != item['filename']:
        raise ValueError('Unsafe manifest filename')
    before = time.monotonic()
    new_bytes = 0
    if not path.exists():
        partial = path.with_suffix(path.suffix+'.part')
        try:
            req = urllib.request.Request(item['url'], headers={'User-Agent': 'Public-research-reanalysis/1.0'})
            with urllib.request.urlopen(req, timeout=60) as response, partial.open('wb') as dest:
                expected = response.headers.get('Content-Length')
                while chunk := response.read(1024*1024):
                    new_bytes += len(chunk)
                    if new_bytes > spec['maximum_bytes_per_file'] or ledger['downloaded_bytes']+new_bytes > spec['maximum_new_download_bytes']:
                        raise ValueError('Frozen acquisition cap exceeded')
                    dest.write(chunk)
                if expected is not None and new_bytes != int(expected):
                    raise ValueError('Incomplete HTTP response')
            partial.replace(path)
        finally:
            if partial.exists():
                partial.unlink()
    if path.stat().st_size > spec['maximum_bytes_per_file']:
        raise ValueError('Cached file exceeds frozen cap')
    sha = sha256_file(path)
    if item.get('previous_sha256') and item['previous_sha256'] != sha:
        raise ValueError('Cached preparation input hash changed')
    ledger['downloaded_bytes'] += new_bytes
    return path, {'url': item['url'], 'filename': path.name, 'bytes': path.stat().st_size,
                  'sha256': sha, 'downloaded_bytes': new_bytes, 'seconds': round(time.monotonic()-before, 3)}


def process_sample(sample, files, spec, out, cache):
    import numpy as np
    import pandas as pd
    import anndata as ad
    start = time.monotonic()
    x, features, barcodes, duplicates = read_counts(files)
    symbols = np.array([f[1] if len(f)>1 else f[0] for f in features])
    total = np.asarray(x.sum(axis=1)).ravel()
    genes = np.asarray((x>0).sum(axis=1)).ravel()
    mt = np.char.startswith(np.char.lower(symbols.astype(str)), 'mt-')
    mt_counts = np.asarray(x[:, mt].sum(axis=1)).ravel()
    pct_mt = np.divide(100.0*mt_counts, total, out=np.zeros(len(total), dtype=float), where=total>0)
    qc = spec['qc']
    keep = (total>0) & (genes>=qc['min_genes']) & (genes<=qc['max_genes']) & (pct_mt<=qc['max_percent_mt'])
    obs = pd.DataFrame({'gsm': sample['gsm'], 'n_counts': total, 'n_genes': genes,
                        'percent_mt': pct_mt, 'qc_pass': keep}, index=[sample['gsm']+':'+b for b in barcodes])
    memberships = []
    counts = []
    for label, panel in spec['candidate_panels'].items():
        # Symbols can map to multiple feature IDs; test detection once per symbol.
        detected = np.zeros(len(total), dtype=int)
        measured = [g for g in panel if np.any(symbols == g)]
        for gene in measured:
            detected += np.asarray(x[:, symbols == gene].sum(axis=1)).ravel()>0
        member = detected >= spec['minimum_detected_panel_genes']
        obs['panel_detected_'+label] = detected
        memberships.append(member)
        counts.append({'gsm': sample['gsm'], 'candidate_compartment': label,
                       'panel_genes': len(panel), 'assayed_panel_genes': len(measured),
                       'qc_pass_panel_positive_cells': int((member & keep).sum())})
    memberships = np.column_stack(memberships)
    positive_n = memberships.sum(axis=1)
    labels = np.full(len(total), 'unassigned', dtype=object)
    labels[positive_n > 1] = 'ambiguous_multiple_panels'
    for j, label in enumerate(spec['candidate_panels']):
        labels[(positive_n == 1) & memberships[:, j]] = label
        counts[j]['qc_pass_exclusive_candidate_cells'] = int((keep & (labels == label)).sum())
    obs['candidate_compartment'] = labels
    for row in counts:
        row['provisional_only'] = True
    var = pd.DataFrame({'gene_symbol': symbols}, index=[f[0] for f in features])
    data = ad.AnnData(x, obs=obs, var=var)
    data.uns['annotation_status'] = 'provisional_marker_coverage_not_validated_cell_types'
    data.write_h5ad(cache / (sample['gsm']+'_qc.h5ad'), compression='gzip')
    summary = {'gsm': sample['gsm'], 'treatment': sample['treatment'], 'endpoint': sample['endpoint'],
               'analysis_unit': sample['analysis_unit'], 'input_cells': len(total), 'qc_pass_cells': int(keep.sum()),
               'qc_excluded_cells': int((~keep).sum()), 'features': x.shape[1], 'mt_features': int(mt.sum()),
               'median_counts': float(np.median(total)), 'median_genes': float(np.median(genes)),
               'median_percent_mt': float(np.median(pct_mt)), 'duplicate_coordinates_collapsed': duplicates,
               'ambiguous_candidate_cells_after_qc': int((keep & (positive_n>1)).sum()),
               'processing_seconds': round(time.monotonic()-start, 3)}
    (out / (sample['gsm']+'_qc.json')).write_text(json.dumps(summary, indent=2)+'\n', encoding='utf-8')
    return summary, counts


def write_table(path, rows):
    if rows:
        with path.open('w', encoding='utf-8', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--spec', type=Path, default=PAPER/'trials/u3_run_spec.json')
    args = parser.parse_args()
    contract = json.loads((PAPER/'analysis_contract.json').read_text(encoding='utf-8'))
    if contract['full_run_authorized'] is not True:
        raise RuntimeError('Owner execution authorization missing')
    spec = json.loads(args.spec.read_text(encoding='utf-8'))
    out = PAPER/'trials/u3_acquire_qc'
    cache = PAPER/'cache/count_inputs'
    processed = PAPER/'cache/u3_qc'
    for directory in [out, cache, processed]: directory.mkdir(parents=True, exist_ok=True)
    status_path = out/'run_record.json'
    archive_existing_record(status_path)
    fingerprint = sha256_file(args.spec)+sha256_file(Path(__file__))
    started = time.monotonic()
    state = {'started_utc': datetime.now(timezone.utc).isoformat(), 'pid': os.getpid(),
             'stage': 'U3_acquisition_QC', 'status': 'running', 'planned_samples': len(spec['samples']),
             'completed_samples': 0, 'downloaded_bytes': 0, 'biological_treatment_effects_estimated': False,
             'spec_sha256': sha256_file(args.spec), 'code': code_identity(ROOT, __file__), 'inputs': [],
             'packages': {p: importlib.metadata.version(p) for p in ['numpy','scipy','pandas','anndata','h5py']}}
    write_json_atomic(status_path, state)
    summaries, coverage = [], []
    try:
        for sample in spec['samples']:
            state['current_sample'] = sample['gsm']
            write_json_atomic(status_path, state)
            checkpoint = out/(sample['gsm']+'_checkpoint.json')
            files = {}
            for item in sample['files']:
                files[item['kind']], receipt = download(item, cache, spec, state)
                state['inputs'].append({'gsm': sample['gsm'], **receipt})
            input_hash = hashlib.sha256(json.dumps({k:sha256_file(v) for k,v in files.items()}, sort_keys=True).encode()).hexdigest()
            existing = json.loads(checkpoint.read_text()) if checkpoint.exists() else {}
            if existing.get('fingerprint') == fingerprint and existing.get('input_hash') == input_hash and (processed/(sample['gsm']+'_qc.h5ad')).exists():
                summary, rows = existing['summary'], existing['coverage']
            else:
                summary, rows = process_sample(sample, files, spec, out, processed)
                write_json_atomic(checkpoint, {'fingerprint': fingerprint, 'input_hash': input_hash,
                                              'summary': summary, 'coverage': rows})
            summaries.append(summary)
            coverage.extend(rows)
            state['completed_samples'] = len(summaries)
            state['elapsed_seconds'] = round(time.monotonic()-started, 2)
            # This estimate concerns this bounded U3 job only, never the full biological pipeline.
            new_files = [r for r in state['inputs'] if r['downloaded_bytes']>0]
            if new_files:
                bps = sum(r['downloaded_bytes'] for r in new_files)/max(sum(r['seconds'] for r in new_files), 0.001)
                mean_bytes = sum(r['bytes'] for r in state['inputs'])/len(summaries)
                remaining = len(spec['samples'])-len(summaries)
                estimate = remaining*(mean_bytes/bps + sum(s['processing_seconds'] for s in summaries)/len(summaries))
                state['estimated_U3_completion_utc'] = (datetime.now(timezone.utc)+timedelta(seconds=estimate)).isoformat()
                state['estimate_scope'] = 'remaining acquisition/QC only; network and library size may vary'
            write_table(out/'sample_qc.csv', summaries)
            write_table(out/'candidate_compartment_coverage.csv', coverage)
            write_json_atomic(status_path, state)
            print(json.dumps({'gsm': sample['gsm'], 'completed': len(summaries), 'planned': len(spec['samples']),
                              'qc_pass_cells': summary['qc_pass_cells'], 'elapsed_seconds': state['elapsed_seconds'],
                              'estimated_U3_completion_utc': state.get('estimated_U3_completion_utc')}), flush=True)
        state['status'] = 'completed_QC_annotation_review_required'
    except Exception as exc:
        state['status'] = 'failed'
        state['error'] = type(exc).__name__+': '+str(exc)
        raise
    finally:
        state['updated_utc'] = datetime.now(timezone.utc).isoformat()
        state['elapsed_seconds'] = round(time.monotonic()-started, 2)
        write_json_atomic(status_path, state)


if __name__ == '__main__':
    main()
