"""Rebuild reusable donor/animal counts without changing frozen trial outputs.

Run with the repository's recovered Python launcher. All cache files are local,
regenerable and ignored. Human matrices are streamed once per cohort, retaining
the original compartment membership and a separate subtype decomposition.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SITE = ROOT / '.venv-x64/Lib/site-packages'
if SITE.is_dir():
    sys.path.insert(0, str(SITE))
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

TRIALS = ROOT / 'Research Article/gate1_01_niethamer_2025/trials'
sys.path.insert(0, str(TRIALS))
import gsea_utils as gu
import g2_gsea_ipf as g2

CACHE = HERE / 'cache'
CACHE.mkdir(exist_ok=True)


class Inputs:
    def __init__(self):
        self.inputs = []
        self.facts = {}

    def add_input(self, path):
        path = Path(path)
        item = {'path': path.relative_to(ROOT).as_posix(), 'bytes': path.stat().st_size,
                'mtime_ns': path.stat().st_mtime_ns}
        # Large raw inputs are identified by stat and count-total parity with the
        # frozen outputs; the historical trial did not provide raw-file hashes.
        if path.stat().st_size < 10_000_000:
            item['sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
        self.inputs.append(item)

    def set(self, key, value):
        self.facts[key] = value


def save(name, counts, genes, units, record):
    genes = np.asarray(genes, dtype=str)
    uniq, inv = np.unique(genes, return_inverse=True)
    collapsed = np.zeros((len(units), len(uniq)), dtype=np.int64)
    for i in range(len(units)):
        np.add.at(collapsed[i], inv, counts[i])
    assert np.all(collapsed >= 0)
    assert np.array_equal(collapsed.sum(1), np.asarray(counts).sum(1))
    units = units.reset_index(drop=True).copy()
    units['unit_id'] = ['u' + str(i) for i in range(len(units))]
    np.savez_compressed(CACHE / f'{name}.npz', counts=collapsed, genes=uniq)
    units.to_csv(CACHE / f'{name}_units.csv', index=False)
    pd.DataFrame(collapsed.T, index=uniq, columns=units.unit_id).to_csv(
        CACHE / f'{name}_counts.csv.gz', compression='gzip', index_label='gene')
    record['n_units'] = len(units)
    record['n_genes'] = len(uniq)
    record['total_counts'] = int(collapsed.sum())
    record['completed_utc'] = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())
    (HERE / f'{name}_preparation.json').write_text(json.dumps(record, indent=2) + '\n')
    print(name, record['n_units'], record['n_genes'], record['total_counts'], flush=True)


def human(name):
    rec = Inputs()
    loader = g2.cohort_136831 if name == 'GSE136831' else g2.cohort_135893
    meta, genes, raw, key = loader(rec)
    meta['comp'] = None
    for comp, labels in g2.COMPARTMENTS.items():
        meta.loc[meta.label.isin(labels[key]), 'comp'] = comp
    meta['disease'] = meta.disease.replace({'Control': 'control'})
    selected = meta.comp.notna() & meta.disease.isin(['IPF', 'control'])
    totals = meta[selected].groupby(['comp', 'donor', 'disease']).size().rename('cells').reset_index()
    totals = totals[totals.cells >= 50].reset_index(drop=True)
    keys = set(zip(totals.comp, totals.donor))
    selected &= pd.Series(list(zip(meta.comp, meta.donor)), index=meta.index).isin(keys)
    subunits = meta[selected].groupby(['comp', 'donor', 'disease', 'label']).size().rename('cells').reset_index()
    subkey = {(r.comp, r.donor, r.label): i for i, r in subunits.iterrows()}
    goc = np.full(len(meta), -1, dtype=np.int64)
    chosen = np.flatnonzero(selected.to_numpy())
    goc[chosen] = [subkey[(r.comp, r.donor, r.label)] for r in meta.loc[chosen].itertuples()]
    print(name, 'streaming', len(subunits), 'subtype-donor groups', flush=True)
    subcounts, seen = gu.pseudobulk_from_mtx(raw, goc, len(subunits), chunk_rows=2_000_000)
    subcounts = subcounts.T
    record = {'inputs': rec.inputs, 'nnz_read': seen, 'preparation_code_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    save(name + '_subtypes', subcounts, genes, subunits, dict(record))
    totalkey = {(r.comp, r.donor): i for i, r in totals.iterrows()}
    counts = np.zeros((len(totals), len(genes)), dtype=np.int64)
    for i, r in subunits.iterrows():
        counts[totalkey[(r.comp, r.donor)]] += subcounts[i]
    frozen = pd.read_csv(TRIALS / 'g2_gsea_ipf/g2_units.csv')
    frozen = frozen[frozen.cohort == name].set_index(['comp', 'donor'])
    observed = totals.set_index(['comp', 'donor'])
    assert set(observed.index) == set(frozen.index)
    expected = frozen.loc[observed.index]
    assert np.array_equal(observed.cells.to_numpy(), expected.cells.to_numpy())
    assert np.array_equal(counts.sum(1), expected['counts'].to_numpy())
    record['historical_cell_and_library_total_parity'] = True
    save(name, counts, genes, totals, record)


def mouse():
    base = ROOT / 'Research Article/gate1_01_niethamer_2025/GSE262927'
    src = base / 'processed/final_clustered.h5ad'
    mypath = base / 'myeloid_focus/tables/myeloid_cell_metadata.csv'
    mpath = base / 'myeloid_focus/batch_sensitivity/tables/sample_infection_round.csv'
    rec = Inputs()
    for p in [src, mypath, mpath]:
        rec.add_input(p)
    a = ad.read_h5ad(src, backed='r')
    obs = a.obs.copy()
    genes = a.var['gene_symbol'].astype(str).to_numpy() if 'gene_symbol' in a.var else a.var_names.to_numpy()
    a.file.close()
    my = pd.read_csv(mypath)
    my = my[my.label != 'other lineage'].copy()
    idx = obs.index.get_indexer(my.cell_id)
    assert np.all(idx >= 0)
    my['obs_index'] = idx
    annotated = obs.has_author_metadata.astype(str).isin(['True', 'true']).to_numpy()
    my = my[annotated[idx]]
    units = my.groupby(['sample_id', 'day', 'label']).size().rename('cells').reset_index()
    key = {(r.sample_id, r.label): i for i, r in units.iterrows()}
    goc = np.full(len(obs), -1, dtype=np.int64)
    goc[my.obs_index] = [key[(r.sample_id, r.label)] for r in my.itertuples()]
    counts = gu.pseudobulk_from_h5ad_counts(src, goc, len(units))
    md = pd.read_csv(mpath)[['sample_id', 'genotype', 'sex', 'round']]
    units = units.merge(md, on='sample_id', validate='many_to_one')
    record = {'inputs': rec.inputs, 'preparation_code_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    save('mouse_subtypes', counts, genes, units, dict(record))
    totalunits = units.groupby(['sample_id', 'day', 'genotype', 'sex', 'round'], sort=True).cells.sum().reset_index()
    totalkey = {r.sample_id: i for i, r in totalunits.iterrows()}
    totalcounts = np.zeros((len(totalunits), len(genes)), dtype=np.int64)
    for i, r in units.iterrows():
        totalcounts[totalkey[r.sample_id]] += counts[i]
    frozen = pd.read_csv(TRIALS / 'g1_gsea_by_phase/g1_units.csv')
    frozen = frozen[frozen.compartment == 'myeloid'].set_index('animal').loc[totalunits.sample_id]
    assert np.array_equal(totalunits.cells.to_numpy(), frozen.cells.to_numpy())
    assert np.array_equal(totalcounts.sum(1), frozen['counts'].to_numpy())
    record['historical_cell_and_library_total_parity'] = True
    save('mouse_myeloid', totalcounts, genes, totalunits, record)
    amac = units.label == 'aMAC'
    au = units[amac].reset_index(drop=True)
    ac = counts[amac]
    frozenw = pd.read_csv(TRIALS / 'w1_amac_pseudobulk_de/w1_units.csv').set_index('animal')
    assert np.array_equal(au.cells.to_numpy(), frozenw.loc[au.sample_id].aMAC_cells.to_numpy())
    save('mouse_amac', ac, genes, au, {'inputs': rec.inputs, 'historical_cell_parity': True})


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('cohort', choices=['mouse', 'GSE136831', 'GSE135893'])
    args = parser.parse_args()
    mouse() if args.cohort == 'mouse' else human(args.cohort)
