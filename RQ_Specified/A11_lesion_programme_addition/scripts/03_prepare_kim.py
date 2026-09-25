"""After discovery reproduction, stream only declared Kim cells into paired pseudobulks."""
import argparse
import gzip
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]


def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def main():
    import numpy as np
    import pandas as pd
    ap = argparse.ArgumentParser()
    ap.add_argument('--data-root', type=Path, required=True)
    args = ap.parse_args()
    out = HERE / 'tables/test_v2'
    cache = HERE / 'cache/test_v2'
    cache.mkdir(parents=True, exist_ok=True)
    if (out / 'prepare_run.json').exists() or (cache / 'counts.tsv.gz').exists():
        raise SystemExit('Refusing to overwrite Kim pseudobulks')
    check = pd.read_csv(out / 'instrument_check.tsv', sep='\t')
    assert len(check) == 46 and check.abs_error.max() <= 1e-6
    gate = json.loads((HERE / 'tables/gates_run.json').read_text())
    config = json.loads((HERE / 'config/kim2020_test_contract.json').read_text())
    source = args.data_root / 'raw_data/GSE131907'
    paths = [source / Path(name).name for name in gate['inputs']]
    for p, meta in zip(paths, gate['inputs'].values()):
        if p.stat().st_size != meta['bytes'] or sha(p) != meta['sha256']:
            raise SystemExit('Kim input differs from original gate: ' + p.name)
    pairs = pd.read_csv(HERE / 'tables/kim_pairing.tsv', sep='\t')
    pairs = pairs[pairs.patient.isin(gate['gate2']['eligible_patients'])]
    counts = pd.read_csv(HERE / 'tables/kim_cell_counts.tsv', sep='\t').set_index('patient')
    annotation = source / 'GSE131907_Lung_Cancer_cell_annotation.txt.gz'
    a = pd.read_csv(annotation, sep='\t')
    cell_column = 'Index' if 'Index' in a else a.columns[0]
    assert a[cell_column].is_unique
    labels = {'normal':config['populations']['normal_arm']['author_labels'],
              'lesion':config['populations']['lesion_arm']['author_labels']}
    groups, cells = [], []
    for _, p in pairs.iterrows():
        for arm, sample_column, n_column in [('normal','normal_sample','normal_type2_cells'),
                                              ('lesion','tumour_sample','lesion_epithelial_cells')]:
            selected = a.loc[(a.Sample==p[sample_column]) & a.Cell_subtype.isin(labels[arm]), cell_column].tolist()
            assert len(selected) == counts.loc[p.patient,n_column]
            groups.append({'unit_id':p.patient+'_'+arm,'patient':p.patient,'arm':arm,'cells':len(selected)})
            cells.append(selected)
    all_cells = [c for x in cells for c in x]
    assert len(all_cells)==len(set(all_cells))
    matrix = source / 'GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz'
    with gzip.open(matrix,'rt') as f:
        header = f.readline().rstrip('\r\n').split('\t')
    header = [h.strip('"') for h in header]
    assert len(set(header))==len(header)
    assert set(all_cells)<=set(header)
    usecols = [header[0], *all_cells]
    blocks, names = [], []
    seen = 0
    for chunk in pd.read_csv(matrix,sep='\t',usecols=usecols,index_col=0,chunksize=256):
        arr = chunk.to_numpy()
        assert np.isfinite(arr).all() and (arr>=0).all() and (arr==np.floor(arr)).all()
        blocks.append(np.column_stack([chunk[cs].sum(axis=1).to_numpy(dtype=np.int64) for cs in cells]))
        names.extend(chunk.index.astype(str))
        seen += len(chunk)
        if seen % 4096 == 0:
            print('Kim gene rows aggregated:',seen,flush=True)
    assert len(names)==gate['kim_gene_index']['n_genes'] and len(set(names))==len(names)
    assert hashlib.sha256('\n'.join(sorted(names)).encode()).hexdigest()==gate['kim_gene_index']['sorted_names_sha256']
    pd.DataFrame(np.vstack(blocks),index=names,columns=[g['unit_id'] for g in groups]).to_csv(cache/'counts.tsv.gz',sep='\t',index_label='gene')
    pd.DataFrame(groups).to_csv(out/'units.tsv',sep='\t',index=False)
    frozen = json.loads((ROOT/'RQ_Specified/A5_A11_shared_component_contract/tables/frozen_modules.json').read_text())
    es1 = {m['name']:m['genes'] for m in json.loads((ROOT/'Research Article/epithelial_state_specificity/modules.json').read_text())['modules']}
    derived = json.loads((HERE/'tables/derived_modules.json').read_text())
    member = pd.read_csv(ROOT/'RQ_Specified/A5_A11_shared_component_contract/tables/module_membership.tsv',sep='\t')
    modules = {'lesion_specific':frozen['modules']['lesion_specific']['genes'],
               'shared_remodelling':frozen['modules']['shared_remodelling']['genes'],
               'stress_excluded':derived['lesion_specific_stress_excluded']['genes'],
               'injury_lesion_pair':member.loc[member.in_I & member.in_L,'gene'].tolist(),
               **{k:es1[k] for k in ['AT1_published_400','AT2_published_400','HALLMARK_P53_PATHWAY',
                                      'HALLMARK_HYPOXIA','HALLMARK_INFLAMMATORY_RESPONSE']}}
    orth = pd.read_csv(ROOT/'Research Article/gate2_C3_yu_lee_choi_min_2026/trials/u4_resources/strict_one_to_one_orthologs.csv')
    m2h = dict(zip(orth.mouse_symbol,orth.human_symbol)); universe=set(names)
    rows=[]
    for k, gs in modules.items():
        hs = [m2h[g] for g in gs if g in m2h and m2h[g] in universe]
        assert len(hs)==len(set(hs)) and len(hs)/len(gs)>=.7
        rows += [{'module':k,'gene':h,'source_genes':len(gs),'assayed_genes':len(hs)} for h in hs]
    pd.DataFrame(rows).to_csv(out/'human_module_genes.tsv',sep='\t',index=False)
    outputs = [cache/'counts.tsv.gz',out/'units.tsv',out/'human_module_genes.tsv',out/'instrument_check.tsv']
    record={'completed_utc':datetime.now(timezone.utc).isoformat(),
            'preregistration_commit':'eb5317e','code_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
            'script_sha256':sha(Path(__file__)),'config_sha256':sha(HERE/'config/kim2020_test_contract.json'),
            'inputs':{str(p):sha(p) for p in paths},'outputs':{str(p.relative_to(HERE)):sha(p) for p in outputs},
            'gene_rows':len(names),'cells':len(all_cells),'patients':len(pairs),
            'pandas':pd.__version__,'numpy':np.__version__}
    (out/'prepare_run.json').write_text(json.dumps(record,indent=2)+'\n')
    print('Kim prepared:',len(pairs),'patients;',len(names),'genes',flush=True)


if __name__=='__main__':
    main()
