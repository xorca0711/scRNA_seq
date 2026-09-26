"""Prepare bounded-memory count inputs and eligibility; do not estimate programme effects."""
from __future__ import annotations
import argparse, gzip, hashlib, json, re
from datetime import datetime, timezone
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix, coo_matrix

BASE = Path(__file__).resolve().parents[1]
ROOT = BASE.parents[1]
CACHE = BASE / 'cache/continuation_v1'
OUT = BASE / 'tables/pilot_v1'
WORK = BASE / 'processed/pilot_v1'

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(4*1024*1024),b''): h.update(b)
    return h.hexdigest()

def write_json(p,x):
    p.write_text(json.dumps(x,indent=2)+'\n',encoding='utf-8')

def choose(meta, cfg):
    counts=pd.crosstab(meta.unit,meta.state).reindex(columns=['start','intermediate','destination'],fill_value=0)
    counts['eligible']=counts.ge(cfg['cell_floor']).all(axis=1)
    return meta[meta.unit.isin(counts.index[counts.eligible])].reset_index(drop=True),counts

def dense_text(path,meta,header_has_gene,output,raw_check=None):
    genes=[]; libraries=np.zeros(len(meta),dtype=np.int64); check_n=0
    with gzip.open(path,'rt') as f,output.open('xb') as out:
        header=f.readline().rstrip('\n\r').split('\t')
        if header_has_gene: header=header[1:]
        assert len(set(header))==len(header)
        lookup={s:i for i,s in enumerate(header)}
        inds=np.array([lookup[s] for s in meta.cell_id])
        if raw_check:
            raw, raw_genes, raw_cols=raw_check
            rid={g:i for i,g in enumerate(raw_genes)}
            check_indices=np.array([lookup[s] for s in raw_cols])
        for line in f:
            gene, values=line.rstrip('\n\r').split('\t',1)
            # Integer-formatted UCSC exports need not parse billions of zeros as
            # floating point. Validate the entire row; long/non-integer tokens
            # use the original parser so accepted numerical inputs do not change.
            encoded=values.encode('ascii')
            plain_integer=not encoded.translate(None,b'0123456789\t') and not re.search(rb'\d{10}',encoded)
            v=np.fromstring(values,sep='\t',dtype=np.int64 if plain_integer else np.float64)
            assert len(v)==len(header) and np.isfinite(v).all() and (v>=0).all() and (v==np.floor(v)).all(),gene
            assert v.max()<np.iinfo(np.int32).max
            if raw_check and gene in rid:
                np.testing.assert_array_equal(v[check_indices],raw[rid[gene],:].toarray().ravel())
                check_n+=len(check_indices)
            row=v[inds].astype(np.int32); row.tofile(out); libraries+=row
            genes.append(gene)
            if len(genes)%5000==0: print(f'{output.stem}: {len(genes)} genes prepared',flush=True)
    return genes,libraries,check_n

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--legacy-root',type=Path,required=True)
    p.add_argument('--a5-source-root',type=Path,required=True)
    a=p.parse_args()
    cfg=json.loads((BASE/'config/pilot_v1.json').read_text())
    if WORK.exists() or any((OUT/name).exists() for name in ['preparation.json','eligibility.tsv','depth_coverage.tsv','ortholog_universe.tsv']):
        raise SystemExit('Refusing to overwrite pilot preparation')
    OUT.mkdir(parents=True,exist_ok=True); WORK.mkdir(parents=True)
    inputs={}; datasets={}; coverage=[]
    def bind(p):
        inputs[str(p)]=dict(bytes=p.stat().st_size,sha256=sha(p)); return p
    # D1: source-QC cells, fixed interval, independent mouse IDs already audited.
    old=a.legacy_root/'RQ_Specified/A0_conserved_epithelial_transition_program'
    d1meta=bind(old/'cache/sources/GSE141259_HighResolution_cellinfo.csv.gz')
    d1genes=bind(a.a5_source_root/'RQ_Specified/A5_developmental_programme_reuse/cache/GSE141259_HighResolution_genes.txt.gz')
    d1matrix=bind(ROOT/'RQ_Specified/A5_developmental_programme_reuse/cache/GSE141259_HighResolution_rawcounts.mtx.gz')
    d1bars=bind(ROOT/'RQ_Specified/A5_developmental_programme_reuse/cache/GSE141259_HighResolution_barcodes.txt.gz')
    original=pd.read_csv(d1meta,sep='\t'); original['day']=original.time_point.str.extract(r'(\d+)')[0].astype(int)
    labels={cfg['D1'][s]:s for s in ['start','intermediate','destination','control']}
    d=original[original.day.between(10,15)&~original.sample_id.str.startswith('NC-')&original.cell_type.isin(labels)].copy()
    d=d.rename(columns={'cell_barcode':'cell_id','sample_id':'unit','cell_type':'author_state'})
    d['state']=d.author_state.map(labels); d=d[['cell_id','unit','state','author_state','day']]
    d,cc=choose(d,cfg); cc['role']='D1'; coverage.append(cc.reset_index())
    genes=[s.strip() for s in gzip.open(d1genes,'rt')]; bars=[s.strip() for s in gzip.open(d1bars,'rt')]
    assert len(set(genes))==len(genes) and len(set(bars))==len(bars)
    bi={s:i for i,s in enumerate(bars)}; ix=np.array([bi[s] for s in d.cell_id])
    mapping=np.full(len(bars),-1,dtype=np.int64); mapping[ix]=np.arange(len(ix))
    out=WORK/'D1_counts.bin'; kept_rows=[]; kept_columns=[]; kept_values=[]
    with gzip.open(d1matrix,'rt') as f:
        line=f.readline()
        while line.startswith('%'): line=f.readline()
        nr,nc,nnz=map(int,line.split()); transpose=(nr,nc)==(len(bars),len(genes))
        assert (nr,nc) in [(len(bars),len(genes)),(len(genes),len(bars))]
        seen=0
        while lines:=f.readlines(4*1024*1024):
            z=np.fromstring(''.join(lines),sep=' ').reshape(-1,3)
            assert np.isfinite(z).all() and (z==np.floor(z)).all()
            rr=z[:,0].astype(np.int64)-1; co=z[:,1].astype(np.int64)-1
            if transpose: rr,co=co,rr
            keep=mapping[co]>=0
            kept_rows.append(rr[keep].astype(np.int32))
            kept_columns.append(mapping[co[keep]].astype(np.int32))
            kept_values.append(z[keep,2].astype(np.int32))
            seen+=len(z)
            if seen%5_000_000<len(z):print(f'D1: {seen} source entries inspected',flush=True)
        assert seen==nnz
    x=coo_matrix((np.concatenate(kept_values),(np.concatenate(kept_rows),np.concatenate(kept_columns))),shape=(len(genes),len(d)),dtype=np.int32).tocsr()
    del kept_rows,kept_columns,kept_values
    libraries=np.asarray(x.sum(axis=0,dtype=np.int64)).ravel()
    with out.open('xb') as f:
        for i in range(0,len(genes),256):x[i:i+256].toarray().tofile(f)
    del x
    datasets['D1']=(d,genes,libraries)
    # D2: donor IDs and author state definitions, not capture counts as n.
    d2meta=bind(CACHE/'Sountoulidis_lungdev_meta.tsv')
    raw=pd.read_csv(d2meta,sep='\t'); assert len(raw)==163236
    assert raw['Unnamed: 0'].is_unique and raw.groupby('name').donor.nunique().max()==1
    labels={cfg['D2'][s]:s for s in ['start','intermediate','destination']}
    d=raw[raw.indiv_clusters.isin(labels)].copy().rename(columns={'Unnamed: 0':'cell_id','donor':'unit','indiv_annotations':'author_state'})
    d['state']=d.indiv_clusters.map(labels); d=d[['cell_id','unit','state','author_state','age','name','nCount_RNA']]
    d,cc=choose(d,cfg); cc['role']='D2'; coverage.append(cc.reset_index())
    # Validate deposited integers against the independent GEO raw-count file at
    # ten fixed, annotation-matched barcodes, across every shared feature.
    check_path=bind(CACHE/'GSM6645126_counts.h5')
    with h5py.File(check_path,'r') as h:
        m=h['matrix']; rx=csc_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=tuple(m['shape'][:]))
        rg=m['features/name'].asstr()[:].tolist(); rb=m['barcodes'].asstr()[:].tolist()
    ci=raw.loc[raw.name=='10X140_1','Unnamed: 0'].iloc[:10].tolist()
    ri={s.rsplit('_',1)[-1]:i for i,s in enumerate(rb)}
    ri=[ri[s.split(':')[-1].removesuffix('x')] for s in ci]
    rx=rx[:,ri]
    genes,libraries,raw_checks=dense_text(bind(CACHE/'Sountoulidis_counts.tsv.gz'),d,True,WORK/'D2_counts.bin',(rx,rg,ci))
    # The UCSC version may filter genes; publish any library-total discrepancy.
    d['deposited_metadata_total']=d.pop('nCount_RNA'); datasets['D2']=(d,genes,libraries)
    # V1: known mouse mappings only. The branch amendment predates all scores.
    vpath=bind(CACHE/'Haber_counts.txt.gz')
    with gzip.open(vpath,'rt') as f: cells=f.readline().strip().split('\t')
    mapping=pd.read_csv(bind(BASE/'tables/v1_haber_verified_mouse_mapping.csv')).set_index('batch').mouse_id.to_dict()
    labels={cfg['V1'][s]:s for s in ['start','intermediate','destination']}
    rows=[]
    for cell in cells:
        batch,barcode,label=cell.split('_',2)
        if batch in mapping and label in labels:
            rows.append(dict(cell_id=cell,unit=mapping[batch],state=labels[label],author_state=label,batch=batch))
    d,cc=choose(pd.DataFrame(rows),cfg); cc['role']='V1'; coverage.append(cc.reset_index())
    genes,libraries,_=dense_text(vpath,d,False,WORK/'V1_counts.bin');datasets['V1']=(d,genes,libraries)
    for role,(d,genes,libs) in datasets.items():
        assert (libs>0).all() and len(d.unit.unique())>=cfg['unit_floor']
        d['library_total']=libs; d.to_csv(WORK/f'{role}_cells.tsv',sep='\t',index=False)
        (WORK/f'{role}_genes.json').write_text(json.dumps(genes)+'\n')
        write_json(WORK/f'{role}_shape.json',dict(genes=len(genes),cells=len(d),dtype='int32'))
    # One-to-one mapping is conservative, with ambiguous symbols removed.
    mgip=bind(CACHE/'MGI_MouseHuman.rpt'); mgi=pd.read_csv(mgip,sep='\t',dtype=str).drop_duplicates(['DB Class Key','NCBI Taxon ID','Symbol'])
    pairs=[]
    for key,g in mgi.groupby('DB Class Key',sort=False):
        mm=g[g['NCBI Taxon ID']=='10090'].Symbol.unique(); hh=g[g['NCBI Taxon ID']=='9606'].Symbol.unique()
        if len(mm)==len(hh)==1:pairs.append(dict(class_key=key,mouse=mm[0],human=hh[0]))
    pairs=pd.DataFrame(pairs); pairs=pairs[~pairs.mouse.duplicated(False)&~pairs.human.duplicated(False)]
    for role,(_,genes,_) in datasets.items():
        symbol=pd.Series(genes); valid=set(symbol[~symbol.duplicated(False)])
        pairs=pairs[pairs['human' if role=='D2' else 'mouse'].isin(valid)]
    pairs=pairs.sort_values('human'); assert len(pairs)>=10000
    pairs.to_csv(OUT/'ortholog_universe.tsv',sep='\t',index=False)
    pd.concat(coverage,ignore_index=True).to_csv(OUT/'eligibility.tsv',sep='\t',index=False)
    q=[]
    for role,(d,genes,libs) in datasets.items():
        for (unit,state),z in d.groupby(['unit','state']):
            q.append(dict(role=role,unit=unit,state=state,cells=len(z),median_library_total=float(z.library_total.median()),cells_at_500=int((z.library_total>=500).sum())))
    pd.DataFrame(q).to_csv(OUT/'depth_coverage.tsv',sep='\t',index=False)
    prepared={p.name:dict(bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(WORK.iterdir())}
    write_json(OUT/'preparation.json',dict(completed_utc=datetime.now(timezone.utc).isoformat(),script_sha256=sha(Path(__file__)),config_sha256=sha(BASE/'config/pilot_v1.json'),inputs=inputs,prepared=prepared,raw_count_comparisons=raw_checks,orthologs=len(pairs),programme_effects_computed=False,roles={r:dict(units=int(d.unit.nunique()),cells=len(d),genes=len(g)) for r,(d,g,_) in datasets.items()}))
    print(json.dumps({'status':'prepared_no_programme_effects','orthologs':len(pairs),'raw_count_comparisons':raw_checks,'roles':{r:dict(units=int(d.unit.nunique()),cells=len(d)) for r,(d,_,_) in datasets.items()}}))

if __name__=='__main__': main()
