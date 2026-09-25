"""Whole-section descriptive measurements, with coordinate/ROI audit.

Source tissue masks are not pathology-region annotations. No neighborhood
significance test or cell-cell communication claim is generated here.
"""
from pathlib import Path
import json,sys,os,time,tarfile,gzip,shutil,gc
from datetime import datetime,timezone
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ.setdefault(k,'2')
PAPER=Path(__file__).resolve().parents[1];ROOT=PAPER.parents[1];sys.path.insert(0,str(ROOT))
from analysis.lib.provenance import write_json_atomic,sha256_file,code_identity

def read_visium(path,cache):
 import numpy as np,pandas as pd
 from scipy.io import mmread
 from threadpoolctl import threadpool_limits
 names=['matrix.mtx.gz','features.tsv.gz','barcodes.tsv.gz','tissue_positions.csv','tissue_positions_list.csv','scalefactors_json.json']
 cache.mkdir(parents=True,exist_ok=True);inventory=[];found={}
 with tarfile.open(path) as tar:
  for member in tar.getmembers():
   if not member.isfile() or '/._' in member.name or member.name.startswith('._'):continue
   inventory.append(dict(name=member.name,bytes=member.size));base=Path(member.name).name
   if base not in names:continue
   assert base not in found and member.size<1_000_000_000
   target=(cache/base).resolve();assert target.is_relative_to(cache.resolve())
   with tar.extractfile(member) as source,target.open('wb') as dest:shutil.copyfileobj(source,dest,1024*1024)
   found[base]=target
 for name in ['matrix.mtx.gz','features.tsv.gz','barcodes.tsv.gz']:assert name in found
 features=pd.read_csv(found['features.tsv.gz'],sep='\t',header=None);barcodes=pd.read_csv(found['barcodes.tsv.gz'],sep='\t',header=None)[0].astype(str)
 with threadpool_limits(limits=2),gzip.open(found['matrix.mtx.gz'],'rb') as f:x=mmread(f).T.tocsr()
 assert x.shape==(len(barcodes),len(features)) and barcodes.is_unique
 cols=['barcode','in_tissue','array_row','array_col','pxl_row_in_fullres','pxl_col_in_fullres']
 pos=None;disagreement=None
 if 'tissue_positions.csv' in found:
  pos=pd.read_csv(found['tissue_positions.csv']);assert set(cols)<=set(pos.columns);pos=pos.set_index('barcode')
 if 'tissue_positions_list.csv' in found:
  older=pd.read_csv(found['tissue_positions_list.csv'],header=None,names=cols).set_index('barcode')
  if pos is None:pos=older
  else:
   common=pos.index.intersection(older.index);disagreement=int((pos.loc[common,cols[1:]].to_numpy()!=older.loc[common,cols[1:]].to_numpy()).any(axis=1).sum())
 assert pos is not None and pos.index.is_unique and set(barcodes)<=set(pos.index)
 pos=pos.loc[barcodes].copy();pos.index.name='barcode'
 return x,features[1].astype(str).to_numpy(),pos,inventory,disagreement

def read_h5(path):
 import numpy as np,pandas as pd,h5py
 from scipy import sparse
 with h5py.File(path,'r') as h:
  assert 'matrix' in h;g=h['matrix'];x=sparse.csc_matrix((g['data'][:],g['indices'][:],g['indptr'][:]),shape=tuple(g['shape'][:])).T.tocsr();genes=g['features']['name'][:].astype(str);barcodes=g['barcodes'][:].astype(str);keys=[];h.visit(keys.append)
 return x,genes,pd.DataFrame(index=pd.Index(barcodes,name='barcode')),keys

def main():
 import numpy as np,pandas as pd
 from scipy import sparse
 out=PAPER/'trials/u5_spatial_context';cache=PAPER/'cache/u5_spatial_measurements';cache.mkdir(exist_ok=True)
 record=out/'processing_run_record.json';assert not record.exists();spec=json.loads((out/'acquisition_specification.json').read_text())
 humans=pd.read_csv(PAPER/'trials/u5_ipf_spec/pathway_genes.tsv',sep='\t');mice=pd.read_csv(PAPER/'trials/u4_resources/pathway_genes.tsv',sep='\t')
 markers=['IL1B','IL1A','IL1R1','IL1RAP','IL1RN','IL1R2','SIGIRR','TGFB1','CCL2','CXCL12','AREG','HBEGF','FGF7','FGF10','COL1A1','C1QA','KRT8','KRT5','KRT17','CLDN4','SFTPC','AGER','AQP5','CD3D','CD8A','TNF','IFNG']
 orth=pd.read_csv(PAPER/'trials/u4_resources/strict_one_to_one_orthologs.csv');print('ortholog columns',list(orth.columns),flush=True)
 # Symbol names are read from the frozen table, not inferred by capitalization.
 hm=dict(zip(orth['human_symbol'],orth['mouse_symbol']))
 profiles=[];quality=[];inventories=[];coverage=[];state=dict(status='running',pid=os.getpid(),started_utc=datetime.now(timezone.utc).isoformat(),code=code_identity(ROOT,__file__),completed_files=[]);write_json_atomic(record,state);t=time.monotonic()
 try:
  for job in spec['jobs']:
   while True:
    acq=json.loads((out/'acquisition_run_record.json').read_text());available={r['gsm']:r for r in acq['files']}
    if job['gsm'] in available:break
    if acq['status']=='failed':raise RuntimeError('Spatial acquisition failed before '+job['gsm'])
    time.sleep(5)
   item=available[job['gsm']];path=ROOT/item['path'];gsm=job['gsm'];here=cache/gsm;here.mkdir(exist_ok=True)
   if job['series']=='GSE307534':
    x,genes,obs,inventory,disagreement=read_visium(path,here);coords=True
    inventories.extend(dict(gsm=gsm,**r) for r in inventory)
   else:
    x,genes,obs,keys=read_h5(path);coords=False;disagreement=None;inventories.extend(dict(gsm=gsm,name=k,bytes=None) for k in keys)
   assert x.data.min(initial=0)>=0 and np.array_equal(x.data,np.round(x.data));x.eliminate_zeros()
   # Collapse duplicate symbols before gene-level means; detection uses symbols.
   unique=sorted(set(genes));mapping={g:i for i,g in enumerate(unique)}
   if len(unique)!=len(genes):
    project=sparse.csr_matrix((np.ones(len(genes),dtype=np.int32),(np.arange(len(genes)),[mapping[g] for g in genes])),shape=(len(genes),len(unique)));x=(x@project).tocsr();genes=np.asarray(unique)
   n=x.shape[0];total=np.asarray(x.sum(axis=1)).ravel();ng=np.diff(x.indptr);mt=np.array([g.upper().startswith('MT-') for g in genes]);mtfrac=np.asarray(x[:,mt].sum(axis=1)).ravel()/np.maximum(total,1)*100
   tissue=obs.in_tissue.eq(1).to_numpy() if coords else np.ones(n,dtype=bool);keep=tissue&(total>=500)&(ng>=200)&(mtfrac<=15)
   obs['total_counts']=total;obs['detected_genes']=ng;obs['pct_counts_mt']=mtfrac;obs['qc_pass']=keep;obs.to_csv(here/'spot_qc_and_coordinates.csv.gz',compression='gzip')
   sets=(mice if job['series']=='GSE267228' else humans).groupby('set').gene.apply(list).to_dict();marker_map={g:(hm.get(g) if job['series']=='GSE267228' else g) for g in markers};sets.update({'GENE_'+g:[v] for g,v in marker_map.items() if v})
   present=set(genes);panel=sorted(set(g for gs in sets.values() for g in gs)&present);positions={g:i for i,g in enumerate(genes)};y=x[keep][:,[positions[g] for g in panel]].astype(np.float64);y=y.multiply((10000/np.maximum(total[keep],1))[:,None]).tocsr();linear=np.asarray(y.mean(axis=0)).ravel() if y.shape[0] else np.full(len(panel),np.nan)
   np.log1p(y.data,out=y.data);means=np.asarray(y.mean(axis=0)).ravel() if y.shape[0] else np.full(len(panel),np.nan);lookup={g:i for i,g in enumerate(panel)};maps={}
   for name,gs in sets.items():
    measured=sorted(set(gs)&present);fraction=len(measured)/len(set(gs));eligible=len(measured)>=1 if name.startswith('GENE_') else (fraction>=.7 and len(measured)>=10)
    coverage.append(dict(gsm=gsm,module=name,assayed_genes=len(measured),source_genes=len(set(gs)),assayed_fraction=fraction,eligible=eligible))
    if not measured:continue
    ix=[lookup[g] for g in measured];value=float(np.mean(means[ix]));profiles.append(dict(series=job['series'],gsm=gsm,patient=job.get('patient',gsm),histology=job.get('histology',job['title']),module=name,spots=int(keep.sum()),assayed_fraction=fraction,eligible=bool(eligible and keep.sum()>=100),mean_log1p_norm10000=value,mean_norm10000=float(np.mean(linear[ix])),coordinates_available=coords))
    if name in ['GENE_IL1B','GENE_COL1A1','GENE_C1QA','GENE_KRT8','GENE_CLDN4','HALLMARK_INFLAMMATORY_RESPONSE','HALLMARK_TNFA_SIGNALING_VIA_NFKB']:maps[name]=np.asarray(y[:,ix].mean(axis=1)).ravel()
   if coords:
    view=obs.loc[keep].copy()
    for key,values in maps.items():view[key]=values
    view.to_csv(here/'spatial_display_values.csv.gz',compression='gzip')
   quality.append(dict(series=job['series'],gsm=gsm,patient=job.get('patient',gsm),histology=job.get('histology',job['title']),input_spots=n,qc_spots=int(keep.sum()),genes=len(genes),coordinates_available=coords,position_version_disagreements=disagreement,source_sha256=item['sha256'],independent_pathology_labels=False))
   pd.DataFrame(quality).to_csv(out/'section_qc.csv',index=False);pd.DataFrame(profiles).to_csv(out/'whole_section_programs.csv',index=False);pd.DataFrame(coverage).to_csv(out/'program_coverage.csv',index=False);pd.DataFrame(inventories).to_csv(out/'deposited_member_inventory.csv',index=False)
   state['completed_files'].append(gsm);state.update(elapsed_seconds=round(time.monotonic()-t,1));write_json_atomic(record,state);print(gsm,'spatial/count summary',len(state['completed_files']),'QC',int(keep.sum()),flush=True)
   del x,y,obs;gc.collect()
  state.update(status='completed_measurements_require_report',sections=len(quality),program_rows=len(profiles))
 except Exception as e:state.update(status='failed',error=repr(e));raise
 finally:state.update(elapsed_seconds=round(time.monotonic()-t,1),updated_utc=datetime.now(timezone.utc).isoformat());write_json_atomic(record,state)

if __name__=='__main__':main()
