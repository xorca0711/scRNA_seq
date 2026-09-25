"""Library-wise full human QC, frozen-reference annotation and pseudobulks.

No full-cohort cell matrix is materialized. All histologies retain deposited
labels; the atlas never supplies a malignant, precursor, DATP or KAC label.
"""
from pathlib import Path
import os,sys,json,time,gc
from datetime import datetime,timezone
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ.setdefault(k,'2')
PAPER=Path(__file__).resolve().parents[1];ROOT=PAPER.parents[1];sys.path.insert(0,str(ROOT))
from analysis.lib.provenance import write_json_atomic,sha256_file,code_identity
from u5_human_count_store import build_store,read_panel_store,aggregate_store


def transfer(classifier,reference_labels,latent):
    import numpy as np,pandas as pd
    from sklearn import config_context
    output=[]
    with config_context(working_memory=128):
        for start in range(0,len(latent),256):
            dist,idx=classifier.kneighbors(latent[start:start+256]);sd=np.std(dist,axis=1);scale=(2/np.where(sd>0,sd,1e-12))**2;w=np.exp(-dist/scale[:,None]);w/=w.sum(axis=1,keepdims=True)
            block={'mean_neighbor_distance':dist.mean(axis=1)}
            for label in reference_labels:
                neighbors=reference_labels[label].to_numpy()[idx];best=[];unc=[]
                for i,row in enumerate(neighbors):
                    choices,inverse=np.unique(row,return_inverse=True);votes=np.bincount(inverse,weights=w[i]);winner=np.argmax(votes);best.append(choices[winner]);unc.append(max(0,1-votes[winner]))
                block[label]=best;block[label+'_uncertainty']=unc
            output.append(pd.DataFrame(block))
    return pd.concat(output,ignore_index=True)


def main():
    import numpy as np,pandas as pd,anndata as ad,h5py,torch,scvi
    from anndata.io import read_elem
    from sklearn.neighbors import NearestNeighbors
    out=PAPER/'trials/u5_human_full';cache=PAPER/'cache/u5_human_full';cache.mkdir(parents=True,exist_ok=True)
    reference=ROOT/'Thesis/gate1_04_sikkema_2023_hlca/trials/s2_reference_mapping/reference';model_path=PAPER/'cache/u5_human_pilot/query_model'
    order=pd.read_csv(reference/'HLCA_reference_model_gene_order_ids_and_symbols.csv');legacy=pd.read_csv(reference/'HLCA_reference_model/var_names.csv',header=None)[0].tolist();order=order.set_index('gene_id').loc[legacy].reset_index()
    pilot_spec=json.loads((PAPER/'trials/u5_human_pilot/input_specification.json').read_text());panel=pilot_spec['panel']
    record=out/'processing_run_record.json';assert not record.exists(),'Existing run requires explicit checkpoint-aware continuation'
    assert json.loads((out/'count_store_validation.json').read_text())['status']=='passed'
    spec=dict(frozen_utc=datetime.now(timezone.utc).isoformat(),scope='all 75 source-deposited libraries, all 23 patients; library counts are not patient replication',qc=pilot_spec['qc'],annotation='frozen three-patient pilot scArches adaptor; exact HLCA-core weighted 50-neighbor transfer; all fine labels with uncertainty <=0.2 primary, <=0.3 sensitivity; original histology retained',annotation_validation='pilot multi-marker review supports major compartments; ambiguous and low-confidence cells excluded from primary subtype pseudobulks; distant/new neoplastic states are not forced into KAC/malignancy classes',background_RNA_limit='surfactant RNA is widespread outside epithelial assignments in the pilot; identity uses multigene/reference support, not one marker; LR remains RNA compatibility',model_sha256=sha256_file(model_path/'model.pt'),reference_gene_map_sha256=sha256_file(reference/'HLCA_reference_model_gene_order_ids_and_symbols.csv'),count_store_code_sha256=sha256_file(PAPER/'trials/u5_human_count_store.py'),normalization='raw counts preserved; TMM/pathway normalization follows donor-by-histology aggregation',unit='patient; repeated same-histology libraries pooled after annotation',primary_uncertainty=.2,sensitivity_uncertainty=.3)
    write_json_atomic(out/'processing_specification.json',spec)
    state=dict(status='running',pid=os.getpid(),started_utc=spec['frozen_utc'],stage='load_reference',code=code_identity(ROOT,__file__),completed_files=[]);write_json_atomic(record,state);start=time.monotonic();model=None;all_units=[];quality=[]
    try:
        with h5py.File(reference/'HLCA_full_v1.1_emb.h5ad','r') as h:
            obs=read_elem(h['obs']);core=obs.core_or_extension.astype(str).str.lower().eq('core').to_numpy();keys=[k for k in h['obsm'] if 'scanvi' in k.lower()] or [k for k in h['obsm'] if 'emb' in k.lower()]
            embedding=h['obsm'][keys[0]][:][core] if len(keys)==1 else h['X'][:][core]
        labels=['ann_level_1','ann_level_2','ann_level_3','ann_level_4','ann_level_5','ann_finest_level'];ref_labels=obs.loc[core,labels].astype(str).reset_index(drop=True);del obs
        classifier=NearestNeighbors(n_neighbors=50,metric='euclidean',n_jobs=2).fit(embedding);torch.set_num_threads(2);scvi.settings.seed=20260924
        planned=json.loads((out/'acquisition_specification.json').read_text())['jobs']
        for job in planned:
            gsm=job['gsm'];state.update(current_gsm=gsm,stage='await_processed_download');write_json_atomic(record,state)
            while True:
                downloaded=json.loads((out/'acquisition_run_record.json').read_text());files={f['gsm']:f for f in downloaded['files']}
                if gsm in files:break
                if downloaded['status']=='failed':raise RuntimeError('Acquisition failed before '+gsm)
                time.sleep(5)
            source=ROOT/files[gsm]['path'];here=cache/gsm;here.mkdir(exist_ok=True);store=here/'raw_counts.csc.h5';t0=time.monotonic()
            state['stage']='stream_full_counts';write_json_atomic(record,state)
            matrix,meta=read_panel_store(store,panel) if store.exists() else build_store(source,store,panel)
            n=len(meta['barcodes']);m=pd.DataFrame({k:meta[k] for k in ['full_library_size','n_genes_by_counts','pct_counts_mt','full_gene_max_count']})
            m.index=[gsm+':'+x for x in meta['barcodes']];m.index.name='cell_id';m['gsm']=gsm;m['patient']=str(job['patient']);m['histology']=job['histology'];m['source_cell_index']=np.arange(n)
            m['qc_pass']=(m.n_genes_by_counts>=500)&(m.full_library_size>=1000)&(m.pct_counts_mt<20)
            qc=np.flatnonzero(m.qc_pass);a=ad.AnnData(matrix[qc].copy(),obs=m.iloc[qc].copy(),var=pd.DataFrame(index=panel));del matrix
            coverage=order.gene_symbol.isin(meta['genes']).mean();assert coverage>=.7
            query=ad.AnnData(a[:,order.gene_symbol].X.astype(np.float32).copy(),obs=a.obs.copy(),var=pd.DataFrame(index=legacy));query.obs['dataset']='GSE308103_FixedRNA';query.obs['scanvi_label']='unlabeled'
            state.update(stage='frozen_reference_annotation',current_qc_cells=len(query));write_json_atomic(record,state)
            if model is None:model=scvi.model.SCANVI.load(str(model_path),adata=query,accelerator='cpu')
            latent=model.get_latent_representation(adata=query,batch_size=256);tr=transfer(classifier,ref_labels,latent);tr.index=a.obs_names;a.obs=a.obs.join(tr,validate='one_to_one')
            a.write_h5ad(here/'annotated_panel.h5ad',compression='gzip');a.obs.to_csv(here/'cell_annotations.csv.gz',compression='gzip',index_label='cell_id')
            np.savez_compressed(here/'query_latent.npz',latent=latent)
            state['stage']='full_gene_pseudobulk';write_json_atomic(record,state)
            for cutoff in [.2,.3]:
                confident=a.obs.ann_finest_level_uncertainty<=cutoff;levels=sorted(a.obs.loc[confident,'ann_finest_level'].unique());code_map={x:i for i,x in enumerate(levels)};codes=np.full(n,-1,dtype=np.int32)
                accepted=a.obs.loc[confident];codes[accepted.source_cell_index]=accepted.ann_finest_level.map(code_map).to_numpy()
                genes,counts,detected=aggregate_store(store,codes,len(levels));tag=f'unc{int(cutoff*100):02d}'
                units=[]
                for i,label in enumerate(levels):
                    group=accepted[accepted.ann_finest_level==label]
                    units.append(dict(unit_id=gsm+'_'+tag+'_'+str(i),gsm=gsm,patient=job['patient'],histology=job['histology'],uncertainty_cutoff=cutoff,label=label,compartment=group.ann_level_1.mode().iloc[0],lineage=group.ann_level_2.mode().iloc[0],cells=len(group),full_library_sum=int(group.full_library_size.sum()),panel_gene_coverage=float(coverage)))
                pd.DataFrame(units).to_csv(here/(tag+'_units.csv'),index=False);np.savez_compressed(here/(tag+'_pseudobulks.npz'),genes=np.asarray(genes,dtype=str),counts=counts,detected_cells=detected)
                all_units.extend(units)
            q=dict(gsm=gsm,patient=job['patient'],histology=job['histology'],input_cells=n,qc_cells=len(a),confident_primary=int((a.obs.ann_finest_level_uncertainty<=.2).sum()),confident_sensitivity=int((a.obs.ann_finest_level_uncertainty<=.3).sum()),source_counts=int(meta['full_library_size'].sum()),model_gene_fraction=float(coverage),median_neighbor_distance=float(a.obs.mean_neighbor_distance.median()),elapsed_seconds=round(time.monotonic()-t0,1),source_sha256=files[gsm]['sha256'])
            quality.append(q);state['completed_files'].append(gsm);state.update(completed_libraries=len(quality),qc_cells_total=sum(r['qc_cells'] for r in quality),elapsed_seconds=round(time.monotonic()-start,1));write_json_atomic(record,state)
            pd.DataFrame(quality).to_csv(out/'library_qc_and_annotation.csv',index=False);pd.DataFrame(all_units).to_csv(out/'library_subtype_units.csv',index=False)
            print(gsm,len(quality),len(planned),'QC',len(a),'seconds',q['elapsed_seconds'],flush=True)
            del a,query,meta,m,tr,latent,counts,detected;gc.collect()
        state.update(status='completed_inputs_ready',stage='all_library_pseudobulks_ready')
    except Exception as exc:state.update(status='failed',error=type(exc).__name__+': '+str(exc));raise
    finally:state.update(updated_utc=datetime.now(timezone.utc).isoformat(),elapsed_seconds=round(time.monotonic()-start,1));write_json_atomic(record,state)

if __name__=='__main__':main()
