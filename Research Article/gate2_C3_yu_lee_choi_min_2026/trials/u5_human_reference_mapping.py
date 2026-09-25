"""Bounded scArches annotation against the existing HLCA reference.

Training is balanced across patient/histology strata, independent of gene
effects. All QC nuclei are then encoded. Atlas labels are annotation candidates,
not malignancy, precursor or KAC diagnoses; marker review remains required.
"""
from pathlib import Path
import os,sys,json,time
from datetime import datetime,timezone
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ.setdefault(k,'2')
PAPER=Path(__file__).resolve().parents[1];ROOT=PAPER.parents[1];sys.path.insert(0,str(ROOT))
from analysis.lib.provenance import write_json_atomic,code_identity,sha256_file

def main():
    import numpy as np,pandas as pd,anndata as ad,torch,scvi,h5py
    from scipy import sparse
    from anndata.io import read_elem
    from sklearn.neighbors import NearestNeighbors
    from sklearn import config_context
    out=PAPER/'trials/u5_human_pilot';cache=PAPER/'cache/u5_human_pilot'
    ref=ROOT/'Research Article/gate1_04_sikkema_2023_hlca/trials/s2_reference_mapping/reference';model_path=ref/'HLCA_reference_model_converted'
    src=cache/'annotation_input.h5ad'
    if (out/'mapping_run_record.json').exists():raise RuntimeError('Mapping run exists; inspect before restarting')
    spec=dict(frozen_utc=datetime.now(timezone.utc).isoformat(),seed=20260924,training_cap_per_patient_histology=2000,max_epochs=500,early_stopping_patience=10,early_stopping_min_delta=.001,train_size=.9,batch_size=128,knn_k=50,uncertainty_primary=.2,uncertainty_sensitivity=.3,scope='HLCA candidate annotation for broad cell compartments, no KAC/malignancy inference',query_format='Fixed-RNA FFPE nuclei; 1731/2000 model genes measured, missing model genes zero-filled',training='uniform seeded nuclei per patient/histology; repeated libraries sampled proportionally inside their biological stratum; encode all QC nuclei after adaptation',input_sha256=sha256_file(src),model_sha256=sha256_file(model_path/'model.pt'),gene_map_sha256=sha256_file(ref/'HLCA_reference_model_gene_order_ids_and_symbols.csv'))
    write_json_atomic(out/'mapping_specification.json',spec)
    state=dict(status='running',stage='prepare',pid=os.getpid(),started_utc=spec['frozen_utc'],code=code_identity(ROOT,__file__));write_json_atomic(out/'mapping_run_record.json',state);start=time.monotonic()
    try:
        a=ad.read_h5ad(src);order=pd.read_csv(ref/'HLCA_reference_model_gene_order_ids_and_symbols.csv')
        legacy=pd.read_csv(ref/'HLCA_reference_model/var_names.csv',header=None)[0].astype(str).tolist()
        order=order.set_index('gene_id').loc[legacy].reset_index()
        assert order.gene_symbol.is_unique
        x=a[:,order.gene_symbol].X.tocsr().astype(np.float32)
        query=ad.AnnData(x,obs=a.obs.copy(),var=pd.DataFrame(index=legacy));del a,x
        query.obs['dataset']='GSE308103_FixedRNA';query.obs['scanvi_label']='unlabeled'
        scvi.settings.seed=spec['seed'];torch.set_num_threads(2)
        selected=[];rng=np.random.default_rng(spec['seed'])
        for _,idx in query.obs.groupby(['patient','histology'],observed=True).indices.items():selected.extend(rng.choice(idx,min(len(idx),2000),replace=False))
        train=query[np.sort(selected)].copy()
        train.obs[['gsm','patient','histology']].to_csv(cache/'mapping_training_cells.csv.gz',index_label='cell_id',compression='gzip')
        train.obs.groupby(['patient','histology','gsm'],observed=True).size().rename('training_nuclei').to_csv(out/'mapping_training_coverage.csv')
        scvi.model.SCANVI.prepare_query_anndata(train,str(model_path))
        model=scvi.model.SCANVI.load_query_data(train,str(model_path),freeze_dropout=True)
        state.update(stage='adapt_reference',training_nuclei=train.n_obs,all_nuclei=query.n_obs);write_json_atomic(out/'mapping_run_record.json',state)
        model.train(max_epochs=spec['max_epochs'],train_size=.9,batch_size=128,plan_kwargs=dict(weight_decay=0.,reduce_lr_on_plateau=True,lr_patience=8,lr_factor=.1),check_val_every_n_epoch=1,early_stopping=True,early_stopping_monitor='elbo_validation',early_stopping_patience=10,early_stopping_min_delta=.001,early_stopping_mode='min',n_samples_per_label=None,accelerator='cpu',enable_progress_bar=False)
        history={k:[float(x) for x in np.asarray(v).ravel()] for k,v in model.history.items() if 'elbo' in k}
        write_json_atomic(out/'mapping_training_history.json',history);model.save(str(cache/'query_model'),overwrite=False)
        scvi.model.SCANVI.prepare_query_anndata(query,str(model_path))
        latent=model.get_latent_representation(adata=query,batch_size=256)
        np.savez_compressed(cache/'query_latent.npz',latent=latent,cell_id=np.asarray(query.obs_names,dtype=str))
        state.update(stage='reference_label_transfer',epochs=len(history['elbo_train']));write_json_atomic(out/'mapping_run_record.json',state)
        labels=['ann_level_1','ann_level_2','ann_level_3','ann_level_4','ann_level_5','ann_finest_level']
        with h5py.File(ref/'HLCA_full_v1.1_emb.h5ad','r') as h:
            obs=read_elem(h['obs']);core=obs.core_or_extension.astype(str).str.lower().eq('core').to_numpy()
            possible=[k for k in h['obsm'] if 'scanvi' in k.lower()] or [k for k in h['obsm'] if 'emb' in k.lower()]
            emb=h['obsm'][possible[0]][:][core] if len(possible)==1 else h['X'][:][core]
            assert emb.shape[1]==latent.shape[1]
        rl=obs.loc[core,labels].astype(str).reset_index(drop=True);del obs
        classifier=NearestNeighbors(n_neighbors=50,metric='euclidean',n_jobs=2).fit(emb)
        output=[]
        with config_context(working_memory=128):
            for start_row in range(0,len(latent),256):
                distance,idx=classifier.kneighbors(latent[start_row:start_row+256])
                sd=np.std(distance,axis=1);scale=(2/np.where(sd>0,sd,1e-12))**2
                weights=np.exp(-distance/scale[:,None]);weights/=weights.sum(axis=1,keepdims=True)
                block={'mean_neighbor_distance':distance.mean(axis=1)}
                for label in labels:
                    neighbors=rl[label].to_numpy()[idx];best=[];unc=[]
                    for i,row in enumerate(neighbors):
                        choices,inverse=np.unique(row,return_inverse=True);votes=np.bincount(inverse,weights=weights[i]);winner=np.argmax(votes)
                        best.append(choices[winner]);unc.append(max(0,1-votes[winner]))
                    block[label]=best;block[label+'_uncertainty']=unc
                output.append(pd.DataFrame(block,index=query.obs_names[start_row:start_row+256]))
        result=pd.concat(output);result=query.obs[['gsm','patient','histology']].join(result,validate='one_to_one')
        result.to_csv(cache/'reference_candidate_annotations.csv.gz',index_label='cell_id',compression='gzip')
        result.groupby(['patient','histology','ann_level_1','ann_level_2','ann_finest_level'],observed=True).size().rename('nuclei').to_csv(out/'candidate_annotation_counts.csv')
        state.update(status='completed_candidate_annotations_require_marker_review',reference_core_cells=len(emb),all_nuclei=len(result),median_finest_uncertainty=float(result.ann_finest_level_uncertainty.median()),finest_uncertain_fraction=float((result.ann_finest_level_uncertainty>.2).mean()))
    except Exception as exc:state.update(status='failed',error=type(exc).__name__+': '+str(exc));raise
    finally:state.update(elapsed_seconds=round(time.monotonic()-start,1),updated_utc=datetime.now(timezone.utc).isoformat());write_json_atomic(out/'mapping_run_record.json',state)

if __name__=='__main__':main()
