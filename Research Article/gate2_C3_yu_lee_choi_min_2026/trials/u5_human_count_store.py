"""Bounded gene-row count parser and column-block pseudobulk utilities.

The cache is an explicitly documented CSC store, not an AnnData file. It
retains full raw counts so annotation and downstream aggregation do not require
re-parsing the dense deposited text or loading a full library into RAM.
"""
from pathlib import Path
import gzip,json,re

def build_store(source,target,selected_genes):
    import numpy as np,h5py
    from scipy import sparse
    target=Path(target)
    if target.exists():raise RuntimeError('Count store already exists')
    temp=target.with_suffix('.partial.h5');gm={g:i for i,g in enumerate(selected_genes)}
    pr=[];pc=[];pv=[];genes=[];ptr=[0];rr=[];vv=[];written=0
    with gzip.open(source,'rt') as text,h5py.File(temp,'w') as h:
        barcodes=text.readline().rstrip('\r\n').split('\t');n=len(barcodes);assert len(set(barcodes))==n
        total=np.zeros(n,dtype=np.int64);ng=np.zeros(n,dtype=np.int32);mt=np.zeros(n,dtype=np.int64);maximum=np.zeros(n,dtype=np.int32)
        data=h.create_dataset('data',shape=(0,),maxshape=(None,),dtype='i4',chunks=(262144,),compression='gzip',compression_opts=1)
        indices=h.create_dataset('indices',shape=(0,),maxshape=(None,),dtype='i4',chunks=(262144,),compression='gzip',compression_opts=1)
        def flush():
            nonlocal written
            if not rr:return
            ix=np.concatenate(rr);value=np.concatenate(vv);size=len(ix);data.resize((written+size,));indices.resize((written+size,));data[written:written+size]=value;indices[written:written+size]=ix;written+=size;rr.clear();vv.clear()
        for line in text:
            gene,values=line.rstrip('\r\n').split('\t',1);genes.append(gene)
            if re.search(r'[^0-9\t]',values):raise ValueError('Counts must be lexical nonnegative integers')
            v=np.fromstring(values,sep='\t',dtype=np.int64)
            assert len(v)==n and np.all(v>=0) and v.max()<=np.iinfo(np.int32).max
            total+=v;ng+=v>0
            np.maximum(maximum,v,out=maximum)
            if gene.startswith('MT-'):mt+=v
            idx=np.flatnonzero(v).astype(np.int32);val=v[idx].astype(np.int32);rr.append(idx);vv.append(val);ptr.append(ptr[-1]+len(idx))
            if gene in gm:pr.append(idx);pc.append(np.full(len(idx),gm[gene],dtype=np.int32));pv.append(val)
            if len(genes)%256==0:flush()
        flush();assert len(set(genes))==len(genes) and written==ptr[-1]
        h.create_dataset('indptr',data=np.asarray(ptr,dtype=np.int64));h.create_dataset('genes',data=np.asarray(genes,dtype=object),dtype=h5py.string_dtype())
        h.create_dataset('barcodes',data=np.asarray(barcodes,dtype=object),dtype=h5py.string_dtype());h.create_dataset('full_library_size',data=total);h.create_dataset('n_genes_by_counts',data=ng);h.create_dataset('mitochondrial_counts',data=mt);h.create_dataset('full_gene_max_count',data=maximum)
        h.attrs['format']='review_gene_row_CSC_v1';h.attrs['shape']=(n,len(genes));h.attrs['total_counts']=int(total.sum())
    temp.replace(target)
    panel=sparse.coo_matrix((np.concatenate(pv),(np.concatenate(pr),np.concatenate(pc))),shape=(n,len(selected_genes))).tocsr()
    return panel,dict(barcodes=barcodes,genes=genes,full_library_size=total,n_genes_by_counts=ng,pct_counts_mt=100*mt/np.maximum(total,1),full_gene_max_count=maximum)


def read_panel_store(path,selected_genes):
    import h5py,numpy as np
    from scipy import sparse
    rows=[];cols=[];values=[];gm={g:i for i,g in enumerate(selected_genes)}
    with h5py.File(path,'r') as h:
        genes=h['genes'].asstr()[:];ptr=h['indptr'][:];n,ng=h.attrs['shape']
        for begin in range(0,ng,256):
            end=min(ng,begin+256);start=int(ptr[begin]);stop=int(ptr[end]);indices=h['indices'][start:stop];data=h['data'][start:stop]
            for j in range(begin,end):
                if genes[j] not in gm:continue
                lo=int(ptr[j]-start);hi=int(ptr[j+1]-start);rows.append(indices[lo:hi]);cols.append(np.full(hi-lo,gm[genes[j]],dtype=np.int32));values.append(data[lo:hi])
        meta=dict(barcodes=h['barcodes'].asstr()[:].tolist(),genes=genes.tolist(),full_library_size=h['full_library_size'][:],n_genes_by_counts=h['n_genes_by_counts'][:],full_gene_max_count=h['full_gene_max_count'][:])
        meta['pct_counts_mt']=100*h['mitochondrial_counts'][:]/np.maximum(meta['full_library_size'],1)
    matrix=sparse.coo_matrix((np.concatenate(values),(np.concatenate(rows),np.concatenate(cols))),shape=(n,len(selected_genes))).tocsr()
    return matrix,meta


def aggregate_store(path,group_codes,n_groups,block_genes=256):
    import h5py,numpy as np
    from scipy import sparse
    with h5py.File(path,'r') as h:
        assert h.attrs['format']=='review_gene_row_CSC_v1'
        n,ng=h.attrs['shape'];assert len(group_codes)==n
        genes=h['genes'].asstr()[:];ptr=h['indptr'][:];keep=np.flatnonzero(group_codes>=0)
        group=sparse.csr_matrix((np.ones(len(keep),dtype=np.int64),(group_codes[keep],keep)),shape=(n_groups,n))
        sums=np.zeros((n_groups,ng),dtype=np.int64);detected=np.zeros((n_groups,ng),dtype=np.int64)
        for begin in range(0,ng,block_genes):
            end=min(ng,begin+block_genes);start=int(ptr[begin]);stop=int(ptr[end])
            block=sparse.csc_matrix((h['data'][start:stop],h['indices'][start:stop],ptr[begin:end+1]-start),shape=(n,end-begin))
            sums[:,begin:end]=(group@block).toarray();detected[:,begin:end]=(group@(block>0)).toarray()
        expected=int(h['full_library_size'][:][keep].sum());assert int(sums.sum())==expected
        return genes,sums,detected
