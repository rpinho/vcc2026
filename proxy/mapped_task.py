"""Regenerable file-backed CSR tasks: bound anonymous RAM on the 16-GB Mini."""
import hashlib
import json
import fcntl
from pathlib import Path
import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

P=Path(__file__).resolve().parent


def source_path(dest):
    return P.parent/'rsi/artifacts/tasks/H1_rich_search.h5ad' if dest=='H1' else P/f'tasks/{dest}.h5ad'


def load(dest):
    source=source_path(dest);backed=ad.read_h5ad(source,backed='r')
    obs=backed.obs.copy();var=backed.var.copy();backed.file.close()
    if dest=='H1':
        role=obs.evaluation_role.astype(str).to_numpy()
        obs['target_gene']=np.where(np.isin(role,['input_control','reference_control']),'non-targeting',obs.target_gene.astype(str).to_numpy())
    obs['target_gene']=obs.target_gene.astype(str)
    labels=obs.target_gene.to_numpy();mask=labels=='non-targeting'
    folder=P/f'v2_cache/mapped/{dest}';folder.mkdir(parents=True,exist_ok=True)
    identity=dict(size=source.stat().st_size,mtime_ns=source.stat().st_mtime_ns,source=str(source),format='stream-binomial-csr-v1')
    receipt_path=folder/'receipt.json'
    with (folder/'build.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        receipt=json.loads(receipt_path.read_text()) if receipt_path.exists() else {}
        if receipt.get('source_identity')!=identity:
            build(source,folder,identity,labels,var.index.to_numpy(str),mask)
            receipt=json.loads(receipt_path.read_text())
    arrays={k:np.load(folder/(k+'.npy'),mmap_mode='r') for k in ('data','indices','indptr')}
    n=receipt['nnz'];x=sparse.csr_matrix((arrays['data'][:n],arrays['indices'][:n],arrays['indptr']),shape=receipt['shape'],copy=False)
    a=ad.AnnData(x,obs=obs,var=var)
    panel=np.array(sorted(set(labels)-{'non-targeting'})) if dest=='H1' else pd.read_csv(P/f'tasks/{dest}_panel.csv').target_gene.astype(str).to_numpy()
    panel=np.array([t for t in panel if (labels==t).any()])
    return a,panel,receipt


def build(source,folder,identity,labels,genes,mask):
    print('BUILD_MAPPED_TASK',source,flush=True)
    with h5py.File(source,'r') as h:
        g=h['X'];ptr=g['indptr'][:];depth=[]
        ctrlrows=np.flatnonzero(mask)
        contiguous=ctrlrows[-1]-ctrlrows[0]+1==len(ctrlrows)
        ranges=[(int(ctrlrows[0]),int(ctrlrows[-1])+1)] if contiguous else [(lo,min(lo+512,len(labels))) for lo in range(0,len(labels),512)]
        for leftrow,rightrow in ranges:
            for lo in range(leftrow,rightrow,512):
                hi=min(lo+512,rightrow);left=int(ptr[lo]);right=int(ptr[hi])
                if not mask[lo:hi].any():continue
                values=g['data'][left:right]
                assert np.all(np.diff(ptr[lo:hi+1])>0),'empty source row'
                sums=np.add.reduceat(values,ptr[lo:hi]-left,dtype=np.int64)
                depth.extend(sums[mask[lo:hi]].tolist())
        median=float(np.median(depth));probability=min(1.,20000./median)
        changes=np.diff(np.r_[False,mask,False].astype(np.int8));starts=np.flatnonzero(changes==1);stops=np.flatnonzero(changes==-1)
        blocks=[];offset=0
        for start,stop in zip(starts,stops):
            blocks.append(dict(data_start=int(ptr[start]),data_stop=int(ptr[stop]),control_start=offset,control_stop=offset+int(stop-start)))
            offset+=int(stop-start)
        by_start={b['data_start']:b for b in blocks}
        digest=hashlib.sha256();arrays={}
        for name in ('data','indices','indptr'):
            dataset=g[name]
            arrays[name]=np.lib.format.open_memmap(folder/(name+'.npy'),mode='w+',dtype=dataset.dtype,shape=dataset.shape)
            bounds=sorted(set([0,len(dataset),*range(0,len(dataset),1_000_000),
                              *((b['data_start'] for b in blocks) if name=='data' and probability<1 else ())]))
            rng=np.random.default_rng(20261018)
            for lo,hi in zip(bounds[:-1],bounds[1:]):
                values=dataset[lo:hi];digest.update(memoryview(np.ascontiguousarray(values)).cast('B'))
                if name=='data' and probability<1:
                    if lo in by_start:by_start[lo]['rng_state']=rng.bit_generator.state
                    values=rng.binomial(values,probability).astype(dataset.dtype)
                arrays[name][lo:hi]=values
            arrays[name].flush()
            print('MAPPED_ARRAY',name,len(dataset),flush=True)
        digest.update('\0'.join(labels).encode());digest.update('\0'.join(genes).encode())
        base=digest.hexdigest()
        x=sparse.csr_matrix((arrays['data'],arrays['indices'],arrays['indptr']),shape=g.attrs['shape'],copy=False)
        if probability<1:x.eliminate_zeros()
        for value in arrays.values():value.flush()
        receipt=dict(source_identity=identity,shape=list(map(int,x.shape)),nnz=int(x.nnz),base_matrix_sha256=base,
                     native_control_median_umi=median,thinning_probability=probability,
                     control_rng_blocks=blocks if probability<1 else [])
        temporary=folder/'receipt.tmp';temporary.write_text(json.dumps(receipt,indent=2));temporary.replace(folder/'receipt.json')
        print('MAPPED_TASK_READY',int(x.nnz),flush=True)
