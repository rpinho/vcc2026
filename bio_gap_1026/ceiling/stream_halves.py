"""One expression pass, whitelist-first CSR slices, five deterministic ID splits.
Lossless byte-shuffled zlib checkpoints; only A is stored, B=full-A.
The retained full sums are used only after a streamed-full equality check.
"""
import runtime as rt
rt.enforce(network=True)
import numpy as np, pandas as pd, scipy.sparse as sp, h5py
from anndata.io import read_elem
from remote import Remote
import json, hashlib, time, os, zlib, struct, gc, subprocess, sys
SEEDS=np.arange(20261027,20261032)
START=time.monotonic(); LIMIT=5.8*3600
SRC=rt.REPO/'model_sidechain/approx/sources/H1'
meta=json.loads((SRC/'meta.json').read_text())
with np.load(rt.ROOT.parent/'results/baseline.npz') as z: targets=z['targets'].astype(str)
labels=list(targets)+['non-targeting']; code={x:i for i,x in enumerate(labels)}
gi=np.arange(len(meta['genes'])); ti=np.array([meta['labels'].index(t) for t in labels])
full={k:np.load(SRC/(k+'.npy'),mmap_mode='r')[ti].copy() for k in ['count_sum','cpm_sum','cpm_sq_sum','n_cells']}
assert len(targets)==168
G=len(gi); L=len(labels)
cp=rt.ROOT/'checkpoints';cp.mkdir(exist_ok=True)
# All live sums are float64. Checkpoints preserve every bit.
a=np.zeros((3,5,L,G),np.float64); total=np.zeros((3,L,G),np.float64)
n=np.zeros((5,L),np.int64); nt=np.zeros(L,np.int64)
statepath=cp/'state.json'
state=json.loads(statepath.read_text()) if statepath.exists() else {'file':0,'row':0,'receipts':[]}

def pack(path,arr):
    raw=arr.view(np.uint8).reshape(-1,8).T.copy().tobytes()
    payload=zlib.compress(raw,6)
    tmp=path.with_suffix('.tmp')
    with tmp.open('wb') as f:f.write(payload)
    tmp.replace(path)

def unpack(path,shape):
    raw=zlib.decompress(path.read_bytes())
    return np.frombuffer(raw,np.uint8).reshape(8,-1).T.copy().view(np.float64).reshape(shape)

# Sharded commits with per-label cursors make interrupted checkpoint writes safe.
cursors=np.zeros((L,2),np.int64)
finalized={p.stem for p in (rt.ROOT/'finalized').glob('*.json')}
for j,t in enumerate(labels):
    if t in finalized:
        nt[j]=full['n_cells'][j];cursors[j]=3,0
        continue
    path=cp/(t+'.bin')
    if path.exists():
        with path.open('rb') as f:
            hsize=struct.unpack('<I',f.read(4))[0]; hdr=json.loads(f.read(hsize)); payload=f.read()
        raw=zlib.decompress(payload)
        vals=np.frombuffer(raw,np.uint8).reshape(8,-1).T.copy().view(np.float64).reshape(18,G)
        a[:,:,j]=vals[:15].reshape(3,5,G);total[:,j]=vals[15:]
        n[:,j]=hdr['n'];nt[j]=hdr['nt'];cursors[j]=hdr['cursor']

def checkpoint(fi,row):
    rt.guard()
    for j,t in enumerate(labels):
        if t in finalized:continue
        if tuple(cursors[j]) >= (fi,row):continue
        vals=np.concatenate([a[:,:,j].reshape(15,G),total[:,j]],0)
        # Counts are integer doubles; byte shuffle compresses their zero mantissas.
        raw=vals.view(np.uint8).reshape(-1,8).T.copy().tobytes()
        payload=zlib.compress(raw,6)
        hdr=json.dumps({'n':n[:,j].tolist(),'nt':int(nt[j]),'cursor':[fi,row]}).encode()
        path=cp/(t+'.bin');tmp=cp/(t+'.tmp')
        size=sum(p.stat().st_size for p in rt.ROOT.rglob('*') if p.is_file())
        if size+len(payload)+len(hdr)+4>195_000_000:
            raise RuntimeError('Lossless checkpoint would exceed disk budget')
        with tmp.open('wb') as f:f.write(struct.pack('<I',len(hdr)));f.write(hdr);f.write(payload)
        tmp.replace(path);cursors[j]=fi,row
    state.update(file=fi,row=row,elapsed_seconds=time.monotonic()-START,pid=os.getpid(),nice=os.getpriority(os.PRIO_PROCESS,0))
    tmp=statepath.with_suffix('.tmp');tmp.write_text(json.dumps(state,indent=2));tmp.replace(statepath)
    rt.status(f'Checkpoint release={fi} row={row}; cells={int(nt.sum())}; range bytes={state.get("current_bytes",0)}; job bytes={rt.guard()["job_bytes"]}.')

def hashes(ids):
    # Keyed BLAKE2b gives independent seeded Bernoulli halves for controls/targets.
    return np.array([[hashlib.blake2b(str(x).encode(),key=int(s).to_bytes(8,'little'),digest_size=1).digest()[0]&1 for x in ids] for s in SEEDS],bool)

rt.status('Attempt 3: three URLs authorized; lossless whitelist-only five-seed streaming starts. One numerical worker; BLAS=1; no raw data saved.')
for fi,url in enumerate(rt.URLS):
    if fi < int(cursors[:,0].min()):continue
    rt.guard()
    with Remote(url) as remote,h5py.File(remote,'r',rdcc_nbytes=32*1024**2) as f:
        obs=read_elem(f['obs']);var=read_elem(f['var'])
        assert list(var.index.astype(str))==meta['genes'];assert obs.index.is_unique
        labs=obs.target_gene.astype(str).to_numpy();ids=obs.index.astype(str).to_numpy()
        controls=set(ids[labs=='non-targeting'])
        digest=hashlib.sha256('\n'.join(sorted(controls)).encode()).hexdigest()
        assert digest=='8a78bbd979caa07d49fe9b88532d35a1be8b58325d5a003d8c4e27c55488b996'
        codes=np.array([code.get(x,-1) for x in labs]);
        assert not any(t in finalized for t in set(labs) if t!='non-targeting'), 'Completed target overlaps later release'
        if fi:codes[labs=='non-targeting']=-1
        allowed=np.flatnonzero(codes>=0)
        half=hashes(ids[allowed]); half_by_row=np.zeros((5,len(obs)),bool);half_by_row[:,allowed]=half
        ptr=f['X/indptr'][:];data=f['X/data'];indices=f['X/indices']
        assert f['X'].attrs['encoding-type']=='csr_matrix'
        receipt=dict(url=url,etag=remote.etag,size=remote.size,shape=list(map(int,f['X'].attrs['shape'])),allowed_cells=len(allowed),control_ids_sha256=digest,logical_expression_rows_read=0,heldout_expression_rows_read=0)
        print('OPEN',receipt,flush=True)
        begin=0 if fi>int(cursors[:,0].min()) else int(cursors[cursors[:,0]==fi,1].min())
        for r0 in range((begin//512)*512,len(obs),512):
            if time.monotonic()-START>LIMIT:checkpoint(fi,r0);raise TimeoutError('Wall limit')
            rt.guard();r1=min(r0+512,len(obs))
            rr=np.arange(r0,r1);cc=codes[rr]
            keep=cc>=0
            for k in np.flatnonzero(keep):
                j=cc[k]
                if (fi,int(rr[k]))<tuple(cursors[j]):keep[k]=False
            rr=rr[keep];cc=cc[keep]
            if len(rr):
                # Slice ONLY contiguous runs of allowed row IDs. Never bridge a forbidden row.
                cuts=np.r_[0,np.flatnonzero(np.diff(rr)>1)+1,len(rr)]
                mats=[]
                for lo,hi in zip(cuts[:-1],cuts[1:]):
                    rows=rr[lo:hi];p0=int(ptr[rows[0]]);p1=int(ptr[rows[-1]+1])
                    vals=data[p0:p1].astype(np.float64);ix=indices[p0:p1]
                    assert np.isfinite(vals).all() and np.all(vals>=0) and np.array_equal(vals,np.floor(vals))
                    mats.append(sp.csr_matrix((vals,ix,ptr[rows[0]:rows[-1]+2]-p0),shape=(len(rows),G)))
                x=sp.vstack(mats,format='csr');lib=np.asarray(x.sum(1)).ravel();assert np.all(lib>0)
                cpm=x.multiply((1e6/lib)[:,None]).tocsr();sq=cpm.copy();sq.data**=2
                for j in np.unique(cc):
                    take=cc==j;nt[j]+=int(take.sum())
                    for k,mat in enumerate([x,cpm,sq]):total[k,j]+=np.asarray(mat[take].sum(0)).ravel()
                    for s in range(5):
                        take=(cc==j)&half_by_row[s,rr];n[s,j]+=int(take.sum())
                        if np.any(take):
                            for k,mat in enumerate([x,cpm,sq]):a[k,s,j]+=np.asarray(mat[take].sum(0)).ravel()
                receipt['logical_expression_rows_read']+=len(rr)
            if ((r0//512)+1)%4==0:
                print('PROGRESS',fi,r1,int(nt.sum()),remote.bytes,flush=True)
            if ((r0//512)+1)%16==0:
                state['current_bytes']=remote.bytes;checkpoint(fi,r1)
        receipt.update(range_bytes=remote.bytes,requests=remote.requests)
        state['receipts']=[r for r in state['receipts'] if r['url']!=url]+[receipt]
        checkpoint(fi+1,0)
        subprocess.run([sys.executable,'-B',str(rt.ROOT/'finalize_release.py')],check=True)
        finalized={p.stem for p in (rt.ROOT/'finalized').glob('*.json')}
        del obs,var,half_by_row;gc.collect()
# Validate full newly streamed whitelist against prior full summaries, never all labels.
checks={}
for k,key in enumerate(['count_sum','cpm_sum','cpm_sq_sum']):
    active=np.array([t not in finalized for t in labels]);err=np.abs(total[k,active]-full[key][active]);checks[key]=dict(max_abs=float(err.max()),max_relative=float(np.max(err/np.maximum(np.abs(full[key][active]),1e-30))))
    np.testing.assert_allclose(total[k,active],full[key][active],rtol=2e-12,atol=1e-7)
assert len(finalized)==168
np.testing.assert_array_equal(nt,full['n_cells'])
state.update(complete=True,full_sum_checks=checks,n_cells=nt.tolist(),seeds=SEEDS.tolist(),seconds=time.monotonic()-START)
statepath.write_text(json.dumps(state,indent=2));rt.save('stream_receipt.json',state)
rt.status(f'Stream complete and full-summary equality verified: {nt.sum()} unique whitelist/control cells. {checks}')
