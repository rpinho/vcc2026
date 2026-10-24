"""All-target public raw-count summaries and matched controls, on canonical genes.
Streaming 512 rows; writes only model_billm. No challenge treated outcomes.
"""
import common as rt
import numpy as np, pandas as pd, anndata as ad
from scipy import sparse
import gc, time
G = pd.read_csv(rt.REPO/'data/gene_names.csv').gene_name.to_numpy(str)
NT = 'non-targeting'

def align(x, genes):
    idx = pd.Index(genes).get_indexer(G); ok = idx >= 0
    out = np.zeros(x.shape[:-1]+(len(G),),np.float32)
    out[...,ok] = x[...,idx[ok]]
    return out, ok

def commit(name, family, targets, genes, counts, control, n, provenance):
    if len(set(genes)) != len(genes):
        unique=np.array(sorted(set(genes))); proj=sparse.csr_matrix((np.ones(len(genes)),(np.arange(len(genes)),pd.Index(unique).get_indexer(genes))),shape=(len(genes),len(unique)))
        counts=counts@proj; control=control@proj; genes=unique
    # counts and control both on native axis; normalize BEFORE mapping (no HVG renormalization).
    counts = counts/np.maximum(counts.sum(-1,keepdims=True),1e-30)
    control = control/np.maximum(control.sum(-1,keepdims=True),1e-30)
    response = np.log1p(50000*counts)-np.log1p(50000*control)
    response, measured = align(response,genes); base,_ = align(control,genes)
    if base.ndim == 1: base = np.broadcast_to(base,response.shape).copy()
    keep = (n>=30)&np.isfinite(response).all(1)&(counts.sum(1)>0)
    p=rt.ROOT/f'data/{name}.npz'
    np.savez_compressed(p,response=response[keep],control=base[keep],measured=measured,
        targets=np.asarray(targets,str)[keep],n_cells=n[keep],family=np.array(family),genes=G)
    rt.save(p.with_suffix('.json'),dict(name=name,family=family,rows=int(keep.sum()),genes=int(measured.sum()),
        sha256=rt.sha(p),bytes=p.stat().st_size,provenance=provenance,
        response='log1p(50000 * native pooled probability) perturbed minus matched control; zeros are masked when unmeasured'))
    print('SAVED',name,int(keep.sum()),p.stat().st_size,flush=True)

def raw(name):
    start=time.monotonic(); path=rt.REPO/f'public_data/{name}.h5ad'; a=ad.read_h5ad(path,backed='r')
    genes=a.var.gene_name.astype(str).to_numpy() if 'gene_name' in a.var else a.var_names.to_numpy(str)
    labels=a.obs['gene' if 'gene' in a.obs else 'target_gene'].astype(str).to_numpy()
    batches=a.obs['gem_group' if 'gem_group' in a.obs else 'batch'].astype(str).to_numpy()
    bc,bu=pd.factorize(batches); targets=np.array(sorted(set(labels)-{NT}))
    tc=pd.Index(np.concatenate(([NT],targets))).get_indexer(labels)
    sums=np.zeros((len(targets)+1,len(genes)),np.float64)
    cs=np.zeros((len(bu),len(genes)),np.float64); cn=np.zeros(len(bu),np.int64)
    mix=np.zeros((len(targets)+1,len(bu)),np.int64); n=np.bincount(tc,minlength=len(targets)+1)
    for lo in range(0,a.n_obs,512):
        hi=min(lo+512,a.n_obs); x=sparse.csr_matrix(a.X[lo:hi],dtype=np.float64)
        assert np.isfinite(x.data).all() and (x.data>=0).all() and np.equal(x.data,np.floor(x.data)).all()
        cols=np.arange(hi-lo); grp=sparse.csr_matrix((np.ones(len(cols)),(tc[lo:hi],cols)),shape=(len(targets)+1,len(cols)))
        sums += (grp@x).toarray(); np.add.at(mix,(tc[lo:hi],bc[lo:hi]),1)
        nt=tc[lo:hi]==0
        if nt.any():
            bg=sparse.csr_matrix((np.ones(nt.sum()),(bc[lo:hi][nt],np.arange(nt.sum()))),shape=(len(bu),nt.sum()))
            cs += (bg@x[nt]).toarray(); np.add.at(cn,bc[lo:hi][nt],1)
        if lo%51200==0: print(name,lo,'/',a.n_obs,flush=True)
    assert (cn>0).all(), 'No unmatched-batch fallback allowed'
    # Match control raw mean to each target's observed batch composition.
    control=(mix[1:]/n[1:,None])@(cs/cn[:,None]);a.file.close()
    commit(name,name,targets,genes,sums[1:],control,n[1:],dict(path=str(path.relative_to(rt.REPO)),
        sha256=rt.sha(path),seconds=time.monotonic()-start,batches=len(bu),controls_per_batch=cn.tolist(),
        matching='raw control mean mixed by target cell-count fraction over batches'))

def banks(name):
    path=rt.REPO/f'public_data/model_v4_signature_banks/{name}.npz'
    if name=='K562_GWPS_raw':path=rt.REPO/'rsi/artifacts/cd4_banks/K562_GWPS_raw.npz'
    z=np.load(path)
    n=z['n_cells']; control=z['matched_control_probability']; family='CD4' if name.startswith('CD4_') else ('K562' if name.startswith('K562') else name)
    commit(name,family,z['targets'],z['genes'],z['target_count_sums'],control,n,
        dict(path=str(path.relative_to(rt.REPO)),sha256=rt.sha(path),scope='all cached targets; HCT/HEK/CD4 caches have 437 targets, not their full original screens'))

def main():
    for name in ['Jurkat','HepG2','RPE1','K562','HCT116','HEK293T','H1','CD4_Rest','CD4_Stim8hr','CD4_Stim48hr','K562_GWPS_raw']:
        if (rt.ROOT/f'data/{name}.npz').exists(): continue
        rt.status('Preparing source summaries: '+name+'. Full training has NOT started. Next: exact MPS and CPU step benchmark, then enforce 6h cutoff.')
        (raw if name in ['Jurkat','HepG2','RPE1','K562'] else banks)(name); gc.collect()
    rt.save(rt.ROOT/'data/manifest.json',{p.stem:rt.json.loads(p.read_text()) for p in (rt.ROOT/'data').glob('*.json') if p.name!='manifest.json'})
    rt.status('Source summaries complete. Full training has NOT started. Run benchmark.py before train.py.')
if __name__=='__main__': main()
