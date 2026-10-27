"""Apply existing public closed-form cores to training168, source-only H1 exclusion."""
import runtime as rt
rt.enforce();rt.guard(True)
import numpy as np,pandas as pd,h5py,time,ast
z=np.load(rt.ROOT/'results/baseline.npz');targets=z['targets'];genes=z['genes'];Y=z['truth'].astype(float);P=z['realized'].astype(float);own=z['own_idx'];bulk=z['bulk'];n=len(targets)
ns=dict(np=np,pd=pd,ROOT=rt.REPO/'model_public',BANKS=rt.REPO/'public_data/model_v4_signature_banks',FAMILIES={'K562':['K562_GWPS_raw'],'H1':['H1'],'HCT116':['HCT116'],'HEK293T':['HEK293T'],'CD4':['CD4_Rest','CD4_Stim8hr','CD4_Stim48hr']},sha=rt.sha)
rt.definitions(rt.REPO/'model_public/linear.py',['source_panel','fit'],ns)
start=time.monotonic();ts,y,coverage,receipts=ns['source_panel'](genes,exclude='H1');read_seconds=time.monotonic()-start
linearfit=ns['fit'];slimns=dict(ns,ast=ast)
rt.definitions(rt.REPO/'model_public/slim.py',['upstream_solve','fit'],slimns)
with h5py.File(rt.REPO/'manifold_1026/inputs/string_embeddings.h5') as h:
    embeddings={str(t):h[str(t)][:].astype(float) for t in set(targets)|set(ts) if str(t) in h}
out={};entries=[]
keep=~np.isin(genes,targets);Ym=Y.copy();Ym[np.arange(n),own]=0;Pm=P.copy();Pm[np.arange(n),own]=0;E0=np.sum((Pm-Ym)**2)
def pds(a):
    aa=a[:,keep];yy=Ym[:,keep];c=aa @ yy.T/np.maximum(np.linalg.norm(aa,axis=1)[:,None]*np.linalg.norm(yy,axis=1)[None,:],1e-30)
    r=(c>np.diag(c)[:,None]).sum(1)+.5*((c==np.diag(c)[:,None]).sum(1)-1)
    return 1-r/(n-1)
for name,fit in [('rank10',linearfit),('SLIM',slimns['fit'])]:
    s=time.monotonic();basis,w,b,info=fit(y,ts,genes,embeddings=embeddings if name=='SLIM' else None)
    if name=='SLIM':
        good=np.array([t in embeddings for t in targets]);feat=np.zeros((n,w.shape[1]));feat[good]=np.stack([embeddings[t] for t in targets[good]])
    else:
        ix=pd.Index(genes).get_indexer(targets);good=ix>=0;feat=np.zeros((n,basis.shape[1]));feat[good]=basis[ix[good]]
    delta=b[None,:]+(feat @ w.T)@basis.T;delta[:,coverage.sum(1)==0]=0
    # Same declared logCP10k inversion as model_public.ResponseBank, using legal
    # destination pooled control because this job diagnoses the bulk channel.
    q=np.maximum((bulk[None,:]*10000+1)*np.exp(np.clip(delta,-8,8))-1,0);q/=q.sum(1,keepdims=True)
    a=np.log1p(50000*q)-np.log1p(50000*bulk);a[np.arange(n),own]=0
    e=np.sum((a-Ym)**2,1);pn=np.sum(a*a,1);yn=np.sum(Ym*Ym,1);cos=np.sum(a*Ym,1)/np.maximum(np.sqrt(pn*yn),1e-30);pdsv=pds(a)
    out[name]=dict(SSE_reduction_vs_B0=float(1-e.sum()/E0),error_ratio_zero=float(e.sum()/yn.sum()),macro_relative_error=float(np.mean(e/yn)),mean_cosine=float(cos.mean()),PDS168=float(pdsv.mean()),available_targets=int(good.sum()),fit_and_predict_seconds=time.monotonic()-s,source_read_seconds=read_seconds,fit=info)
    for i,t in enumerate(targets):entries.append(dict(method=name,target=t,sse=e[i],cosine=cos[i],pds=pdsv[i]))
    print(name,out[name],flush=True)
rt.save('results/cpu_baselines.json',dict(results=out,sources=receipts,training_targets=targets.tolist(),excluded_family='H1',destination_treated_rows_used_for_fit=0,note='Frozen source-only numerical cores, no H1 fitting. This is the existing cross-line adapter, not a published benchmark reproduction. Full STRING vectors already local. No downloaded assets. Bulk-control rather than arithmetic-control inversion explicitly matches this job channel.'))
pd.DataFrame(entries).to_csv(rt.ROOT/'results/cpu_baselines_per_target.csv',index=False)
rt.status('SLIM and rank10 source-only adapters measured on training168; no destination-treated fitting.')
