"""168-only shared-response diagnostic with nested target-LOO coefficients."""
import runtime as rt
rt.enforce();rt.guard(True)
import numpy as np, pandas as pd, time, gc
from functools import lru_cache
from threadpoolctl import threadpool_info
start=time.monotonic();n=168
rt.save('results/run.json',dict(pid=rt.os.getpid(),requested_nice=10,effective_nice=rt.os.getpriority(rt.os.PRIO_PROCESS,0),start_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),free_gib=rt.guard()))
pins={}
def pin(p):
    p=rt.Path(p);pins[str(p)]=rt.sha(p);return p
with np.load(pin(rt.PARENT/'results/baseline.npz')) as z:
    targets=z['targets'];genes=z['genes'];Y=z['truth'].astype(float);P=z['realized'].astype(float);own=z['own_idx'];bulk=z['bulk']
assert Y.shape==P.shape==(168,18077) and len(set(targets))==168
valid=own>=0;ii=np.arange(n)
def mask(a):
    a=np.array(a,copy=True);a[ii[valid],own[valid]]=0;return a
Y=mask(Y);P=mask(P);g=len(genes)
old=pd.read_csv(pin(rt.PARENT/'results/per_target.csv'))
assert old.target.tolist()==targets.tolist()
assert np.allclose(np.sum((P-Y)**2,1),old.sse,rtol=0,atol=1e-9)
old_analysis=rt.json.loads(pin(rt.PARENT/'results/analysis.json').read_text())
base=np.log1p(50000*bulk)
def project(a):
    q=np.maximum(np.expm1(base+a),0);q/=q.sum(1,keepdims=True)
    return mask(np.log1p(50000*q)-base)

def score(a):
    err=np.sum((a-Y)**2,1);yn=np.sum(Y*Y,1);pn=np.sum(a*a,1)
    cos=np.sum(a*Y,1)/np.maximum(np.sqrt(yn*pn),1e-30)
    keep=~np.isin(genes,targets);aa=a[:,keep];yy=Y[:,keep]
    c=aa @ yy.T/np.maximum(np.linalg.norm(aa,axis=1)[:,None]*np.linalg.norm(yy,axis=1)[None,:],1e-30)
    diag=np.diag(c)[:,None]
    rank=np.sum(c>diag,1)+.5*(np.sum(c==diag,1)-1)
    pds=1-rank/(n-1)
    B=np.sum((a-a.mean(0))**2)
    return dict(sse=err,yn=yn,pn=pn,cos=cos,pds=pds,B=B,similarity=c)

bst=score(P)
assert abs(bst['pds'].mean()-old_analysis['summary']['B0']['PDS168'])<1e-12
rt.status(f'B0 reproduced: E={bst["sse"].sum()/bst["yn"].sum():.9f}, cosine={bst["cos"].mean():.9f}, PDS={bst["pds"].mean():.9f}. Requested nice10, inherited effective nice20.')
# All new fitting uses only this masked training168 truth. Every row's own output
# gene is excluded in fit/evaluation, including each inner fold's shared mean.
YY=Y @ Y.T;PY=P @ Y.T;yy=np.diag(YY);py=np.diag(PY);pp=np.sum(P*P,1)
ys=Y.sum(0);Ysum=YY.sum(1);Psum=PY.sum(1);total=float(YY.sum())
Yown=np.zeros((n,n));Yown[:,valid]=Y[:,own[valid]]
ysown=Yown.sum(0)
@lru_cache(None)
def beta(excluded):
    """Fit on S, with inner predictor mean excluding j and excluded rows."""
    ex=np.array(excluded,dtype=int);tr=np.ones(n,bool);tr[ex]=False;L=int(tr.sum());d=L-1
    y_s=Ysum-YY[:,ex].sum(1);p_s=Psum-PY[:,ex].sum(1)
    snorm=total-2*Ysum[ex].sum()+YY[np.ix_(ex,ex)].sum()
    sown=ysown-Yown[ex].sum(0)
    mm=(snorm-2*y_s+yy-sown*sown)/(d*d)
    pm=(p_s-py)/d;ym=(y_s-yy)/d
    assert np.min(mm[tr])>0
    b=float(np.sum((py-pm)[tr])/np.sum(pp[tr]))
    c=float(np.sum((py-pm*ym/mm)[tr])/np.sum((pp-pm*pm/mm)[tr]))
    return b,c

M=mask((ys[None,:]-Y)/(n-1));mm=np.sum(M*M,1);pm=np.sum(P*M,1);ym=np.sum(Y*M,1)
Q=P-(pm/mm)[:,None]*M
Yperp=Y-(ym/mm)[:,None]*M
betas=np.array([beta((t,)) for t in range(n)])
preds={'B0':P,'a_mean':M,'b_mean_plus_B0':mask(M+betas[:,0,None]*P),'c_mean_plus_orthogonal':mask(M+betas[:,1,None]*Q)}
shared=dict(projection_energy_fraction=float(np.sum(ym*ym/mm)/np.sum(yy)),mean_truth_cosine=float(np.mean(ym/np.sqrt(yy*mm))),mean_B0_cosine=float(np.mean(pm/np.sqrt(pp*mm))),B0_projection_energy_fraction=float(np.sum(pm*pm/mm)/np.sum(pp)),orthogonal_E=float(np.sum((Q-Yperp)**2)/np.sum(Yperp*Yperp)),orthogonal_mean_cosine=float(np.mean(np.sum(Q*Yperp,1)/np.sqrt(np.sum(Q*Q,1)*np.sum(Yperp*Yperp,1)))),orthogonal_optimal_global_scalar_descriptive=float(np.sum(Q*Yperp)/np.sum(Q*Q)),truth_energy=float(yy.sum()),orthogonal_truth_energy=float(np.sum(Yperp*Yperp)))
rt.status(f'Shared projection explains {100*shared["projection_energy_fraction"]:.2f}% of truth energy; B0/shared mean cosine={shared["mean_B0_cosine"]:.4f}; unscaled orthogonal E={shared["orthogonal_E"]:.4f}. Reconstructing non-H1 family features.')
# Exact previous bank aggregation/realization; restrict selected rows to 168.
manifest=rt.json.loads(pin(rt.REPO/'model_billmB/data/manifest.json').read_text())
family_sum={};family_den={};coverage={};source_rows=[]
for name,m in sorted(manifest.items()):
    if m['family']=='H1':continue
    rt.guard(True)
    path=pin(rt.REPO/f'model_billmB/data/{name}.npz')
    assert pins[str(path)]==m['sha256']
    with np.load(path) as z:
        ind=pd.Index(z['targets']).get_indexer(targets);good=np.flatnonzero(ind>=0)
        axis=pd.Index(z['genes']).get_indexer(genes);assert (axis>=0).all()
        meas=z['measured'][axis];r=z['response']
        f=m['family'];family_sum.setdefault(f,np.zeros((n,g),np.float32));family_den.setdefault(f,np.zeros((n,g),np.float32));coverage.setdefault(f,np.zeros(n,bool))
        for i in good:
            family_sum[f][i]+=r[ind[i],axis]*meas;family_den[f][i]+=meas
        coverage[f][good]=True
        source_rows.append(dict(bank=name,family=f,selected_training_rows=int(len(good))))
    del r;gc.collect()
    print('bank',name,'selected',len(good),flush=True)
family_names=sorted(family_sum);Fs=[]
for f in family_names:
    a=project(family_sum[f]/np.maximum(family_den[f],1e-8)).astype(np.float32).astype(float)
    a[~coverage[f]]=P[~coverage[f]]
    Fs.append(a)
R=np.stack(Fs)-P[None,:,:]
del Fs,family_sum,family_den;gc.collect()
K=len(family_names)
G=np.einsum('kig,lig->ikl',R,R,optimize=True)
RY=np.einsum('kig,ig->ik',R,Y,optimize=True);RP=np.einsum('kig,ig->ik',R,P,optimize=True)
RYall=np.einsum('kig,jg->ikj',R,Y,optimize=True)
w_d=[];w_e=[];out_d=np.empty_like(P);out_e=np.empty_like(P);inner_betas=np.full((n,n),np.nan)
for t in range(n):
    tr=ii!=t;gram=G[tr].sum(0);reg=max(np.trace(gram)/K,1e-12)
    wd=np.linalg.solve(gram+reg*np.eye(K),np.sum((RY-RP)[tr],0))
    # Base c used to fit family residuals is inner cross-fitted: its beta is
    # trained without BOTH t and j; its own training means also exclude k.
    ib=np.zeros(n)
    for j in ii[tr]:ib[j]=beta(tuple(sorted((int(t),int(j)))))[1]
    inner_betas[t,tr]=ib[tr]
    pmj=(Psum-PY[:,t]-py)/(n-2)
    s2=total-2*Ysum[t]+yy[t]
    ysj=Ysum-YY[:,t];ownsum=ysown-Yown[t]
    mmj=(s2-2*ysj+yy-ownsum**2)/(n-2)**2
    RM=(RYall.sum(2)-RYall[:,:,t]-RY)/(n-2)
    RQ=RP-(pmj/mmj)[:,None]*RM
    rhs=RY-RM-ib[:,None]*RQ
    we=np.linalg.solve(gram+reg*np.eye(K),rhs[tr].sum(0))
    out_d[t]=P[t]+wd @ R[:,t]
    out_e[t]=preds['c_mean_plus_orthogonal'][t]+we @ R[:,t]
    w_d.append(wd);w_e.append(we)
preds['d_family_correction']=mask(out_d);preds['e_c_plus_family']=mask(out_e)
# Exact pre-projection vector equations are primary. Previous-job feasibility
# projection is a separately reported sensitivity, never count rendering.
projected={name:project(a) for name,a in preds.items() if name!='B0'}
projected['B0']=P
prior_d=score(projected['d_family_correction'])
for key,val in [('sse',prior_d['sse'].sum()),('PDS168',prior_d['pds'].mean()),('centered_energy',prior_d['B'])]:
    wanted=old_analysis['summary']['LOO_family_transfer'][key]
    assert np.isclose(val,wanted,rtol=2e-8,atol=2e-7),(key,val,wanted)
rt.status('PASS previous family-correction replay, including feasibility projection. All five new arms generated with nested LOO; computing 2000 paired target bootstraps.')
rng=np.random.default_rng(20261026);boots=rng.integers(n,size=(2000,n))
summary={};rows=[];similarities={}
for space,arms in [('vector',preds),('projected',projected)]:
    summary[space]={}
    for name,a in arms.items():
        st=score(a);e=st['sse'];T=st['yn'];pdsv=st['pds'];b=bst
        boot_E=e[boots].sum(1)/T[boots].sum(1)
        boot_PDS=pdsv[boots].mean(1)
        dE=(e-b['sse'])[boots].sum(1)/T[boots].sum(1)
        dP=(pdsv-b['pds'])[boots].mean(1)
        result=dict(E=float(e.sum()/T.sum()),mean_cosine=float(st['cos'].mean()),PDS=float(pdsv.mean()),B=float(st['B']),B_relative_B0=float(st['B']/b['B']),delta_E=float((e.sum()-b['sse'].sum())/T.sum()),delta_PDS=float((pdsv-b['pds']).mean()),E_CI95=np.quantile(boot_E,[.025,.975]).tolist(),PDS_CI95=np.quantile(boot_PDS,[.025,.975]).tolist(),delta_E_CI95=np.quantile(dE,[.025,.975]).tolist(),delta_PDS_CI95=np.quantile(dP,[.025,.975]).tolist())
        summary[space][name]=result;similarities[space+'__'+name]=st['similarity']
        for i in range(n):rows.append(dict(space=space,arm=name,target=targets[i],sse=e[i],truth_norm2=T[i],prediction_norm2=st['pn'][i],cosine=st['cos'][i],pds=pdsv[i]))
        print(space,name,result,flush=True)
# Only small summaries/cross-products are retained; no full challenger vectors.
pd.DataFrame(rows).to_csv(rt.ROOT/'results/per_target.csv',index=False)
np.savez_compressed(rt.ROOT/'results/sufficient.npz',YY=YY,PY=PY,PP=pp,G=G,RY=RY,RP=RP,RYall=RYall,Yown=Yown,targets=targets,genes=genes,**similarities)
rt.save('results/coefficients.json',dict(beta_b=betas[:,0].tolist(),beta_c=betas[:,1].tolist(),inner_beta_c=[[None if np.isnan(v) else float(v) for v in row] for row in inner_betas],family_names=family_names,family_d=np.asarray(w_d).tolist(),family_e=np.asarray(w_e).tolist()))
rt.save('results/analysis.json',dict(shared=shared,summary=summary,bootstrap=dict(replicates=2000,seed=20261026,unit='paired targets',method='percentile; resample held-out per-target losses/ranks; fixed 168-target retrieval candidate panel and fitted predictions; no bootstrap refits'),equations=dict(a='m_-t',b='m_-t + beta_b[-t] P_t',c='m_-t + beta_c[-t] (P_t - proj_m P_t)',d='P_t + sum_f w_d[-t,f] (F_ft-P_t)',e='c_t + sum_f w_e[-t,f] (F_ft-P_t)'),source_rows=source_rows,training_rows=168,heldout_79_inputs_opened=False,rendering=False,official_MSE='not measurable; vcc_mse_instr_1026 found no instrument',seconds=time.monotonic()-start,pins=pins,free_gib=rt.guard(),blas=threadpool_info()))
rt.status(f'Analysis complete in {time.monotonic()-start:.1f}s; primary vector and projected sensitivity tables saved. Independent verification next.')
