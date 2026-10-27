"""Bounded follow-up: maximize nested alignment subject to inner PDS loss <=.005.

All original arms retained. Finer small weights are specified before computing
this follow-up's outer results. No outer-target PDS/truth enters selection.
"""
import analyze as c
import numpy as np, pandas as pd, time
rt=c.rt
CONSERVATIVE='--conservative' in rt.sys.argv
WEIGHTS=(0.,.005,.01,.02,.03,.05,.075,.1,.15,.2,.35,.5,.75,1.)
start=time.monotonic()
rt.save('results/constrained_run.json',dict(pid=rt.os.getpid(),free_gib=rt.guard(True)))
prior=rt.json.loads((rt.ROOT/'results/analysis.json').read_text())
oldchoices=rt.json.loads((rt.ROOT/'results/choices.json').read_text())
c.load(); c.PINS.update(prior['pins'])
with np.load(rt.ROOT/'results/source_kernel.npz') as z: K=z['K']
saved=np.load(rt.ROOT/'results/coefficients.npz')
aa={arm:saved[arm+'__a'] for arm in oldchoices}
bb={arm:saved[arm+'__b'] for arm in oldchoices}
newnames=('B0_knn_PDS_CI','B0_map_PDS_CI') if CONSERVATIVE else ('B0_knn_PDS','B0_map_PDS')
for arm in newnames:
    aa[arm]=np.zeros(c.N); bb[arm]=np.zeros((c.N,c.N)); oldchoices[arm]=[]
YYr=c.Y[:,c.keep]@c.Y[:,c.keep].T
PYr=c.P[:,c.keep]@c.Y[:,c.keep].T
truthnorm=np.linalg.norm(c.Y[:,c.keep],axis=1)
ii=np.arange(c.N); newrows=[]

def pds(scores):
    diag=np.diag(scores)[:,None]
    ranks=(scores>diag).sum(1)+.5*((scores==diag).sum(1)-1)
    return 1-ranks/(len(scores)-1)

def select(inner,rows):
    B0=c.baseline_B(rows); basecross=c.py[rows].sum()
    base_scores=PYr[np.ix_(rows,rows)]/truthnorm[rows][None,:]
    basepds=pds(base_scores)
    best=None; log=[]
    for h,b in enumerate(inner):
        cross,_,Bq=c.cross_energy(np.zeros(len(rows)),b,rows)
        s=np.sqrt(B0/Bq)
        removed=np.sum(b*c.Yown[:,rows].T,axis=1)
        sumP=np.zeros(c.N); sumP[rows]=1
        Bcross=float(np.sum(b*c.PY[rows])-(sumP @ c.PY @ b.sum(0)-(sumP @ c.Pown[:,rows])@removed)/len(rows))
        source_scores=(b @ YYr)[:,rows]/truthnorm[rows][None,:]
        for w in WEIGHTS:
            spread=(1-w)**2*B0+w*w*s*s*Bq+2*(1-w)*w*s*Bcross
            fs=float(np.sqrt(B0/max(spread,1e-20)))
            objective=float(((1-w)*basecross+w*s*cross.sum())*fs/c.yy[rows].sum())
            differences=pds((1-w)*base_scores+w*s*source_scores)-basepds
            dp=float(differences.mean())
            lower=float(dp-1.96*differences.std(ddof=1)/np.sqrt(len(rows)))
            feasible=bool((lower if CONSERVATIVE else dp)>=-.005-1e-15)
            log.append(dict(h=h,weight=w,objective=objective,base_scale=float(s),final_scale=fs,inner_delta_PDS=dp,inner_PDS_lower_N95=lower,feasible=feasible))
            if feasible and (best is None or objective>best['objective']+1e-14):
                best=dict(h=h,weight=w,base_scale=float(s),final_scale=fs,objective=objective,inner_delta_PDS=dp,inner_PDS_lower_N95=lower)
    assert best is not None
    return best,log

rt.status('Final bounded follow-up fixed before its evaluation: '+('lower normal 95% bound on inner delta PDS' if CONSERVATIVE else 'inner delta PDS')+' >=-.005; weights '+str(WEIGHTS)+'. All earlier arms retained; no further search after conservative follow-up.')
for t in range(c.N):
    rt.guard(); rows=ii[ii!=t]
    ik=np.zeros((len(c.KS),len(rows),c.N)); im=np.zeros((len(c.HYPERS),len(rows),c.N))
    for pos,j in enumerate(rows):
        train=ii[(ii!=t)&(ii!=j)]
        ik[:,pos]=c.knn(K,train,[j])[:,0]
        im[:,pos]=c.maps(K,train,[j])[:,0]
    for arm,inner,outer in ((newnames[0],ik,c.knn(K,rows,[t])),(newnames[1],im,c.maps(K,rows,[t]))):
        ch,log=select(inner,rows)
        aa[arm][t],bb[arm][t]=c.apply(ch,outer)
        ch['target']=str(c.targets[t])
        ch['hyper']=dict(k=c.KS[ch['h']]) if 'knn' in arm else dict(zip(('lambda','rank'),c.HYPERS[ch['h']]))
        oldchoices[arm].append(ch)
        for item in log: newrows.append(dict(outer=t,arm=arm,**item))
    if t%8==0 or t==c.N-1:
        rt.status(f'PDS-constrained nested LOO {t+1}/{c.N}, elapsed {time.monotonic()-start:.1f}s.')
oldlog=pd.read_csv(rt.ROOT/'results/inner_objectives.csv',low_memory=False)
pd.concat([oldlog,pd.DataFrame(newrows)],ignore_index=True).to_csv(rt.ROOT/'results/inner_objectives.csv',index=False)
rt.save('results/choices.json',oldchoices)
np.savez_compressed(rt.ROOT/'results/coefficients.npz',**{name+'__a':aa[name] for name in aa},**{name+'__b':bb[name] for name in bb})
c.evaluate(aa,bb)
updated=rt.json.loads((rt.ROOT/'results/analysis.json').read_text())
updated.update({key:prior[key] for key in ('constrained_blends',) if key in prior})
updated['conservative_blends' if CONSERVATIVE else 'constrained_blends']=dict(arms=list(newnames),weights=list(WEIGHTS),
    inner_PDS_constraint=-.005,inner_candidate_panel='167 outer training targets, never outer t',
    constraint='lower normal 95% confidence bound' if CONSERVATIVE else 'point mean',
    motivation='Earlier blends fail final PDS; bounded same-family follow-up. Additional arms exploratory.',seconds=time.monotonic()-start)
updated['original_analysis_seconds']=prior.get('original_analysis_seconds',prior['seconds'])
rt.save('results/analysis.json',updated)
rt.status('Constrained blend follow-up complete; independent verification next.')
