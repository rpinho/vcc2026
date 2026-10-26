import local_runtime as j
j.enforce_io();j.guard()
import numpy as np,polars as pl
from cell_eval2.score import _replicate_entries
from cell_eval2.scoring import score_one
import score_official as sc
SIX=sc.SIX;ARMS=['B0','B1','B2','B3','B4','B5'];N=20000;SEED=20261026

def aggregate(cols,mat,indices=None):
 x=mat if indices is None else mat[indices];axis=0 if indices is None else 1
 with np.errstate(invalid='ignore',divide='ignore'):
  out={m:np.nanmean(x[...,i],axis=axis) for i,m in enumerate(cols)}
  a=x[...,cols.index('expr_mse_unbiased_capped')];b=x[...,cols.index('expr_distance_unbiased')];ok=np.isfinite(a)&np.isfinite(b)
  out[SIX[0]]=np.where(ok,a,0).sum(axis=axis)/np.where(ok,b,0).sum(axis=axis)
 return out

def main():
 indices=np.random.default_rng(SEED).integers(79,size=(N,79));sample={};records={};draws={};diag={};missing={}
 with np.load(j.ROOT/'results/heldout_truth.npz') as z:y=z['response'];targets=z['targets']
 for arm in ARMS:
  rec=j.json.loads((j.ROOT/f'results/official_{arm}.json').read_text());records[arm]=rec
  with np.load(j.ROOT/f'results/target_metrics_{arm}.npz') as z:cols=z['columns'].tolist();mat=z['values'];assert np.array_equal(z['targets'],targets)
  raw=aggregate(cols,mat);assert all(abs(raw[m]-rec['raw'][m])<1e-12 for m in SIX)
  anchors=rec['anchors'];entries=_replicate_entries(anchors['baseline'],pl.DataFrame({'metric':list(anchors['replicate']),'replicate':list(anchors['replicate'].values())}))
  # Small batches avoid 20k x79 xmetric temporary arrays and swap.
  vectors={m:np.empty(N) for m in SIX}
  for lo in range(0,N,200):
   val=aggregate(cols,mat,indices[lo:lo+200])
   for m in SIX:vectors[m][lo:lo+200]=val[m]
  scaled=np.column_stack([[score_one(float(v),entries[m].base,entries[m].scoring) for v in vectors[m]] for m in SIX]);sample[arm]={SIX[3]:vectors[SIX[3]],SIX[1]:vectors[SIX[1]],'score_avg':scaled.mean(1)}
  for m,v in sample[arm].items():draws[arm+'__'+m]=v
  with np.load(j.ROOT/f'results/predictions_{arm}.npz') as z:p=z['response']
  dot=np.einsum('ij,ij->i',p,y);cos=dot/np.maximum(np.linalg.norm(p,axis=1)*np.linalg.norm(y,axis=1),1e-30);sse=np.sum((p-y)**2,axis=1)
  diag[arm]=dict(cosine_mean=float(cos.mean()),summed_sse=float(sse.sum()),targets=targets.tolist(),cosine=cos.tolist(),sse=sse.tolist())
  missing[arm]=dict(finite_targets={m:int(np.isfinite(mat[:,cols.index(m)]).sum()) for m in SIX if m in cols},missing_STRING={str(t):{m:None if not np.isfinite(mat[i,cols.index(m)]) else float(mat[i,cols.index(m)]) for m in SIX if m in cols} for i,t in enumerate(targets) if t in ['KDR','SOX2']})
 comparisons={}
 for a,b in [('B1','B0'),('B2','B0'),('B1','B4'),('B1','B5'),('B2','B4'),('B2','B5'),('B1','B3'),('B2','B3'),('B4','B0'),('B5','B0')]:
  c={}
  for m in [SIX[3],SIX[1],'score_avg']:
   point=records[a]['score_avg']-records[b]['score_avg'] if m=='score_avg' else records[a]['raw'][m]-records[b]['raw'][m]
   delta=sample[a][m]-sample[b][m];lo,hi=np.quantile(delta,[.025,.975],method='linear');c[m]=dict(delta=float(point),ci95=[float(lo),float(hi)],pass_positive=bool(point>0 and lo>0))
  c['all_three_pass']=all(v['pass_positive'] for v in c.values());comparisons[a+' minus '+b]=c
 passes={a:all(comparisons[a+' minus '+b]['all_three_pass'] for b in ['B0','B4','B5']) for a in ['B1','B2']}
 verdict='ADOPT' if any(passes.values()) else 'REJECT — permanently close manifold/STRING'
 np.savez_compressed(j.ROOT/'results/bootstrap_draws.npz',indices=indices,**draws)
 j.save('results/analysis.json',dict(verdict=verdict,arms_pass=passes,comparisons=comparisons,arms={a:dict(raw=r['raw'],scaled=r['scaled'],score_avg=r['score_avg']) for a,r in records.items()},response_diagnostics=diag,finite_cohorts_and_missing_STRING=missing,bootstrap=dict(draws=N,seed=SEED,method='paired targets; fixed panel PDS and anchors; official finite means and ratio-of-sums then nonlinear scaling')))
 print(verdict);print(j.json.dumps(comparisons,indent=2))
if __name__=='__main__':main()
