"""Frozen A/B mean construction and paired, dispersed integer-count emission."""
import pathlib,numpy as np,scipy.sparse as sp
ROOT=pathlib.Path(__file__).resolve().parent

def mean_vector(ctrl,shift,genes,target,residual=None):
 effect=shift.copy()
 if residual is not None:
  cpm=ctrl/ctrl.sum()*1e6
  effect+=np.where(cpm>=5,residual,0)
 mu=ctrl*np.exp2(effect)
 own=np.flatnonzero(genes==target)
 if len(own):mu[own[0]]=0.1*ctrl[own[0]]
 return mu

def emit(X,ctrl,mu,context_index,target_index,n=400):
 rng=np.random.default_rng(np.random.SeedSequence([2026,context_index,target_index]))
 ix=rng.integers(0,X.shape[0],n)
 out=X[ix].tocsr().astype(float)
 scale=np.divide(mu,ctrl,out=np.zeros_like(mu),where=ctrl>0)
 out.data*=scale[out.indices];low=np.floor(out.data)
 out.data=low+(rng.random(len(low))<out.data-low)
 out.eliminate_zeros();out.sort_indices()
 return out.astype(np.int32)
