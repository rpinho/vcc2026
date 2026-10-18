"""Run downloaded AtlasShift numerical kernels only; no downloads or challenge inference."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
import importlib.util,sys,pathlib,json,hashlib
import numpy as np
root=pathlib.Path(__file__).parent
path=root/'vendor/atlas_model.py'
spec=importlib.util.spec_from_file_location('atlas_public',path);m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m)
out=root/'artifacts';out.mkdir(exist_ok=True)
rng=np.random.default_rng(7331)
n,g=400,40
depths=np.rint(np.linspace(11000,31000,n)).astype('int64')
p0=rng.dirichlet(np.full(g,3.0))
z=depths/depths.mean();loading=np.linspace(-.5,.5,g)
composition=p0[None,:]*np.exp((z[:,None]-1)*loading[None,:]);composition/=composition.sum(1,keepdims=True)
p=composition.mean(0);bulk=depths@composition/depths.sum()
template=rng.gamma(5,1,(n,g))*p0[None,:]
counts=m.dual_moment_counts(template,p,bulk,depths=depths,seed=20260910)
repeat=m.dual_moment_counts(template,p,bulk,depths=depths,seed=20260910)
one=template/template.sum(1,keepdims=True)
for _ in range(100):
 one*=p/one.mean(0);one/=one.sum(1,keepdims=True)
np.savez(out/'atlas_fixture.npz',counts=counts,repeat=repeat,depths=depths,p=p,bulk=bulk,template=template)
# Exercise aligned effects and weighted source fusion without target leakage.
genes=np.array(['T1','T2','G1','G2']);targets=np.array(['T1','T2'])
s=m.Source('fixture',targets,genes,np.array([[.05,.25,.4,.3],[.2,.05,.3,.45]]),np.array([.2,.2,.3,.3]),np.array([100,100]))
a,mask=m.aligned_effect(s,targets,genes,space='bulk_delta',common_subtract=1)
res={'source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'shape':list(counts.shape),'seed':20260910,'row_totals_exact':bool(np.array_equal(counts.sum(1),depths)),'deterministic':bool(np.array_equal(counts,repeat)),'mean_cpm_L1':float(abs((counts/depths[:,None]).mean(0)-p).sum()),'bulk_L1':float(abs(counts.sum(0)/counts.sum()-bulk).sum()),'one_moment_mean_cpm_L1':float(abs(one.mean(0)-p).sum()),'one_moment_bulk_L1':float(abs(depths@one/depths.sum()-bulk).sum()),'bulk_effect_fixture':a.tolist(),'nonnegative_integer':bool(np.issubdtype(counts.dtype,np.integer) and counts.min()>=0),'status':'synthetic software check; no biological score'}
(out/'atlas_fixture.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))
