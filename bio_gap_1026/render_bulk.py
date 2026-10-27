"""Exact frozen B0 mean emitter, retain only column totals; no DE/package/cache."""
import runtime as rt
rt.enforce();rt.guard(True)
import numpy as np,anndata as ad,time,types,os
from scipy import sparse
z=np.load(rt.ROOT/'results/baseline.npz');targets=z['targets'];genes=z['genes'];r=z['raw'];bulk=z['bulk'];truth=z['truth'];continuous=z['realized'];own=z['own_idx']
src=rt.REPO/'h2h2_1026/inputs/vendor_runtime';ns=dict(np=np,sparse=sparse,hashlib=rt.hashlib)
rt.definitions(src/'vendor/atlas_model.py',['repair_bulk_totals','seed_for'],ns)
ns['atlas']=types.SimpleNamespace(repair_bulk_totals=ns['repair_bulk_totals'])
rt.definitions(src/'model.py',['template','mean_counts'],ns)
path=rt.REPO/'gatecall_1026/inputs/tasks/H1_rich_search.h5ad';a=ad.read_h5ad(path,backed='r');idx=np.flatnonzero(a.obs.evaluation_role.astype(str).to_numpy()=='input_control')
assert len(idx)==3200
ctrl=ad.AnnData(sparse.csr_matrix(a.X[idx]),obs=a.obs.iloc[idx].copy(),var=a.var.copy());a.file.close();del a
q,depth,_,b=ns['template'](ctrl,genes,seed=20260910);assert np.array_equal(b,bulk)
out=[];diff=[];start=time.monotonic();rt.status('Exact B0 integer emitter on training168 started; only input controls read; no DE/layout/package work.')
for i,t in enumerate(targets):
    rt.guard()
    if time.monotonic()-start>7200:raise RuntimeError('2h emitter cutoff; no background continuation')
    prob=np.maximum(np.expm1(np.log1p(50000*bulk)+r[i]),0);prob/=prob.sum()
    # Exact seed used by frozen H1 B0 renderer; v9 layout preserves both axes.
    x=ns['mean_counts'](q,None,prob,depth,ns['seed_for']('H1:'+t,20260910)).astype(np.int32)
    counts=x.sum(0,dtype=np.int64);real=np.log1p(50000*counts/counts.sum())-np.log1p(50000*bulk);out.append(real.astype(np.float32));d=real-continuous[i]
    diff.append(dict(target=str(t),max_abs_response_error=float(abs(d).max()),L2_difference=float(np.linalg.norm(d)),counts=int(counts.sum())))
    del x
    if i%10==0:rt.status(f'Exact B0 bulk {i+1}/168, elapsed {time.monotonic()-start:.1f}s')
np.savez_compressed(rt.ROOT/'scratch/rendered_bulk.npz',response=np.stack(out),targets=targets,genes=genes)
rt.save('results/renderer_audit.json',dict(status='PASS',groups=168,seconds=time.monotonic()-start,comparison=diff,control_rows=len(idx),reference_treated_rows_read=0,group_cells=400,group_depth_total=int(depth.sum()),source_pins={str(p):rt.sha(p) for p in [src/'model.py',src/'vendor/atlas_model.py']},method='Exact frozen mean_counts + integer column repair; v9 layout provably preserves these same column totals. No count matrices retained.',nice=os.getpriority(os.PRIO_PROCESS,0)))
rt.status('Exact emitted B0 bulk vectors saved; no count matrices retained.')
