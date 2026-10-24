"""Fail-closed feasibility probes BEFORE a full 900-group build or vcc prep.
This script never creates a package. It compares layout against its own learned-profile baseline.
"""
import os,sys,types,gc,time,json,hashlib
from pathlib import Path
REPO=Path(__file__).resolve().parent.parent;R=REPO/'model_billmB'
sys.path.insert(0,str(R));import common as job
for k in ['RAYON_NUM_THREADS','TOKIO_WORKER_THREADS']:os.environ[k]='2'
sys.path.insert(0,str(REPO/'model_billmB'));import io_guard;io_guard.install()
import shutil
rt=types.ModuleType('runtime');rt.ROOT=R;rt.json=json;rt.sha=job.sha;rt.save=job.save;rt.shutil=shutil
rt.guard=lambda *a,**k:None
sys.modules['runtime']=rt
sys.path.insert(0,str(REPO/'model_v7'));sys.path.insert(0,str(REPO/'model_v5'))
import numpy as np,pandas as pd,anndata as ad,torch
from scipy import sparse
import model,v5,build_v5 as old
assert job.sha(R/'artifacts/promoter_pairs.csv')==job.sha(REPO/'model_v5/artifacts/promoter_pairs.csv')
from infer import Ensemble
from network import blend
sys.path.insert(0,str(REPO/'model_v9'))
from de import prepare,Reference
from layout import emit
from audit_worker import audit
from fractions import Fraction

torch.set_num_threads(2);torch.set_num_interop_threads(1)
def bounded_template(a,genes,cells=400,pool_k=4,seed=20260910):
 assert np.array_equal(a.var_names,genes)
 X=sparse.csr_matrix(a.X);depth=np.asarray(X.sum(1,dtype=np.float64)).ravel()
 p=np.asarray(X.T@(1/depth)).ravel()/X.shape[0]
 b=np.asarray(X.sum(0,dtype=np.float64)).ravel();b/=b.sum()
 selected=np.random.default_rng(seed).choice(len(depth),cells*pool_k,replace=False);selected=selected[np.argsort(depth[selected],kind='stable')]
 q=X[selected].toarray().astype(np.float64)/depth[selected,None];q=q.reshape(cells,pool_k,len(genes)).mean(1)
 d=np.rint(depth[selected].reshape(cells,pool_k).mean(1)).astype(np.int64)
 return q,d,p,b
old.template=bounded_template;old.ARGS=types.SimpleNamespace(k=18533)
def main():
 assert json.loads((R/'artifacts/acceptance.json').read_text())['status']=='PASS', 'SEARCH acceptance required'
 full='--full' in sys.argv
 writer=None
 if full:
  pre=json.loads((REPO/'model_billmB/renderer_preflight.json').read_text())
  assert len(pre['records'])==6 and all(r['all_gates_pass'] for r in pre['records']), 'Six probes must pass before full construction'
  from io_counts import Writer
  (R/'out').mkdir(exist_ok=True)
 genes=pd.read_csv(REPO/'data/gene_names.csv').gene_name.to_numpy(str);targets=pd.read_csv(REPO/'data/pert_counts.csv').target_gene.to_numpy(str)
 saved=json.loads((REPO/'model_v5/artifacts/v5b_build.json').read_text());records=[]
 assert job.sha(REPO/'model_v5/model.py')==saved['model_sha256']
 assert job.sha(REPO/'model_v5/v5.py')==saved['v5_sha256']
 fixed=[REPO/'data/gene_names.csv',REPO/'data/pert_counts.csv',REPO/'model_v5/artifacts/v5b_build.json',REPO/'model_v8/artifacts/candA_gates.json']
 fixed += [REPO/f'data/context_{c}.h5ad' for c in 'ABC']
 fixed += [REPO/f for f in ['model_v5/model.py','model_v5/v5.py','model_v5/build_v5.py','model_v9/layout.py','model_v9/audit_worker.py','model_v9/de.py']]
 fixed += [R/n for n in ['renderer_preflight.py','common.py','dataset.py','infer.py','network.py','data/manifest.json','artifacts/weights_verified.json','artifacts/promoter_pairs.csv']]
 fixed += sorted((R/'data').glob('*.npz'))+sorted((R/'banks').glob('*.npz'))+sorted((R/'checkpoints/production').glob('*_weights.pt'))
 signature={str(p.relative_to(REPO)):job.sha(p) for p in fixed}
 input_sha=hashlib.sha256(json.dumps(signature,sort_keys=True).encode()).hexdigest()
 cache=R/'artifacts/validation_groups'/input_sha;cache.mkdir(parents=True,exist_ok=True)
 job.save(cache/'inputs.json',signature)
 if full:writer=Writer(R/'out/billmB.h5ad',genes)
 # Fixed first-two-target probes in each context; no fitness-based choice.
 for ci,c in enumerate('ABC'):
  print('RENDERER CONTEXT',c,flush=True);a=ad.read_h5ad(REPO/f'data/context_{c}.h5ad');assert a.shape==(18400,18533)
  out=R/f'artifacts/ref_{c}'
  if not (out/'cdf_scale.npy').exists():prepare(a.X,out,np.random.default_rng(20261018+ci))
  ref=Reference(c)
  q,d,L,Lc,mask,info=old.profile('b',targets,genes,a,c,ci,print)
  bf=v5.from_logbulk(v5.rescale(L,saved['contexts'][ci]['alpha'],Lc))
  b=np.asarray(a.X.sum(0,dtype=np.float64)).ravel();b/=b.sum();ensemble=Ensemble('production');assert np.array_equal(genes,ensemble.src.genes)
  pred=ensemble.predict(targets,b);new=v5.from_logbulk(Lc+blend(v5.logbulk(bf)-Lc,pred['learned'],.3))
  print('LEARNED PROFILE READY',c,flush=True)
  sources=[model.aligned_bank(n,targets,genes,'log2fc') for n in model.load_banks(None,True)];cd=model.cd4_family(targets,genes,'log2fc')
  if cd is not None:sources.append(cd)
  net=np.zeros(L.shape);zs=np.zeros(L.shape);sign=np.sign(L-Lc)
  for e,m,v in sources:
   z=np.where(m&np.isfinite(v),e/np.sqrt(np.maximum(v,1e-12)),0);sig=np.abs(z)>1.96
   agree=sig*(np.sign(z)==sign);oppose=sig*(np.sign(z)==-sign)
   net+=agree.astype(float)-oppose.astype(float);zs+=np.minimum(np.abs(z),100)*(np.sign(z)==sign)
  confidence=net*1e6+zs
  del ensemble,a,sources,cd;gc.collect()
  budget_records=json.loads((REPO/'model_v8/artifacts/candA_gates.json').read_text())['records'];budgets={(r['context'],r['target']):r['n_calls'] for r in budget_records}
  for i in (range(300) if full else [0,4]):
   counts=cache/f'{c}_{i:03d}.npz';checkpoint=cache/f'{c}_{i:03d}.json'
   if counts.exists() and checkpoint.exists():
    cached=json.loads(checkpoint.read_text());assert cached['inputs_sha256']==input_sha and cached['npz_sha256']==job.sha(counts)
    rec=cached['record'];assert rec['all_gates_pass'] and rec['context']==c and rec['index']==i and rec['target']==targets[i]
    x=sparse.load_npz(counts).toarray();assert hashlib.sha256(x.astype('<i4').tobytes()).hexdigest()==rec['counts_sha256']
    print('REUSE RENDERER GROUP',c,i,flush=True)
   else:
    seed=model.atlas.seed_for(c+':'+targets[i],20260910)
    original=model.emit(q,None,new[i],d,seed,'anchor_single_mean').astype(np.int32)
    assert np.array_equal(original.sum(1),d)
    own=int(np.flatnonzero(genes==targets[i])[0]);x,called,layout_info=emit(original,ref,confidence[i],budgets[c,targets[i]],own)
    rec,details=audit(x,original,c,confidence[i],own,called)
    rec.update(context=c,index=i,target=str(targets[i]),same_profile_comparator='v5b integer renderer on identical learned pooled profile',JK_dominance=Fraction(rec['jk'])>=Fraction(rec['v5b_jk']),v5b_preservation='INAPPLICABLE: learned pooled profile changes',**layout_info)
    rec['all_gates_pass']=bool(rec['JK_dominance'] and rec['mean_certificate']['mean_L1_upper']<.015 and rec['mean_certificate']['bulk_L1_upper']<.0015)
    rec['baseline_counts_sha256']=hashlib.sha256(original.astype('<i4').tobytes()).hexdigest()
    if rec['all_gates_pass']:
     partial=counts.with_suffix('.partial.npz');sparse.save_npz(partial,sparse.csr_matrix(x));partial.replace(counts)
     job.save(checkpoint,dict(inputs_sha256=input_sha,npz_sha256=job.sha(counts),record=rec))
   records.append(rec)
   if writer is not None and rec['all_gates_pass']:writer.append(x,c,targets[i])
   job.save(REPO/('model_billmB/full_gates.json' if full else 'model_billmB/renderer_preflight.json'),dict(status='PROGRESS' if all(r['all_gates_pass'] for r in records) else 'FAIL: do not build or ship a package',records=records,groups_checked=len(records),full900=full and len(records)==900,v5b_certificate='INAPPLICABLE'))
   print('RENDERER GATE',c,i,rec['all_gates_pass'],'JK ratio',rec['jk_ratio'],flush=True)
   # One failed gate already prohibits packaging. Keep the failed inputs for review.
   if not rec['all_gates_pass']:
    np.savez_compressed(REPO/'model_billmB/failed_renderer_probe.npz',original=original,candidate=x,genes=genes,target=targets[i],depth=d)
    if writer is not None:writer.close()
    raise SystemExit('FAILED renderer feasibility gate: no package will be created')
 if writer is not None:writer.close()
 job.save(REPO/('model_billmB/full_gates.json' if full else 'model_billmB/renderer_preflight.json'),dict(status='PASS all900' if full else 'PASS probes only; full900 still required',records=records,groups_checked=len(records),full900=full and len(records)==900,v5b_certificate='INAPPLICABLE',cache_directory=str(cache.relative_to(R)),cache_inputs_sha256=input_sha))
 print('PASS all900' if full else 'PASS all six probes; full900 still required',flush=True)
if __name__=='__main__':main()
