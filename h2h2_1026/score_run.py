import local_runtime as j
j.enforce_io();j.guard()
import os,sys,gc,hashlib,argparse
import numpy as np,pandas as pd,polars as pl,anndata as ad
from scipy import sparse
from scipy.special import ndtri
import score_official as score
import memory_policy
TASKS={'H1':'H1_rich_search'}
NATIVE_MATRIX_ID=score.matrix_id
# Only named immutable functions/classes are compiled; predecessor top-level is never run.
names=['Backed','matrix_id','load_task','modules','reference','conf_array','install_scorer','setup','arrsha']
j.definitions('gatecall_experiment.py',names,globals());score.matrix_id=matrix_id

def status(s):
 (j.ROOT/'STATUS.md').write_text('# STATUS\n\n'+j.utc()+' '+s+'\n\nSix response matrices retained; no tuning/package/network.\n');print(j.utc(),s,flush=True)

def lock():
 assert j.sha(j.ROOT/'PREREGISTRATION.md')==(j.ROOT/'results/preregistration.sha256').read_text().split()[0]
 return j.sha(j.ROOT/'PREREGISTRATION.md')

def confidence_all(model,targets,genes):
 path=j.ROOT/'scratch/confidences.npz'
 if path.exists():return
 responses={a:np.load(j.ROOT/f'results/predictions_{a}.npz')['response'] for a in ['B0','B1','B2','B3','B4','B5']}
 net={a:np.zeros(r.shape) for a,r in responses.items()};zs={a:np.zeros(r.shape) for a,r in responses.items()}
 def add(e,m,v):
  z=np.where(m&np.isfinite(v),e/np.sqrt(np.maximum(v,1e-12)),0);sig=abs(z)>1.96
  for a,r in responses.items():
   same=np.sign(z)==np.sign(r);opp=np.sign(z)==-np.sign(r)
   net[a]+=(sig*same).astype(float)-(sig*opp).astype(float);zs[a]+=np.minimum(abs(z),100)*same
 model.configure(j.ROOT/'inputs/banks',exclude=['H1'])
 for name in model.load_banks(None,True):
  status('Frozen sign-dependent confidence: '+name);e,m,v=model.aligned_bank(name,np.array(targets),genes,'log2fc');add(e,m,v);del e,m,v;gc.collect()
 cd=model.cd4_family(np.array(targets),genes,'log2fc')
 if cd is not None:add(*cd)
 del cd;gc.collect()
 conf={a:net[a]*1e6+zs[a] for a in responses}
 expected=j.json.loads((j.ROOT/'inputs/baseline/H1.json').read_text())['renderer']['confidence_sha256'];assert arrsha(conf['B0'])==expected
 np.savez_compressed(path,**conf)

def run(arms):
 lock();status('Loading frozen control/reference identities and renderer')
 data=setup('H1');a,ctrl,real,base,anchors,targets,genes,response,bulk,conf,q,depth,ref,model,layout=data
 confidence_all(model,targets,genes)
 bounded=install_scorer()
 # Capture the official DE tables; permits independent target-metric audits.
 import cell_eval2.run as cr
 native_prepare=cr._prepare_de_cached
 current=['']
 def capture(pred,truth,*aa,**kk):
  pred.write_parquet(j.ROOT/f'artifacts/de_pred_{current[0]}.parquet')
  p=j.ROOT/'artifacts/de_real.parquet'
  if not p.exists():truth.write_parquet(p)
  return native_prepare(pred,truth,*aa,**kk)
 cr._prepare_de_cached=capture
 for arm in arms:
  dest=j.ROOT/f'results/official_{arm}.json'
  if dest.exists():print('Retaining completed',arm,flush=True);continue
  current[0]=arm;start=j.utc();lock();j.guard();records=[];mom=[]
  with np.load(j.ROOT/f'results/predictions_{arm}.npz') as z:
   assert z['targets'].tolist()==targets and np.array_equal(z['genes'],genes);r=z['response']
  with np.load(j.ROOT/'scratch/confidences.npz') as z:confidence=z[arm]
  def predictions():
   for i,t in enumerate(targets):
    j.guard();status(f'{arm}: render + official score {i+1}/79 {t}')
    own=np.flatnonzero(genes==t);own=int(own[0]) if len(own) else np.empty(0,dtype=np.intp)
    probability=np.maximum(np.expm1(np.log1p(50000*bulk)+r[i]),0);probability/=probability.sum()
    original=model.emit(q,None,probability,depth,model.atlas.seed_for('H1:'+t,20260910),'anchor_single_mean').astype(np.int32)
    x,called,detail=layout.emit(original,ref,confidence[i],int(base['renderer']['budget']),own)
    assert np.array_equal(x.sum(0,dtype=np.int64),original.sum(0,dtype=np.int64))
    assert np.array_equal(x.sum(1,dtype=np.int64),depth) and np.all(x>=0)
    assert np.array_equal(x[:,~ref.keep],original[:,~ref.keep])
    if isinstance(own,int):assert np.array_equal(x[:,own],original[:,own])
    records.append(dict(target=t,selected=called.tolist(),selected_genes=genes[called].tolist(),detail=detail,counts_sha256=arrsha(x),original_sha256=arrsha(original)))
    # Small gene-level count moments, not a full per-cell cache.
    mom.append(x.sum(0,dtype=np.int64))
    yield t,sparse.csr_matrix(x)
    del x,original;gc.collect()
  frame,digest=bounded.score_iter(real,predictions(),score.config())
  raw=score.aggregate(frame);scaled,avg=score.scale(raw,anchors)
  if arm=='B0':
   assert digest==base['arms']['incumbent']['prediction_sha256'],('digest',digest)
   assert all(abs(raw[k]-base['arms']['incumbent']['raw'][k])<1e-12 for k in score.SIX),('raw',raw)
  frame.write_csv(j.ROOT/f'results/tidy_{arm}.csv')
  cols=sorted(set(frame['metric'].to_list()));mat=np.full((79,len(cols)),np.nan);ti={t:i for i,t in enumerate(targets)};mi={m:i for i,m in enumerate(cols)}
  for t,m,value in frame.iter_rows():
   if t in ti:mat[ti[t],mi[m]]=np.nan if value is None else value
  np.savez_compressed(j.ROOT/f'results/target_metrics_{arm}.npz',targets=np.array(targets),columns=np.array(cols),values=mat)
  np.savez_compressed(j.ROOT/f'results/rendered_bulk_{arm}.npz',targets=np.array(targets),genes=genes,count_sums=np.array(mom),depth=depth)
  j.save(f'results/render_{arm}.json',dict(records=records,budget=885,template_seed=20260910,confidence_sha256=arrsha(confidence)))
  j.save(f'results/official_{arm}.json',dict(arm=arm,targets=targets,raw=raw,scaled=scaled,score_avg=avg,prediction_sha256=digest,config=score.config().to_dict(),anchors=anchors,real_sha256=matrix_id(real),control_sha256=matrix_id(ctrl),preregistration_sha256=lock(),started=start,ended=j.utc(),pid=os.getpid(),nice=os.getpriority(os.PRIO_PROCESS,0)))
  # Preserve warm DE cache across all authorized sequential arms.
  status(f'{arm} complete: score_avg={avg:.9f}, fidelity={raw[score.SIX[3]]:.9f}, PDS={raw[score.SIX[1]]:.9f}')
 a.file.close()
if __name__=='__main__':
 assert os.getpriority(os.PRIO_PROCESS,0)>=10
 parser=argparse.ArgumentParser();parser.add_argument('arms',nargs='*',default=['B0','B1','B2','B3','B4','B5']);args=parser.parse_args()
 # One compute worker. A reserve-monitor thread only samples disk and releases heap.
 free=j.shutil.disk_usage(j.ROOT).free
 assert free>=31_000_000_000,("Launch requires 31 GB decimal",free)
 cache=list((j.ROOT/"artifacts/de_cache").glob("*.parquet"));assert cache,"Warm DE cache missing"
 pins=j.json.loads((j.ROOT/"results/input_pins.json").read_text())
 receipt=j.ROOT/"results/operational_resume.json"
 assert j.sha(receipt)==pins["operational_resume_receipt_sha256"]
 prior=j.json.loads(receipt.read_text())
 for arm,digest in prior["predictions"].items():
  path=j.ROOT/f"results/predictions_{arm}.npz"
  assert path.is_file(),f"Missing prediction: {path}"
  actual=j.sha(path);assert actual==digest,f"Prediction hash mismatch: {path}; expected {digest}; got {actual}"
 j.save(f"results/authorized_start_{os.getpid()}.json",dict(utc=j.utc(),pid=os.getpid(),nice=os.getpriority(os.PRIO_PROCESS,0),free_bytes=free,launch_headroom_bytes=31_000_000_000,proactive_stop_bytes=30_500_000_000,warm_cache_files=len(cache),arms=args.arms,preregistration_sha256=lock()))
 finish=memory_policy.install()
 try:run(args.arms)
 finally:finish()
