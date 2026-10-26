"""Emit two fixed-total candidates and retain independent per-group gate receipts."""
import runtime,sys,os,gc,time
import numpy as np,pandas as pd
from concurrent.futures import ProcessPoolExecutor
from io_counts import Writer
from build_worker import group
from cell_eval2.metrics.delta import _across_pert_budget,_row_weights
sys.path.append(str(runtime.ROOT.parent/'validation'))
import gates_core
R=runtime.ROOT
def main():
 assert runtime.json.loads((R/'artifacts/recovered_v5b.json').read_text())['status'].startswith('PASS')
 (R/'artifacts/de').mkdir(exist_ok=True)
 z=np.load(R/'artifacts/stage1_v5b.npz');genes=z['genes'].astype(str);targets=z['targets'].astype(str);totals=z['totals']
 cols={i:int(np.flatnonzero(genes==t)[0]) for i,t in enumerate(targets)};weights=_row_weights(300,len(genes),cols)
 files={n:Writer(R/f'out/{n}.h5ad',genes) for n in ['v9a','v9b']};records={n:[] for n in files};start=time.time()
 hashes={p.name:runtime.sha(p) for p in R.glob('*.py')};runtime.save(R/'artifacts/build_manifest.json',dict(source_sha256=hashes,scorer_sha256=runtime.scorer_hashes(),baseline=runtime.sha(R/'artifacts/recovered_v5b.json'),v8a_budget_receipt=runtime.sha(R.parent/'model_v8/artifacts/candA_gates.json')))
 try:
  with ProcessPoolExecutor(max_workers=int(os.environ.get('V9_WORKERS','3'))) as pool:
   for ci,c in enumerate('ABC'):
    pending={};nxt=0
    for i,t in enumerate(targets):
     runtime.guard()
     while nxt<min(i+5,300):pending[nxt]=pool.submit(group,(c,ci,nxt));nxt+=1
     xa,xb,ra,rb=pending.pop(i).result()
     for name,x,rec in [('v9a',xa,ra),('v9b',xb,rb)]:
      assert np.array_equal(x.sum(0,dtype=np.int64),totals[ci,i]);files[name].append(x,c,t);records[name].append(rec)
     if i%15==0:
      print('BUILD',c,i,'seconds',round(time.time()-start),'calls',ra['n_calls'],rb['n_calls'],'JK ratios',round(ra['jk_ratio'],4),round(rb['jk_ratio'],4),flush=True)
      runtime.save(R/'artifacts/build_progress.json',dict(records=records,seconds=time.time()-start,pids=[p.pid for p in pool._processes.values()]))
 finally:
  for w in files.values():w.close()
 for name in files:
  assert len(records[name])==900
  summary=[]
  for ci,c in enumerate('ABC'):
   rows=[r for r in records[name] if r['context']==c];T=totals[ci];L=np.log1p(50000*T/T.sum(1,keepdims=True));B=float(_across_pert_budget(L,cols,weights));independent=gates_core.across_pert_budget(L,cols,gates_core.row_weights(300,len(genes),cols));assert np.isclose(B,independent,atol=1e-12)
   claim=float(weights @ np.array([r['jk'] for r in rows]));oldclaim=float(weights @ np.array([r['v5b_jk'] for r in rows]));calls=[r['n_calls'] for r in rows]
   summary.append(dict(context=c,B=B,claim=claim,v5b_claim=oldclaim,B_over_claim=B/claim,uncapped_usable_credit=min(B,claim),v5b_uncapped_usable_credit=min(B,oldclaim),calls_min_median_max=[min(calls),float(np.median(calls)),max(calls)],confidence_rho_min_median=[float(min(r['confidence_rho'] for r in rows)),float(np.median([r['confidence_rho'] for r in rows]))],max_mean_CPM_L1=max(r['mean_CPM_L1'] for r in rows),totals_depths_nonDE_exact=True))
  runtime.save(R/f'artifacts/{name}_gates.json',dict(status='COMPLETE independent prediction-side diagnostics; hidden scores unverified',contexts=summary,records=records[name]))
 runtime.save(R/'artifacts/build_complete.json',dict(shape=[360000,18533],groups_per_candidate=900,seconds=time.time()-start))
 print('BUILD COMPLETE',flush=True)
if __name__=='__main__':main()
