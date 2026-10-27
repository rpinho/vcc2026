"""Frozen B0 inference, H1 entirely excluded; labels restricted to training168."""
import runtime as rt
rt.enforce();rt.guard(True)
import numpy as np,pandas as pd,torch,gc,types,time
from torch import nn
from collections import Counter
torch.set_num_threads(1);torch.set_num_interop_threads(1)
start=time.monotonic();pins={}
def bind(p):
    p=rt.Path(p);pins[str(p)]=rt.sha(p);rt.guard();return p
ns=dict(np=np,pd=pd,torch=torch,nn=nn,Counter=Counter,gc=gc,j=rt)
rt.definitions(bind(rt.REPO/'gate_1025/streamed_inference.py'),['features'],ns)
rt.definitions(bind(rt.REPO/'model_billmB/network.py'),['ResponseMLP','blend','objective'],ns)
with np.load(bind(rt.REPO/'manifold_1026/inputs/data/H1_full.npz')) as z:
    targets=z['targets'];canonical=z['genes'];mask=z['measured'];truth=z['response'][:,mask].astype(float)
    training_control=z['control'];nc=z['n_cells']
assert len(targets)==168 and len(set(targets))==168
with np.load(bind(rt.REPO/'gatecall_1026/inputs/prediction/inference_H1.npz')) as z:
    # Predictions/control/axis only, never held-out reference outcomes.
    check_targets=z['targets'];genes=z['genes'];bulk=z['bulk'];check_response=z['response']
assert not set(targets)&set(check_targets)
ix=pd.Index(canonical).get_indexer(genes);assert np.array_equal(canonical[mask],genes)
control=np.zeros(len(canonical));control[ix]=bulk
rt.status('Reconstructing exact frozen features on 168 training targets; H1 source family excluded.')
features,native=ns['features']({'H1'},targets,control,bind);assert np.array_equal(native,canonical)
check_features,_=ns['features']({'H1'},check_targets,control,bind)
def predict(x):
    outputs=[]
    for lo in range(0,len(x),32):
        with torch.no_grad():learned=torch.stack([net(torch.from_numpy(x[lo:lo+32])) for net in models]).mean(0).numpy()
        outputs.append(ns['blend'](x[lo:lo+32,:len(canonical)],learned,.3))
    return np.concatenate(outputs)
models=[]
for seed in [1019,2027,4093]:
    meta=rt.json.loads(bind(rt.REPO/f'model_billmB/checkpoints/H1/seed{seed}.json').read_text());assert meta['excluded_families']==['H1'] and meta['steps']==80000
    net=ns['ResponseMLP'](len(canonical)).cpu();net.load_state_dict(torch.load(bind(rt.REPO/f'model_billmB/checkpoints/H1/seed{seed}_weights.pt'),map_location='cpu',weights_only=True));net.eval();models.append(net)
recheck=predict(check_features)[:,ix]
assert np.array_equal(recheck,check_response),np.max(abs(recheck-check_response))
raw=predict(features)[:,ix]
del models,net,recheck,check_features,check_response;gc.collect()
rt.status('PASS: reconstructed frozen inference is bitwise identical on all 79 existing predictions (no held-out truth accessed).')
# Exact continuous pseudobulk mapping used immediately before integer emitter.
# Layout preserves gene totals; integer rounding is explicitly outside this diagnostic.
def realized(r):
    probability=np.maximum(np.expm1(np.log1p(50000*bulk)+r),0);probability/=probability.sum(1,keepdims=True)
    return np.log1p(50000*probability)-np.log1p(50000*bulk)
pooled=features[:,:len(canonical)][:,ix];weighted=features[:,len(canonical):2*len(canonical)][:,ix]
family_sum={};family_den={};coverage=[[] for _ in targets];kd=[[] for _ in targets];source_stats=[]
manifest=rt.json.loads((rt.REPO/'model_billmB/data/manifest.json').read_text())
for name,m in sorted(manifest.items()):
    if m['family']=='H1':continue
    with np.load(rt.REPO/f'model_billmB/data/{name}.npz') as z:
        ts=z['targets'];ind=pd.Index(ts).get_indexer(targets);good=np.flatnonzero(ind>=0);meas=z['measured'];r=z['response'];ctrl=z['control'];n=z['n_cells']
        f=m['family'];family_sum.setdefault(f,np.zeros((168,len(genes)),np.float32));family_den.setdefault(f,np.zeros((168,len(genes)),np.float32))
        for i in good:
            row=ind[i];family_sum[f][i]+=r[row,ix]*meas[ix];family_den[f][i]+=meas[ix]
            coverage[i].append(dict(bank=name,family=f,n_cells=int(n[row])))
            own=pd.Index(canonical).get_indexer([targets[i]])[0]
            if own>=0 and meas[own]:
                cc=float(ctrl[row,own]);pert=max(float(np.expm1(np.log1p(50000*cc)+r[row,own])/50000),0)
                if cc>1e-7:kd[i].append(float(1-pert/cc))
        source_stats.append(dict(bank=name,family=f,total_rows=len(ts),training168_direct=int(len(good))))
    gc.collect()
family_predictions={f:realized(a/np.maximum(family_den[f],1e-8)) for f,a in family_sum.items()}
for f,a in family_predictions.items():
    for i,t in enumerate(targets):
        own=np.flatnonzero(genes==t)
        if len(own):a[i,own[0]]=0
    a[np.array([not any(v['family']==f for v in row) for row in coverage])]=np.nan
essential_path=rt.REPO/'research/sources/forrestsheldon__virtual-cell-challenge-2026/metadata/datasets/source_targets/replogle_k562_essential_counts.csv'
essential_genes=set(pd.read_csv(bind(essential_path)).source_target_gene)
essential=np.array([t in essential_genes for t in targets])
source_kd=np.array([np.median(a) if a else 0 for a in kd])
own_idx=pd.Index(genes).get_indexer(targets);truth_kd=[]
for i,j in enumerate(own_idx):
    cc=float(training_control[i,ix[j]]) if j>=0 else 0
    pp=max(float(np.expm1(np.log1p(50000*cc)+truth[i,j])/50000),0) if j>=0 else 0
    truth_kd.append(1-pp/cc if cc>1e-7 else None)
rt.save('results/training_metadata.json',dict(targets=targets.tolist(),coverage=coverage,source_knockdown=source_kd.tolist(),observed_H1_knockdown=truth_kd,essential_proxy=essential.tolist(),essential_definition='membership in Replogle K562 essential screen, NOT independently assayed H1 essentiality',sources=source_stats))
rt.guard()
np.savez_compressed(rt.ROOT/'results/baseline.npz',targets=targets,genes=genes,truth=truth.astype(np.float32),raw=raw,realized=realized(raw).astype(np.float32),pool=realized(pooled).astype(np.float32),weighted=realized(weighted).astype(np.float32),bulk=bulk,n_cells=nc,essential=essential,source_kd=source_kd,own_idx=own_idx,**{'family_'+f:a.astype(np.float32) for f,a in family_predictions.items()})
# Small realistic batch for CPU timing; temporary, deleted after timing.
by=np.zeros((64,len(canonical)),np.float32);by[:,ix]=truth[:64].astype(np.float32)
np.savez_compressed(rt.ROOT/'scratch/benchmark_batch.npz',x=features[:64],y=by,mask=np.tile(mask.astype(np.float32),(64,1)))
rt.save('results/provenance.json',dict(status='PASS bitwise frozen inference',pins=pins,training_targets=targets.tolist(),heldout_targets_metadata_only=check_targets.tolist(),heldout_outcomes_loaded=False,H1_labels_used_by_frozen_weights=False,LOO_explanation='Frozen outer-H1-excluded predictor is label-independent for every H1 target, hence is identical in all 168 target-LOO folds. Any new calibration must exclude row i.',continuous_pseudobulk=True,integer_renderer_not_run=True,seconds=time.monotonic()-start,nice=os_nice if (os_nice:=__import__('os').getpriority(__import__('os').PRIO_PROCESS,0)) else 0,free_gib=rt.guard()))
rt.status(f'Baseline and source-family diagnostics saved; {time.monotonic()-start:.1f} seconds. No fitting or held-out outcomes.')
