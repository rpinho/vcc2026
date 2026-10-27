"""Rebuild missing sufficient bulk totals, without integer cells or v9 layouts.

Existing v5b cached bulks are the exact all-target verification fixture. Frozen
alpha is read from the original record; no re-estimation or tuning takes place.
"""
import runtime_local as j
j.guard(launch=True)
import sys, types, importlib.util, gc
import numpy as np, anndata as ad, pandas as pd
from scipy import sparse

def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m)
    return m

def totals(q,b,depth):
    # Exact fitting prefix and column-total rule of frozen mean_counts and
    # repair_bulk_totals. Integer cell allocation cannot change these totals.
    x=q.astype(float).copy();x=.999*x+.001*b;desired=depth.sum()*b
    for it in range(100):
        col=depth @ x
        x*=np.divide(desired,col,out=np.zeros_like(col),where=col>0)
        x/=np.maximum(x.sum(1,keepdims=True),1e-30)
        if it%5==4 and np.abs(depth @ x/depth.sum()-b).sum()<2e-4:break
    expected=x*depth[:,None]
    columns=expected.sum(0)
    target=np.floor(columns).astype(np.int64)
    remaining=int(depth.sum()-target.sum())
    if remaining:target[np.argsort(columns-target,kind='stable')[-remaining:]]+=1
    assert target.sum()==depth.sum()
    return np.log1p(50000*target/target.sum())

if __name__=='__main__':
    rsi=j.VCC/'rsi';j.os.chdir(rsi)
    for p in [rsi,rsi/'vendor',rsi/'vendor/deps']:sys.path.insert(0,str(p))
    rt=types.ModuleType('runtime');rt.ROOT=rsi;rt.VCC=j.VCC;rt.json=j.json;rt.sha=j.sha;rt.guard=j.guard
    def forbidden(*a,**kw):raise RuntimeError('External write prohibited')
    rt.save=forbidden;rt.status=forbidden;sys.modules['runtime']=rt
    v4=load('v4_model',rsi/'vendor/v4_model.py')
    recipes=types.ModuleType('recipes');recipes.v4=v4;sys.modules['recipes']=recipes
    baseline=load('frozen_baseline',rsi/'v5b_loco_baseline.py')
    spec=j.json.loads((rsi/'artifacts/v5b_loco_spec.json').read_text())
    a=ad.read_h5ad(rsi/'artifacts/tasks/H1_rich_search.h5ad',backed='r')
    mask=a.obs.evaluation_role.astype(str).to_numpy()=='input_control'
    ctrl=a[np.flatnonzero(mask)].to_memory();a.file.close()
    ts=np.array(spec['tasks']['H1_rich_search']['targets']);genes=ctrl.var_names.to_numpy(str)
    j.status('Reconstructing frozen v5b source profile; alpha taken from existing score, no claim probe.')
    q,depth,L,Lc,_,_=baseline.profile(ctrl,'H1_rich_search',spec)
    original=j.json.loads((j.VCC/'final_push/artifacts/de_check/H1_rich_search/v5b/score.json').read_text())
    bf=baseline.helpers.from_logbulk(baseline.helpers.rescale(L,original['alpha_v5b'],Lc))
    cached=np.load(j.ROOT/'cached_v5b.npz');ci={t:i for i,t in enumerate(cached['targets'])}
    err=0
    for i,t in enumerate(ts):
        j.guard();m=totals(q,bf[i],depth)
        err=max(err,float(abs(m-cached['means'][ci[t]]).max()))
    assert err<2e-15,err
    j.status(f'Exact bulk-total reconstruction verified on all 79 cached v5b groups: max error {err:.3g}.')
    amp=baseline.helpers.from_logbulk(baseline.helpers.rescale(baseline.helpers.logbulk(bf),1.5,Lc))
    arrays={}
    arrays['v8b_like']=amp
    # Load source-only learned response with a safe stub instead of executing
    # model_billm/common.py (which changes threads and creates external caches).
    import torch
    torch.set_num_threads(1);torch.set_num_interop_threads(1)
    common=types.ModuleType('common');common.ROOT=j.VCC/'model_billm';sys.modules['common']=common
    dataset=load('frozen_billm_dataset',common.ROOT/'dataset.py')
    network=load('frozen_billm_network',common.ROOT/'network.py')
    src=dataset.Sources(['H1']);lookup=pd.Index(src.genes).get_indexer(genes);assert (lookup>=0).all()
    bulk=np.asarray(ctrl.X.sum(0,dtype=np.float64)).ravel();bulk/=bulk.sum()
    canonical=np.zeros(src.g);canonical[lookup]=bulk
    features=np.stack([src.features(t,canonical) for t in ts]);outputs=[]
    for seed in [1019,2027,4093]:
        j.guard();net=network.ResponseMLP(src.g).cpu()
        net.load_state_dict(torch.load(common.ROOT/'checkpoints/H1'/f'seed{seed}_weights.pt',map_location='cpu',weights_only=True));net.eval()
        with torch.no_grad():outputs.append(net(torch.from_numpy(features)).numpy())
        del net;gc.collect()
    learned=torch.stack([torch.from_numpy(o) for o in outputs]).mean(0).numpy()[:,lookup]
    arrays['levers_knob_none']=baseline.helpers.from_logbulk(Lc+network.blend(baseline.helpers.logbulk(bf)-Lc,learned,.3))
    del src,features,outputs,learned;gc.collect()
    for name,p in arrays.items():
        means=[]
        for i,t in enumerate(ts):
            j.guard();means.append(totals(q,p[i],depth))
            if i%20==0:j.status(f'{name}: reconstructing bulk totals {i+1}/79')
        # These are sufficient statistics, not per-cell render arms/packages.
        np.savez_compressed(j.ROOT/f'bulk_{name}.npz',targets=ts,means=means)
    j.save('bulk_reconstruction.json',dict(method='Frozen mean_counts fit + repair_bulk_totals target sums; integer/v9 allocation skipped',
        v5b_all79_max_abs_error=err,alpha=original['alpha_v5b'],no_calibration_refit=True,
        reason='No corresponding moments in 9551 retained cache files; original analogs retain tidy scores and B/claim only.',
        kernels={str(p):j.sha(p) for p in [rsi/'vendor/v4_model.py',rsi/'v5b_loco_baseline.py',common.ROOT/'dataset.py',common.ROOT/'network.py']}))
    j.status('Missing bulk sufficient statistics reconstructed; no cell layouts rerendered.')
