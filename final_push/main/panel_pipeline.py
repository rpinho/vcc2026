"""Panel/context arguments, existing Bill M blend + v5b amplitude + unchanged v9 layout/audit."""
import common as job
import io_guard;io_guard.install()
import sys,os,argparse,time,gc,hashlib
from pathlib import Path
from source_runtime import install,guard
install()
sys.path.insert(0,str(job.REPO/'model_v6'));sys.path.insert(1,str(job.REPO/'model_v4'))
import numpy as np,pandas as pd,torch
from scipy import sparse
import model,v5
from pipeline import read_genes,read_targets,load_controls,Context,Writer
from verify_counts import audit_h5
sys.path.insert(0,str(job.REPO/'model_v9'))
from de import prepare,Reference
from layout import emit
from audit_worker import audit
from infer import Ensemble
from network import blend
from fractions import Fraction as F
torch.set_num_threads(2);torch.set_num_interop_threads(1)

def generic_budget():
    # Reuse public source projections, collapse all context labels. No destination identity or panel fitting.
    r=job.json.loads((job.REPO/'model_v8/artifacts/public_scalars.json').read_text());values=[]
    for src in r['sources']:
        responders=[v for v in src['records'] if v['calls_observed']>=5]
        if len(responders)<8:responders=[v for v in src['records'] if v['calls_observed']>0]
        assert len(responders)>=4
        values.append(float(np.median([v['calls400_estimate'] for v in responders])))
    budget=int(round(np.exp(np.mean(np.log(np.maximum(values,1))))))
    return budget,dict(rule='geometric mean of public source responder-conditioned 400-cell call medians; independent of destination labels/target panel',
        source_medians=values,budget=budget,input_sha256=job.sha(job.REPO/'model_v8/artifacts/public_scalars.json'))

def apply_direction_knobs(a,ctx,targets,genes,ctrl):
    """Selector-best direction knobs (final_push/FINAL_PUSH.md 10-21): expression-compatibility mask, confidence-weighted
    per-target amplitude (mean 1, so the context budget solve is unchanged in level), floored control-CPM silence gate.
    All label-free: the context's own controls plus the public source banks. Applied to the deflated, gated log-bulk
    response ctx.L about ctx.Lc BEFORE the B = claim amplitude solve."""
    rec=dict(compat_cpm=a.compat_cpm,conf_alloc=bool(a.conf_alloc),silence_floor=a.silence_floor)
    if a.compat_cpm is None and not a.conf_alloc and a.silence_floor is None:return rec
    bulk=np.asarray(ctrl.X.sum(0,dtype=np.float64)).ravel();bulk/=bulk.sum()
    delta=ctx.L-ctx.Lc
    if a.compat_cpm is not None:
        keep=bulk*1e6>=a.compat_cpm;delta=delta*keep[None,:];rec['compat_genes_kept']=int(keep.sum())
    factor=np.ones(len(targets))
    if a.conf_alloc:
        sources=[model.aligned_bank(n,targets,genes,'log2fc') for n in model.load_banks(None,True)]
        cd=model.cd4_family(targets,genes,'log2fc')
        if cd is not None:sources.append(cd)
        sign=np.sign(delta);net=np.zeros(delta.shape)
        for e,m,v in sources:
            z=np.where(m&np.isfinite(v),e/np.sqrt(np.maximum(v,1e-12)),0);sig=np.abs(z)>1.96
            net+=(sig*(np.sign(z)==sign)).astype(float)-(sig*(np.sign(z)==-sign)).astype(float)
        support=(net>=1).sum(1).astype(float);conf=support/max(support.max(),1.0)
        f=np.maximum(conf,1e-3);f=f/max(f.mean(),1e-12);factor*=f
        rec.update(conf_factor_min=float(f.min()),conf_factor_max=float(f.max()),support_mean=float(support.mean()))
        del sources,cd
    if a.silence_floor is not None:
        gi={g:i for i,g in enumerate(genes)};lo,hi=1.0,10.0;fs=np.ones(len(targets))
        for i,t in enumerate(targets):
            j=gi.get(t)
            if j is not None:
                cpm=float(bulk[j]*1e6);fs[i]=float(np.clip((np.log(max(cpm,1e-9))-np.log(lo))/(np.log(hi)-np.log(lo)),0,1))
        fs=a.silence_floor+(1-a.silence_floor)*fs;factor*=fs
        rec.update(silence_factor_mean=float(fs.mean()),silence_below_one=int((fs<1).sum()))
    ctx.L=ctx.Lc+delta*factor[:,None];rec['row_factor_mean']=float(factor.mean())
    return rec

def context_cache(a,c,ci,targets,genes,ctrlpath,out):
    path=out/f'profile_{ci}.npz';done=out/f'profile_{ci}.json'
    if done.exists() and path.exists():
        r=job.json.loads(done.read_text());assert job.sha(path)==r['sha256'];return path
    ctrl=load_controls(ctrlpath,genes)
    # Full native control template, source fusion, PC deflation, call gate; no coverage-fill change to billmA arm.
    ctx=Context('v5b',targets,genes,ctrl,c,ci,workers=1)
    knobs=apply_direction_knobs(a,ctx,targets,genes,ctrl)
    probe=np.sort(np.random.default_rng(ci).choice(len(targets),min(40,len(targets)),replace=False));jks=[]
    for i in probe:
        guard(floor=25)
        x=model.emit(ctx.q,None,ctx.b[i],ctx.d,ctx.seed(targets[i]),'anchor_single_mean').astype(np.int32)
        jks.append(float(v5.jackknife_correction(sparse.csr_matrix(x),np.zeros(400,np.intp),1,50000.,chunk=32)[0]))
    weights=v5._row_weights(len(targets),len(genes),v5.excl_cols(targets,genes))
    claim=float(np.mean(jks)*weights.sum());alpha,B=v5.solve_alpha(ctx.L,targets,genes,claim,center=ctx.Lc)
    bf=v5.from_logbulk(v5.rescale(ctx.L,alpha,ctx.Lc))
    ensemble=Ensemble('production',data_dir=a.response_banks,weights_dir=a.weights)
    canonical=ensemble.src.genes;ix=pd.Index(canonical).get_indexer(genes);assert (ix>=0).all()
    bulk=np.asarray(ctrl.X.sum(0,dtype=np.float64)).ravel();bulk/=bulk.sum()
    controls=np.zeros(len(canonical));controls[ix]=bulk
    learned=[]
    for lo in range(0,len(targets),32):learned.append(ensemble.predict(targets[lo:lo+32],controls)['learned'][:,ix])
    learned=np.concatenate(learned);new=v5.from_logbulk(ctx.Lc+blend(v5.logbulk(bf)-ctx.Lc,learned,.3))
    del ensemble;gc.collect()
    # Same v9 confidence calculation as billmA, using argument-selected banks.
    net=np.zeros(ctx.L.shape);zs=np.zeros(ctx.L.shape);sign=np.sign(ctx.L-ctx.Lc)
    sources=[model.aligned_bank(n,targets,genes,'log2fc') for n in model.load_banks(None,True)]
    cd=model.cd4_family(targets,genes,'log2fc')
    if cd is not None:sources.append(cd)
    for e,m,v in sources:
        z=np.where(m&np.isfinite(v),e/np.sqrt(np.maximum(v,1e-12)),0);sig=np.abs(z)>1.96
        net+=(sig*(np.sign(z)==sign)).astype(float)-(sig*(np.sign(z)==-sign)).astype(float)
        zs+=np.minimum(np.abs(z),100)*(np.sign(z)==sign)
    refdir=job.ROOT/f'artifacts/ref_{a.name}_{ci}'
    if not ((refdir/'complete.json').exists() and (refdir/'cdf_scale.npy').exists() and (refdir/'values.npy').exists()):
        prepare(ctrl.X,refdir,np.random.default_rng(20261018+ci))
        job.save(refdir/'complete.json',dict(control_sha256=job.sha(ctrlpath),genes=len(genes)))
    np.savez_compressed(path,q=ctx.q,d=ctx.d,new=new,confidence=net*1e6+zs,genes=genes,targets=targets)
    job.save(done,dict(context=c,alpha=alpha,B=B,claim=claim,sha256=job.sha(path),profile_info=ctx.info,direction_knobs=knobs,
        adapter='panel-specific label-free alpha and context-free public DE call budget; v9 kernels unchanged'))
    del ctrl,ctx,sources,cd;gc.collect();return path

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--targets',required=True);ap.add_argument('--contexts',required=True)
    ap.add_argument('--genes',required=True);ap.add_argument('--banks',required=True);ap.add_argument('--response-banks',required=True)
    ap.add_argument('--weights',required=True);ap.add_argument('--name',required=True);ap.add_argument('--probes',action='store_true')
    ap.add_argument('--compat-cpm',type=float,default=None,help='zero the bulk response on genes under this control CPM')
    ap.add_argument('--conf-alloc',action='store_true',help='per-target amplitude x cross-source confidence (mean 1)')
    ap.add_argument('--silence-floor',type=float,default=None,help='floored control-CPM silence gate on the target gene (floor value)')
    a=ap.parse_args();assert '/' not in a.name and '..' not in a.name
    start=time.time();contexts=dict(kv.split('=',1) for kv in a.contexts.split(','));labels=list(contexts)
    genes=read_genes(a.genes);panels=read_targets(a.targets,labels)
    assert len(set(genes))==len(genes)
    for c,tt in panels.items():assert len(tt)>0 and len(set(tt))==len(tt),(c,'invalid targets')
    out=job.ROOT/f'artifacts/panels/{a.name}';out.mkdir(parents=True,exist_ok=True)
    signature=dict(targets=job.sha(a.targets),genes=job.sha(a.genes),contexts={c:job.sha(p) for c,p in contexts.items()},
        banks={p.name:job.sha(p) for p in Path(a.banks).glob('*.npz')},response_banks={p.name:job.sha(p) for p in Path(a.response_banks).glob('*.npz')},
        weights={p.name:job.sha(p) for p in (Path(a.weights)/'production').glob('*_weights.pt')},
        kernels={f:job.sha(job.REPO/f) for f in ['model_v6/pipeline.py','model_v6/model.py','model_v5/v5.py','model_v9/layout.py','model_v9/audit_worker.py']},
        job_sources={n:job.sha(job.ROOT/n) for n in ['common.py','dataset.py','infer.py','network.py','source_runtime.py']},
        public_call_prior=job.sha(job.REPO/'model_v8/artifacts/public_scalars.json'),
        direction_knobs=dict(compat_cpm=a.compat_cpm,conf_alloc=bool(a.conf_alloc),silence_floor=a.silence_floor),
        builder=job.sha(__file__))
    assert len(signature['weights'])==3
    bound=out/'inputs.json'
    if bound.exists():assert job.json.loads(bound.read_text())==signature,'Inputs changed; choose a fresh name'
    else:job.save(bound,signature)
    model.configure(a.banks,(),Path(a.banks)/'promoter_pairs.csv')
    nconf,policy=generic_budget();job.save(out/'budget_policy.json',policy)
    records=[]
    if not a.probes:
        pre=job.json.loads((out/'probes.json').read_text());assert pre['status']=='PASS' and len(pre['records'])==sum(min(2,len(tt)) for tt in panels.values())
    for ci,c in enumerate(labels):
        tt=panels[c];cache=context_cache(a,c,ci,tt,genes,contexts[c],out)
        refname=f'{a.name}_{ci}';ref=Reference(refname)
        with np.load(cache) as z:q=z['q'];d=z['d'];new=z['new'];confidence=z['confidence']
        probe_ids=sorted(set([0,min(4,len(tt)-1)]))
        for i in (probe_ids if a.probes else range(len(tt))):
            guard(floor=25);receipt=out/f'group_{ci}_{i:03d}.json';counts=out/f'group_{ci}_{i:03d}.npz'
            if receipt.exists() and counts.exists():
                r=job.json.loads(receipt.read_text());assert r['all_gates_pass'] and r['counts_sha256']==job.sha(counts)
                records.append(r);continue
            seed=model.atlas.seed_for(c+':'+tt[i],20260910)
            original=model.emit(q,None,new[i],d,seed,'anchor_single_mean').astype(np.int32)
            own_columns=np.flatnonzero(genes==tt[i])
            own=int(own_columns[0]) if len(own_columns) else np.empty(0,dtype=np.intp)
            x,called,info=emit(original,ref,confidence[i],nconf,own)
            # audit's context key is a filesystem reference name, independent of destination label.
            rec,_=audit(x,original,refname,confidence[i],own,called)
            rec['own_target_on_gene_axis']=bool(len(own_columns))
            if not np.isfinite(rec['confidence_rho']):
                rec['confidence_rho']=None;rec['confidence_rho_undefined']='constant or insufficient confidence ranks; diagnostic only'
            rec.update(context=c,index=i,target=str(tt[i]),JK_dominance=F(rec['jk'])>=F(rec['v5b_jk']),
                baseline_counts_sha256=hashlib.sha256(original.astype('<i4').tobytes()).hexdigest(),
                dense_counts_sha256=rec['counts_sha256'],v5b_preservation='INAPPLICABLE: same learned profile only',**info)
            rec['all_gates_pass']=bool(rec['JK_dominance'] and np.nextafter(rec['mean_certificate']['mean_L1_upper'],np.inf)<.015 and np.nextafter(rec['mean_certificate']['bulk_L1_upper'],np.inf)<.0015)
            if not rec['all_gates_pass']:
                job.save(out/f'failed_{ci}_{i:03d}.json',rec)
                np.savez_compressed(out/f'failed_{ci}_{i:03d}.npz',original=original,candidate=x,genes=genes,target=tt[i])
                raise RuntimeError(f'Fail closed renderer gate {c}/{tt[i]}')
            sparse.save_npz(counts,sparse.csr_matrix(x));rec['counts_sha256']=job.sha(counts)
            job.save(receipt,rec);records.append(rec)
            print('PANEL GATE',c,i,tt[i],rec['jk_ratio'],'elapsed',round(time.time()-start),flush=True)
        del q,d,new,confidence;gc.collect()
    job.save(out/('probes.json' if a.probes else 'full_gates.json'),dict(status='PASS',records=records,inputs_sha256=job.sha(bound),seconds=time.time()-start))
    if a.probes:return
    dest=job.ROOT/f'out/{a.name}.h5ad';dest.parent.mkdir(parents=True,exist_ok=True);writer=Writer(dest,genes)
    try:
        for ci,c in enumerate(labels):
            for i,t in enumerate(panels[c]):writer.append(sparse.load_npz(out/f'group_{ci}_{i:03d}.npz'),c,t)
    finally:writer.close()
    structure=audit_h5(dest,genes,panels,labels)
    job.save(out/'build.json',dict(status='PASS all groups, independent integer/CSR/axes audit',structure=structure,
        seconds=time.time()-start,groups=len(records),input_sha256=job.sha(bound)))
    print('PANEL BUILD PASS',structure,flush=True)
if __name__=='__main__':main()
