"""Bounded local ablations with official cell-eval2; no submission.
Run: ../.venv/bin/python model_de/harness.py --dest H1 --targets 12
Destination family is removed from EVERY source channel. Target sampling uses
cell availability + SHA256 ordering, never observed responses. All predictions
have 400 cells; real observations are never duplicated. Raw metrics are primary.
"""
import os,sys
if '--threads' in sys.argv:
    os.environ['VCC_DE_THREADS']=sys.argv[sys.argv.index('--threads')+1]
import runtime as rt
import sys, argparse, time, gc, hashlib, logging,resource
sys.path.extend([str(rt.VCC/'model_v6'),str(rt.VCC/'proxy/scorer')])
import numpy as np, pandas as pd, anndata as ad
from scipy import sparse
from dataclasses import replace
import model, v5, official_scoring as scoring, bounded_score, bounded_de
import cell_eval2.run as run
from cell_eval2 import aggregate_metrics
from strategy import rank_sources,select,selective_template,redistribute
logging.basicConfig(level=logging.ERROR)
LOADED_HARNESS_SHA=rt.sha(__file__)
LOADED_STRATEGY_SHA=rt.sha(rt.ROOT/'strategy.py')

def load_task(dest,n,mapped=False):
    path=rt.VCC/('rsi/artifacts/tasks/H1_rich_search.h5ad' if dest=='H1' else f'proxy/tasks/{dest}.h5ad')
    a=ad.read_h5ad(path,backed='r');labs=a.obs.target_gene.astype(str).to_numpy()
    if dest=='H1':
        role=a.obs.evaluation_role.astype(str).to_numpy()
        labs[np.isin(role,['input_control','reference_control'])]='non-targeting'
    counts=pd.Series(labs[labs!='non-targeting']).value_counts()
    panel=sorted(counts[counts>=40].index,key=lambda t:hashlib.sha256(('de1018:'+t).encode()).hexdigest())[:n]
    rows=list(np.flatnonzero(labs=='non-targeting'))
    if len(rows)>3200: rows=sorted(np.random.default_rng(1018).choice(rows,3200,replace=False))
    for t in panel:
        ix=np.flatnonzero(labs==t)
        if len(ix)>400: ix=np.sort(np.random.default_rng(1018).choice(ix,400,replace=False))
        rows.extend(ix)
    rows=np.sort(rows);mapping=None
    if mapped:
        folder=rt.VCC/f'proxy/v2_cache/mapped/{dest}'
        receipt_path=folder/'receipt.json';receipt_text=receipt_path.read_text();mapping=rt.json.loads(receipt_text)
        identity=mapping['source_identity'];st=path.stat()
        assert identity['size']==st.st_size and identity['mtime_ns']==st.st_mtime_ns,'mapped input changed'
        arrays={k:np.load(folder/(k+'.npy'),mmap_mode='r') for k in ('data','indices','indptr')}
        nnz=mapping['nnz'];full=sparse.csr_matrix((arrays['data'][:nnz],arrays['indices'][:nnz],arrays['indptr']),shape=mapping['shape'],copy=False)
        real=ad.AnnData(full[rows].copy(),obs=a.obs.iloc[rows].copy(),var=a.var.copy())
        assert receipt_path.read_text()==receipt_text,'mapped receipt changed during read'
        mapping=dict(mapping,source_identity=dict(identity,source=str(path.relative_to(rt.VCC))))
    else:real=a[rows].to_memory()
    real.obs['target_gene']=labs[rows];a.file.close()
    ctrl=real[real.obs.target_gene.to_numpy()=='non-targeting'].copy()
    return real,ctrl,np.array(panel),dict(path=str(path.relative_to(rt.VCC)),
        targets=panel,real_counts={t:int((real.obs.target_gene==t).sum()) for t in panel},
        controls=ctrl.n_obs,sha256=rt.sha(path),selection='>=40 cells, SHA256 de1018:target order, max400 without replacement',
        mapping=mapping,median_control_umi=float(np.median(np.asarray(ctrl.X.sum(1)).ravel())))

def verify_adapter():
    b=rt.json.loads((rt.VCC/'rsi/artifacts/bounded_fixture.json').read_text())
    s=rt.json.loads((rt.VCC/'rsi/artifacts/bounded_score_fixture.json').read_text())
    assert b['status']=='PASS bitwise equality' and s['status']=='PASS bitwise equality'
    assert rt.sha(bounded_de.__file__)==b['adapter_sha256']
    assert rt.sha(bounded_score.__file__)==s['driver_sha256']
    assert scoring.provenance()['official_files']==b['cell_eval2']['official_files']
    bounded_de.install_unchecked()

def main(args):
    rt.guard();verify_adapter()
    acceleration=None
    if args.fast:
        import fast_de
        acceleration=fast_de.install()
    start=time.time();real,ctrl,targets,task=load_task(args.dest,args.targets,args.mapped)
    print('TASK',args.dest,len(targets),'targets',ctrl.n_obs,'controls',flush=True)
    genes=real.var_names.to_numpy(str);model.configure(model.SHARED_BANKS,[args.dest],rt.VCC/'model_v6/artifacts/promoter_pairs_v4.csv')
    q4,d4,mean,bulk=model.template(ctrl,genes,pool_k=4)
    q1,d1,_,_=model.template(ctrl,genes,pool_k=1)
    (p,b),cov=model.heads(targets,genes,mean,bulk,'anchor_single_mean')
    print('HEADS',args.dest,'ready',round(time.time()-start,1),'seconds',flush=True)
    lfc=np.log2((b+1e-15)/(bulk[None,:]+1e-15));sources=[];names=model.load_banks(None,True)
    for name in names:
        rt.guard();sources.append(model.aligned_bank(name,targets,genes,'log2fc'))
    cd=model.cd4_family(targets,genes,'log2fc')
    if cd is not None: sources.append(cd);names.append('CD4')
    print('SOURCE RANKING',args.dest,names,round(time.time()-start,1),'seconds',flush=True)
    confidence,agreement,direction,budgets=rank_sources(sources,lfc,mean*1e6>5,genes,targets)
    del sources,cd;gc.collect()
    masks={str(k):select(confidence,agreement,k) for k in (100,300,1000)}
    masks['adaptive']=select(confidence,agreement,budgets)
    magnitude=np.abs(lfc).copy();magnitude[:,mean*1e6<=5]=-np.inf
    for i,t in enumerate(targets):
        own=np.flatnonzero(genes==t)
        if own.size:magnitude[i,own[0]]=-np.inf
    masks['magnitude']=select(magnitude,np.ones_like(agreement),300)
    # Freeze the source-only design before reading evaluation DE or scoring any arm.
    design=dict(task=task,source_families=names,destination_excluded=args.dest,
        source_budgets=budgets.tolist(),selected={k:m.sum(1).tolist() for k,m in masks.items()},
        caveat='small-panel public-line diagnostic; uncalibrated proxy; raw MSE budget is panel-dependent')
    label=args.dest+('_matched' if args.mapped else '')
    out=rt.ROOT/'artifacts'/label;out.mkdir(parents=True,exist_ok=True);rt.save(out/'design.json',design)
    arms=args.arms.split(',')
    for arm in arms:
        rt.guard();t0=time.time();diagnostics=[]
        original=run._materialize_de_sides
        def capture(*a,**kw):
            rd,pd_=original(*a,**kw)
            for target in pd_['target'].unique().to_list():
                r=rd.filter(rd['target']==target);p_=pd_.filter(pd_['target']==target)
                diagnostics.append(dict(target=target,real_calls=int((r['p_adj']<.05).sum()),
                    predicted_calls=int((p_['p_adj']<.05).sum())))
            return rd,pd_
        run._materialize_de_sides=capture
        invariants=[]
        def predictions():
            for i,t in enumerate(targets):
                print('EMIT',args.dest,arm,i+1,'/',len(targets),t,flush=True)
                rt.guard();seed=model.atlas.seed_for(args.dest+':'+t,20260910)
                q=q4;d=d4;bb=b[i];mask=masks['300'][i]
                if arm=='pool1': q=q1;d=d1
                elif arm.startswith('select'):
                    key=arm.removeprefix('select');mask=masks[key][i]
                    q=selective_template(q1,q4,mask);d=d4
                elif arm=='amplitude09':
                    effect=np.log2((bb+1e-15)/(bulk+1e-15))
                    bb=bulk*np.exp2(np.clip(effect*1.5,-6,6));bb/=bb.sum()
                elif arm in ('quant1000','quant1000_exchange'):
                    d=np.maximum(1000,np.rint(d4/1000)*1000).astype(np.int64)
                elif arm not in ('pool4','compress300','shaped300','shapedadaptive'): raise ValueError(arm)
                if arm.startswith('shaped'):
                    key=arm.removeprefix('shaped');mask=masks[key][i]
                    x0=np.rint(q1*d1[:,None]).astype(np.int32)
                    x=v5.shaped_emit(x0,bb,np.random.default_rng(seed),called=mask).astype(np.int32)
                else:
                    x=model.mean_counts(q,None,bb,d,seed)
                if arm=='compress300':
                    before=x.sum(0);x=redistribute(x,mask)
                    assert np.array_equal(before,x.sum(0))
                moment_certificate=None
                filtered_support=None
                if arm=='quant1000_exchange':
                    from exchange import exchange_equal_depth
                    from exchange_verify import certify
                    own=np.flatnonzero(genes==t)
                    observed=np.einsum('ij,i->j',x,1/np.asarray(x.sum(1),float),optimize=False)/x.shape[0]
                    mask=mask & (np.sign(observed-mean)==direction[i]);filtered_support=int(mask.sum())
                    y=exchange_equal_depth(x,mask,preferred_buffer=int(own[0]) if own.size else None)
                    moment_certificate=certify(x,y);x=y
                bounded_score.validate(x)
                invariants.append(dict(target=t,depth_min=int(x.sum(1).min()),depth_max=int(x.sum(1).max()),
                                       two_moment_certificate=moment_certificate,support_after_actual_sign_filter=filtered_support))
                yield str(t),sparse.csr_matrix(x)
        try:
            cfg=replace(scoring.config(),num_threads=args.threads,gather_threads=1)
            tidy,digest=bounded_score.score_iter(real,predictions(),cfg)
        finally: run._materialize_de_sides=original
        agg=aggregate_metrics(tidy);raw=dict(zip(agg['metric'].to_list(),agg['mean'].to_list()))
        rec=dict(arm=arm,dest=args.dest,raw=raw,prediction_sha256=digest,seconds=time.time()-t0,
            de_counts=diagnostics,validity=invariants,design='design.json',status='LOCAL DIAGNOSTIC; not official validation',
            harness_sha256=LOADED_HARNESS_SHA,strategy_sha256=LOADED_STRATEGY_SHA,de_threads=args.threads)
        rt.save(out/(arm+'.json'),rec);tidy.write_csv(out/(arm+'_per_target.csv'))
        print(args.dest,arm,rt.json.dumps(raw),round(rec['seconds'],1),flush=True);gc.collect()
    rt.save(out/'provenance.json',dict(scorer=scoring.provenance(),seconds=time.time()-start,
        harness_sha256=LOADED_HARNESS_SHA,strategy_sha256=LOADED_STRATEGY_SHA,
        acceleration=acceleration,fast_verification=rt.json.loads((rt.ROOT/'artifacts/verification.json').read_text()) if args.fast else None,
        effective_nice=os.getpriority(os.PRIO_PROCESS,0),
        peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--dest',default='H1');ap.add_argument('--targets',type=int,default=12)
    ap.add_argument('--threads',type=int,choices=(1,2),default=1)
    ap.add_argument('--mapped',action='store_true',help='Read existing feedback-job depth-matched mmap task, without building or editing it')
    ap.add_argument('--fast',action='store_true',help='Use independently bitwise-verified sorted-reference MWU arithmetic')
    ap.add_argument('--arms',default='pool4,pool1,select100,select300,select1000,selectadaptive,selectmagnitude,compress300,amplitude09,shaped300,shapedadaptive')
    main(ap.parse_args())
