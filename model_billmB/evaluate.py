"""Official cell-eval2 on frozen RSI LOCO SEARCH tasks after GPU weights return.
No sealed HOLDOUT reads, no RSI status/score/cache writes, no submission.
This tests a one-head renderer; Bill's separate unchanged DE head is NOT reproduced.
"""
import common as job
import torch
torch.set_num_threads(2)
torch.set_num_interop_threads(1)
sys_guard_path = str(job.REPO/"model_billmB")
import sys
sys.path.insert(0,sys_guard_path)
import io_guard
io_guard.install()
import sys,os,types,importlib.util,gc
from pathlib import Path
import numpy as np
from infer import Ensemble
from network import blend

def main():
    assert job.json.loads((job.ROOT/'artifacts/weights_verified.json').read_text())['status']=='PASS'
    # Preserve the original RSI import context, required by frozen relative-path contracts.
    rsi=job.REPO/'rsi';os.chdir(rsi);sys.path.insert(0,str(rsi))
    # Read original frozen contracts, but do not execute RSI's mutating runtime.
    runtime=types.ModuleType('runtime');runtime.ROOT=rsi;runtime.json=job.json
    runtime.sha=job.sha;runtime.save=job.save
    def guard(extra=0):
        import shutil
        if shutil.disk_usage(job.ROOT).free-extra<25*1024**3:raise RuntimeError('25 GiB disk reserve')
    runtime.guard=guard
    def forbidden_status(*a,**k):raise RuntimeError('RSI status write forbidden')
    runtime.status=forbidden_status;sys.modules['runtime']=runtime
    sys.path.insert(0,str(rsi/'vendor'));sys.path.insert(0,str(rsi/'vendor/deps'))
    os.environ['NUMBA_CACHE_DIR']=str(job.ROOT/'artifacts/numba_cache');os.environ['TMPDIR']=str(job.ROOT/'artifacts/tmp')
    import benchmark,official_scoring as scoring,bounded_score,bounded_de,v5b_loco_baseline as baseline
    # Validate adapter against its original fixture, then route only its cache to this job.
    bounded_de.install();bounded_de.rt=types.SimpleNamespace(**{**vars(runtime),'ROOT':job.ROOT})
    contract=baseline.read_spec();results={}
    dependencies=[job.ROOT/n for n in ['evaluate.py','common.py','dataset.py','infer.py','network.py','data/manifest.json']]
    dependencies += [rsi/n for n in ['benchmark.py','official_scoring.py','bounded_score.py','bounded_de.py','v5b_loco_baseline.py','artifacts/v5b_loco_spec.json']]
    dependencies += sorted((job.ROOT/'data').glob('*.npz'))+sorted((job.ROOT/'checkpoints').glob('*/*_weights.pt'))
    for folder in [job.REPO/'public_data/model_v4_signature_banks',rsi/'artifacts/cd4_banks']:dependencies+=sorted(folder.glob('*.npz'))
    shared_inputs={str(p.relative_to(job.REPO)):job.sha(p) for p in dependencies}
    for task in ['H1_rich_search','CD4_depth192_search']:
        model,ctrl,real=benchmark.setup(task);assert model['split']=='search'
        output=job.ROOT/f'artifacts/eval/{task}';output.mkdir(parents=True,exist_ok=True)
        anchors_path=rsi/f'artifacts/scores/{task}/anchors.json';anchors=job.json.loads(anchors_path.read_text())
        assert anchors['real_sha256']==scoring.matrix_id(real)
        inputs=dict(shared=shared_inputs,task=task,targets=list(map(str,model['targets'])),genes=ctrl.var_names.astype(str).tolist(),
            control_sha256=scoring.matrix_id(ctrl),real_sha256=anchors['real_sha256'],anchors_sha256=job.sha(anchors_path))
        complete=output/'completed.json'
        if complete.exists():
            cached=job.json.loads(complete.read_text());assert cached['inputs']==inputs,'completed SEARCH inputs changed; retain old receipt and review before recomputing'
            for name,h in cached['files'].items():assert job.sha(output/name)==h,'completed SEARCH output changed'
            results[task]=job.json.loads((output/'score.json').read_text())
            print('REUSE COMPLETED SEARCH',task,flush=True);del ctrl,real;gc.collect();continue
        fold=model['context'];print('LOAD ENSEMBLE',task,flush=True);ensemble=Ensemble(fold)
        assert fold not in {b['family'] for b in ensemble.src.banks}
        print('ENSEMBLE LOADED',task,flush=True)
        # Frozen SEARCH panels have measured subsets (18077 / 17772 genes),
        # whereas training/inference uses the canonical 18533-gene axis.
        # Keep tasks/anchors/emitter native; zero-fill input controls for the
        # canonical network, then gather its response by gene name for scoring.
        genes=ctrl.var_names.to_numpy(str)
        assert len(set(genes))==len(genes)
        canonical=ensemble.src.genes
        lookup={g:i for i,g in enumerate(canonical)}
        assert len(lookup)==len(canonical) and all(g in lookup for g in genes)
        indices=np.array([lookup[g] for g in genes],dtype=np.int64)
        assert np.array_equal(canonical[indices],genes)
        b=np.asarray(ctrl.X.sum(0,dtype=np.float64)).ravel();b/=b.sum()
        canonical_control=np.zeros(len(canonical),dtype=np.float64)
        canonical_control[indices]=b
        assert np.array_equal(canonical_control[indices],b) and np.isclose(canonical_control.sum(),1)
        target=np.array(model['targets']);pred=ensemble.predict(target,canonical_control)
        learned=pred['learned'][:,indices]
        assert learned.shape==(len(target),len(genes)) and np.isfinite(learned).all()
        job.save(job.ROOT/f'artifacts/eval/{task}/gene_axis.json',dict(
            task_genes=len(genes),canonical_genes=len(canonical),
            mapping='zero-fill measured input controls to canonical axis; gather learned responses by name onto unchanged task axis',
            task_genes_sha256=__import__('hashlib').sha256('\n'.join(genes).encode()).hexdigest(),
            canonical_genes_sha256=__import__('hashlib').sha256('\n'.join(canonical).encode()).hexdigest(),
            canonical_indices=indices.tolist(),native_control_roundtrip_exact=True))
        # Frozen v5b counts supply the amplitude and emitter for this controlled sanity check.
        q,depth,bf,info=baseline.calibrated_profile(ctrl,task,contract)
        job.save(job.ROOT/f'artifacts/eval/{task}/amplitude.json',info)
        L=baseline.helpers.logbulk(bf);Lc=baseline.helpers.logbulk(b)[None,:]
        mixed=blend(L-Lc,learned,.3);new=baseline.helpers.from_logbulk(Lc+mixed)
        cfg=scoring.config(cache_real=str(job.ROOT/f'artifacts/eval/{task}/real_cache'))
        print('SCORING',task,flush=True)
        def iterator():
            for i,t in enumerate(target):
                if i%10==0:print('EMIT',task,i,'/',len(target),flush=True)
                yield t,baseline.emit(q,depth,new[i],model['context'],t)
        tidy,digest=bounded_score.score_iter(real,iterator(),cfg)
        raw=benchmark.aggregate(tidy);scaled,overall=scoring.scale(raw,anchors['baseline'],anchors['replicate'])
        tidy.write_csv(output/'tidy.csv')
        print('SCORED',task,overall,scaled,flush=True)
        rec=dict(task=task,raw=raw,scaled=scaled,overall=overall,prediction_sha256=digest,reference_sha256=anchors['real_sha256'],
            baseline=job.json.loads((rsi/f'artifacts/baseline_scores/{baseline.BASELINE_ID}/{task}/s20260910.json').read_text())['overall'],
            status='Official scorer, public LOCO sanity check ONLY; not an official challenge score; DE head coupling unresolved')
        job.save(output/'score.json',rec);results[task]=rec
        job.save(complete,dict(status='PASS completed SEARCH task checkpoint',inputs=inputs,
            files={name:job.sha(output/name) for name in ['score.json','tidy.csv','gene_axis.json','amplitude.json']}))
        del ensemble,ctrl,real;gc.collect()
    job.save(job.ROOT/'artifacts/local_scores.json',dict(tasks=results,equal_context_macro=float(np.mean([r['overall'] for r in results.values()])),official_score=None))
if __name__=='__main__':main()
