"""Independent projection placement, PCA, unchanged mean and emitter audit.

Reads SEARCH input controls and source summaries only; no response statistics.
"""
from pathlib import Path
from datetime import datetime, timezone
from itertools import islice
import importlib.util
import gc
import hashlib
import numpy as np
from scipy import sparse
import runtime as rt
import atlas_loco_baseline as atlas
import deflate50_loco_baseline as baseline


def independent_project(effect, axes):
    removed = np.zeros(effect.shape, dtype=float)
    for column in axes.T:
        coefficients = np.einsum('ij,j->i', effect, column)
        removed += np.multiply.outer(coefficients, column)
    return effect-removed


def fixture():
    rng = np.random.default_rng(45)
    counts = rng.poisson(np.geomspace(1,35,12),size=(90,12))
    axes, _ = baseline.control_axes(counts,seed=1,k=4)
    y = np.log1p(50000*counts/counts.sum(1,keepdims=True)).astype(np.float32)
    y -= y.mean(0,keepdims=True)
    _,_,vt = np.linalg.svd(y,full_matrices=False)
    np.testing.assert_allclose(axes@axes.T,vt[:4].T@vt[:4],rtol=0,atol=3e-6)
    delta = rng.normal(size=(5,12))
    got = baseline.project(delta,axes)
    np.testing.assert_allclose(got,independent_project(delta,axes),rtol=1e-12,atol=1e-12)
    np.testing.assert_allclose(got@axes,0,rtol=0,atol=3e-6)
    np.testing.assert_allclose(baseline.project(got,axes),got,rtol=0,atol=3e-6)
    # Demonstrate why pre-amplitude/clipping placement is part of the contract.
    raw = np.array([[8.,-2.,1.]])
    v = np.array([[1.],[1.],[0.]])/np.sqrt(2)
    assert not np.allclose(np.clip(baseline.project(raw,v)*.3,-1,1),
                           baseline.project(np.clip(raw*.3,-1,1),v))
    return {'exact_SVD_small_projector':'PASS','scalar_projection':'PASS',
        'orthogonality_idempotence':'PASS','placement_noncommutativity':'PASS'}


def main():
    spec = baseline.read_spec()
    checks = fixture()
    frozen = rt.json.loads(Path('artifacts/frozen_metadata_index.json').read_text())
    for path,item in frozen['entries'].items(): assert rt.sha(path)==item['sha256']
    assert rt.sha(frozen['archive']) == frozen['archive_sha256']
    module_spec=importlib.util.spec_from_file_location('_audit_archived_axes','artifacts/reproduction_inputs/v5b_round33/v5.py')
    archived=importlib.util.module_from_spec(module_spec); module_spec.loader.exec_module(archived)
    records={}
    for task in baseline.TASKS:
        print('CONTROL ONLY full-panel audit',task,flush=True)
        ctrl=baseline.input_controls(task)
        axes,singular=baseline.control_axes(ctrl.X,spec['axis_seed'])
        old_axes,old_s=archived.control_axes(ctrl.X,k=50,seed=spec['axis_seed'])
        assert np.array_equal(axes,old_axes) and np.array_equal(singular,old_s)
        np.testing.assert_allclose(axes.T@axes,np.eye(50),rtol=0,atol=4e-6)
        model,_=atlas.release_modules()
        q0,d0,unmodified,cov0=atlas.heads(ctrl,task,spec)
        actual=baseline.heads(ctrl,task,spec)
        assert np.array_equal(actual[0],q0) and np.array_equal(actual[1],d0)
        assert np.array_equal(actual[2][0],unmodified[0])
        original=model.desired_mean
        stage={}
        def independent(control,effect,*,space,amplitude,clip):
            if space=='bulk_delta':
                after=independent_project(effect,axes)
                np.testing.assert_allclose(after,baseline.project(effect,axes),rtol=2e-9,atol=2e-12)
                residual=float(np.max(np.abs(after@axes)))
                assert residual < max(1e-8,float(np.max(np.abs(effect@axes)))*1e-5)
                assert np.linalg.norm(after) <= np.linalg.norm(effect)*(1+1e-10)
                stage.update({'removed_norm':float(np.linalg.norm(effect-after)),
                    'max_projection_residual':residual,'amplitude':amplitude,'clip':clip})
                effect=after
            return original(control,effect,space=space,amplitude=amplitude,clip=clip)
        model.desired_mean=independent
        try: expected=atlas.heads(ctrl,task,spec)
        finally: model.desired_mean=original
        assert stage['amplitude']==.3 and stage['clip']==3
        np.testing.assert_allclose(actual[2][1],expected[2][1],rtol=2e-11,atol=2e-15)
        assert np.array_equal(actual[2][0],expected[2][0])
        assert all(np.array_equal(a,b) for a,b in zip(actual[3],cov0))
        for p in actual[2]:
            assert np.isfinite(p).all() and np.all(p>=0)
            np.testing.assert_allclose(p.sum(1),1,rtol=1e-12,atol=1e-12)
        emissions=[]
        iterator=baseline.predictions(ctrl,task,spec)
        for i,(target,x) in enumerate(islice(iterator,2)):
            seed=model.seed_for(spec['tasks'][task]['destination']+':'+target,spec['seed'])
            replay=model.dual_moment_counts(q0,expected[2][0][i],actual[2][1][i],depths=d0,seed=seed)
            # Replay uses the production profile to isolate emitter identity;
            # separate full-profile arithmetic was independently tested above.
            assert np.array_equal(x.toarray(),replay)
            assert np.array_equal(replay.sum(1),d0) and replay.dtype.kind in 'iu' and np.all(replay>=0)
            emissions.append({'target':target,'sha256_dense':hashlib.sha256(replay.tobytes()).hexdigest(),'exact_replay':True})
        iterator.close()
        records[task]={'status':'PASS','input_controls':ctrl.n_obs,'axis_shape':list(axes.shape),
            'archived_axis_replay_bitwise_equal':True,'mean_channel_bitwise_unchanged':True,
            'full_bulk_head_max_abs_error':float(np.max(np.abs(actual[2][1]-expected[2][1]))),
            'coverage_bitwise_unchanged':True,'stage':stage,'emissions':emissions}
        del ctrl,axes,old_axes,actual,expected,unmodified,q0,replay
        gc.collect()
    baseline.read_spec()
    rt.save(baseline.AUDIT,{'utc':datetime.now(timezone.utc).isoformat(),'status':'PASS',
        'spec_sha256':rt.sha(baseline.SPEC),'verifier_sha256':rt.sha(__file__),
        'tasks':records,'fixture':checks,'frozen_metadata_members_unchanged':len(frozen['entries']),
        'external_official_label':None,'holdout_scored':False,
        'scope':'Description-based available-input projection, numerical fidelity and control-only isolation; no exact Becker release or biological utility claim.'})
    print('PASS deflate50 adapter',flush=True)


if __name__=='__main__': main()
