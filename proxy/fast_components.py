"""Exact official PDS/MSE from streamed cells and cached reference sufficient statistics.

No DE approximation: DE components are absent. Numerical acceptance comes from the
seven-arm full scorer, not from this script. --prepare is one-off reference setup.
"""
import os
os.environ.setdefault('PROXY_DE_THREADS','1');os.environ.setdefault('NUMBA_NUM_THREADS','1')
os.environ.setdefault('POLARS_MAX_THREADS','1');os.environ.setdefault('OMP_NUM_THREADS','1')
import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import time
import numpy as np
import pandas as pd
import polars as pl
import anndata as ad
from scipy import sparse
import run_v2 as runner
old=runner.old
from cell_eval2.moments import GroupMoments
from cell_eval2.prep import bulk_lognorm_means
from cell_eval2.catalog import CATALOG
from cell_eval2.scoring import score_one

P=Path(__file__).resolve().parent
KEYS=('pds_cosine','expr_mse_unbiased_capped_norm')
SHORT={'pds_cosine':'pds','expr_mse_unbiased_capped_norm':'mse'}


def pack(stats):
    order=np.array(sorted(stats));records=[stats[t] for t in order]
    return dict(perts=order,means=np.stack([r[0] for r in records]),counts=np.array([r[1] for r in records]),
                sumsq=np.array([r[2] for r in records]),jk=np.array([r[3] for r in records]))


def unpack(values):
    means=values['means'] if 'means' in values else bulk_lognorm_means(values['raw_sums'],50000.)
    return ({'bulk_lognorm':(values['perts'],means)},
            {'bulk_lognorm':GroupMoments(values['perts'],values['counts'],values['sumsq'],values['jk'])})


def prepare(dest):
    old.activate_scorer();runner.install_bounded_hash()
    real,ctrl,panel,note=runner.task(dest);cfg=old.scoring.config()
    anchor=json.loads((P/f'artifacts/scores_v2/{dest}/anchors.json').read_text())
    assert anchor['status']=='VALID' and anchor['real_sha256']==note['real_sha256'],'reference changed; rebuild calibration'
    folder=P/'artifacts/fast_reference';folder.mkdir(parents=True,exist_ok=True)
    stats={};sums={};labels=real.obs.target_gene.astype(str).to_numpy()
    for t in sorted(set(labels)):
        old.runtime.guard()
        part=real[labels==t].copy();stats[t]=old.bounded_score.stat(part,cfg)
        sums[t]=np.asarray(part.X.sum(0,dtype=np.int64)).ravel()
    arrays=pack(stats)
    arrays['raw_sums']=np.stack([sums[t] for t in arrays['perts']])
    assert np.array_equal(arrays['means'],bulk_lognorm_means(arrays['raw_sums'],cfg.bulk_target_sum))
    del arrays['means']  # Integer count sums compress much better than transformed fp64 means.
    arrays['genes']=real.var_names.to_numpy(str)
    np.savez_compressed(folder/(dest+'.npz'),**arrays)
    note.update(reference_npz_sha256=old.runtime.sha(folder/(dest+'.npz')),
                scorer_version=old.scoring.provenance()['version'])
    old.runtime.save(folder/(dest+'.json'),note)
    print('FAST_REFERENCE_READY',dest,flush=True)


def controls(dest,genes):
    path=old.RSI_H1 if dest=='H1' else P/f'tasks/{dest}.h5ad'
    note=json.loads((P/f'artifacts/fast_reference/{dest}.json').read_text())
    identity=note['source_identity'];current=path.stat()
    assert current.st_size==identity['size'] and current.st_mtime_ns==identity['mtime_ns'],'source task changed; regenerate reference and recheck calibration'
    a=ad.read_h5ad(path,backed='r');labels=a.obs.target_gene.astype(str).to_numpy()
    if dest=='H1':mask=a.obs.evaluation_role.astype(str).isin(['input_control','reference_control']).to_numpy()
    else:mask=labels==old.NT
    rows=np.flatnonzero(mask);parts=[]
    for lo in range(0,len(rows),400):parts.append(sparse.csr_matrix(a.X[rows[lo:lo+400]]))
    ctrl=ad.AnnData(sparse.vstack(parts,format='csr'),obs=a.obs.iloc[rows].copy(),var=a.var.copy());a.file.close()
    ctrl.obs['target_gene']=old.NT
    probability=note['thinning_probability']
    if probability<1:
        ctrl.X=sparse.csr_matrix(ctrl.X,dtype=np.int32)
        for block in note['control_rng_blocks']:
            rng=np.random.default_rng();rng.bit_generator.state=block['rng_state']
            lo=int(ctrl.X.indptr[block['control_start']]);hi=int(ctrl.X.indptr[block['control_stop']])
            for start in range(lo,hi,1_000_000):
                end=min(start+1_000_000,hi)
                ctrl.X.data[start:end]=rng.binomial(ctrl.X.data[start:end],probability).astype(np.int32)
        ctrl.X.eliminate_zeros()
    assert np.array_equal(ctrl.var_names.to_numpy(str),genes)
    return ctrl


def evaluate(dest,recipe,diagnostic,prediction_path=None):
    start=time.monotonic();old.activate_scorer();runner.install_bounded_hash()
    folder=P/'artifacts/fast_reference';values=dict(np.load(folder/(dest+'.npz')))
    note=json.loads((folder/(dest+'.json')).read_text());genes=values.pop('genes')
    assert old.runtime.sha(folder/(dest+'.npz'))==note['reference_npz_sha256']
    reference,rm=unpack(values);panel=np.array([t for t in values['perts'] if t!=old.NT]);ctrl=controls(dest,genes)
    old.model.configure(old.model.SHARED_BANKS,[old.FAMILY[dest]],old.runtime.VCC/'model_v6/artifacts/promoter_pairs_v4.csv')
    backed_prediction=None
    if prediction_path:
        backed_prediction=ad.read_h5ad(prediction_path,backed='r')
        assert np.array_equal(backed_prediction.var_names.to_numpy(str),genes),'prediction gene axis mismatch'
        labels=backed_prediction.obs.target_gene.astype(str).to_numpy()
        assert set(labels)-{old.NT}==set(panel),'prediction panel mismatch'
        for t in panel:assert int((labels==t).sum())==400,'400 predicted cells per target required'
        stream=((t,sparse.csr_matrix(backed_prediction.X[np.flatnonzero(labels==t)])) for t in panel)
        recipe=prediction_path.stem
    elif recipe in ('v1','v2b','v3b'):
        import sys
        saved=sys.modules.get('pipeline');sys.path.insert(0,str(P/'legacy'));import legacy_recipes
        if saved:sys.modules['pipeline']=saved
        # Preserve the frozen task's prediction order, including hash-ordered CD4 panel.
        if dest!='H1':panel=pd.read_csv(P/f'tasks/{dest}_panel.csv').target_gene.astype(str).to_numpy()
        stream=legacy_recipes.predict(ctrl,list(panel),dest,recipe)
    else:
        if dest!='H1':panel=pd.read_csv(P/f'tasks/{dest}_panel.csv').target_gene.astype(str).to_numpy()
        stream=old.pipeline.Context(recipe,panel,genes,ctrl,dest,0,workers=1).stream()
    cfg=old.scoring.config();stats={old.NT:old.bounded_score.stat(ctrl,cfg)}
    prediction_digest=hashlib.sha256()
    for t,x in stream:
        old.runtime.guard();old.bounded_score.validate(x)
        for value in (x.data,x.indices,x.indptr):prediction_digest.update(memoryview(np.ascontiguousarray(value)).cast('B'))
        prediction_digest.update(t.encode())
        stats[t]=old.bounded_score.stat(old.bounded_score.make(x,t,genes),cfg)
    if backed_prediction is not None:backed_prediction.file.close()
    predicted,pm=unpack(pack(stats));names=list(old.bounded_score.resolve_metrics(cfg.metrics,version=cfg.version)[0])
    tidy=pl.DataFrame(old.bounded_score.run.dispatch_anndata_metrics(names,predicted,reference,genes,cfg,
        comparator='bulk_lognorm',pred_moments=pm,real_moments=rm,driver='Proxy exact reference sufficient statistics'),schema=old.bounded_score.SCHEMA)
    raw=old.aggregate(tidy);anchor=json.loads((P/f'artifacts/scores_v2/{dest}/anchors.json').read_text())
    assert anchor['real_sha256']==note['real_sha256'] and anchor['status']=='VALID'
    scaled={m:float(score_one(raw[m],anchor['baseline'][m],replace(CATALOG[m].scoring,anchor=anchor['replicate'][m]))) for m in KEYS}
    calibration_path=P/'artifacts/calibration_v2.json'
    calibration=json.loads(calibration_path.read_text()) if calibration_path.exists() else {}
    verification_path=P/'artifacts/calibration_v2_verification.json'
    verified=json.loads(verification_path.read_text()) if verification_path.exists() else {}
    checked=verified.get('calibration_sha256')==old.runtime.sha(calibration_path) if calibration_path.exists() else False
    trusted=[m for m in KEYS if checked and SHORT[m] in calibration.get('trusted_components',[])]
    if not trusted and not diagnostic:raise SystemExit('No fast component has passed calibration; use --diagnostic for untrusted numerical checks')
    result=dict(destination=dest,recipe=recipe,scaled=scaled,raw=raw,macro_eligible_components=[SHORT[m] for m in trusted],
                seconds=time.monotonic()-start,task_sha256=note['real_sha256'],
                anchor_sha256=old.runtime.sha(P/f'artifacts/scores_v2/{dest}/anchors.json'),
                prediction_sha256=prediction_digest.hexdigest(),
                note='PDS/MSE only; destination values alone are diagnostic. Only the fixed two-destination macro has been calibrated; no official-total prediction.')
    out=P/f'artifacts/fast_scores/{dest}';out.mkdir(parents=True,exist_ok=True)
    old.runtime.save(out/(recipe+'.json'),result);print(json.dumps(result,indent=2))


def main():
    p=argparse.ArgumentParser();p.add_argument('--dest',required=True,choices=('H1','CD4'))
    p.add_argument('--recipe',choices=('v1','v2b','v3b','v4a','v4b','v5a','v5b'));p.add_argument('--prepare',action='store_true');p.add_argument('--diagnostic',action='store_true')
    p.add_argument('--prediction',type=Path,help='arbitrary candidate counts on this public destination/panel, 400 cells per target')
    a=p.parse_args()
    from execution_v2 import install_guard
    install_guard(old.runtime)
    if a.prepare:prepare(a.dest)
    else:
        assert bool(a.recipe)!=bool(a.prediction),'provide exactly one of --recipe/--prediction'
        evaluate(a.dest,a.recipe,a.diagnostic,a.prediction)


if __name__=='__main__':main()
