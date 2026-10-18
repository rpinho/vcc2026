"""Official five-split/depth-matched proxy. Stream cells, retain tiny reductions only.

One invocation per arm keeps frozen legacy import aliases isolated.
"""
import os
os.environ.setdefault('PROXY_DE_THREADS','1')
os.environ.setdefault('NUMBA_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
os.environ.setdefault('POLARS_MAX_THREADS','1')
import argparse
from dataclasses import replace
import gc
import json
import logging
import hashlib
import fcntl
from pathlib import Path
import time
import run as old
import numpy as np
import pandas as pd
import polars as pl
from scipy import sparse
from cell_eval2.anchor import compute_replicate_anchor
from cell_eval2.competition import competition_digest, derived_seeds

P=Path(__file__).resolve().parent


def thin_counts_with_control_states(x,probability,control_mask):
    """Same full-matrix binomial stream; checkpoint only control block RNG states."""
    changes=np.diff(np.r_[False,control_mask,False].astype(np.int8))
    starts=np.flatnonzero(changes==1);stops=np.flatnonzero(changes==-1)
    blocks=[];offset=0
    for start,stop in zip(starts,stops):
        blocks.append(dict(data_start=int(x.indptr[start]),data_stop=int(x.indptr[stop]),
                           control_start=offset,control_stop=offset+int(stop-start)))
        offset+=int(stop-start)
    by_start={b['data_start']:b for b in blocks}
    boundaries=sorted(set([0,len(x.data),*range(0,len(x.data),1_000_000),
                           *(b['data_start'] for b in blocks),*(b['data_stop'] for b in blocks)]))
    rng=np.random.default_rng(20261018)
    for lo,hi in zip(boundaries[:-1],boundaries[1:]):
        if lo in by_start:by_start[lo]['rng_state']=rng.bit_generator.state
        x.data[lo:hi]=rng.binomial(x.data[lo:hi],probability).astype(np.int32)
    return blocks


def bounded_matrix_id(a):
    """Same byte digest as the frozen adapter, without multi-GB tobytes copies."""
    h=hashlib.sha256();x=a.X.tocsr()
    for value in (x.data,x.indices,x.indptr):
        h.update(memoryview(np.ascontiguousarray(value)).cast('B'))
    h.update('\0'.join(a.obs.target_gene.astype(str)).encode())
    h.update('\0'.join(a.var_names).encode())
    return h.hexdigest()


def install_bounded_hash():
    from anndata import AnnData
    probe=AnnData(sparse.csr_matrix([[0,1,2],[3,0,1]],dtype=np.int32),
                  obs=pd.DataFrame({'target_gene':['AA','non-targeting']},index=['a','b']),
                  var=pd.DataFrame(index=['AA','BB','CC']))
    assert bounded_matrix_id(probe)==old.scoring.matrix_id(probe)
    old.scoring.matrix_id=bounded_matrix_id


def task(dest):
    print('LOAD_TASK',dest,flush=True)
    from mapped_task import load
    real,panel,mapping=load(dest)
    ctrl=real[real.obs.target_gene.to_numpy()==old.NT].copy()
    labels=real.obs.target_gene.astype(str).to_numpy()
    note=dict(destination=dest,base_matrix_sha256=mapping['base_matrix_sha256'],seed=20261018,
              source_identity=mapping['source_identity'],
              controls=ctrl.n_obs,authentic_ntc_ids=int(ctrl.obs.ntc_id.nunique()) if 'ntc_id' in ctrl.obs else None,
              residual_mismatch='Fewer than 18400 independent controls; authentic 46-NTC identity structure unavailable. No fabricated cells or guide IDs.')
    if dest=='HCT116':
        count=pd.Series(labels).value_counts()
        panel=np.array(sorted(t for t in panel if count[t]==400))
        real=real[np.isin(labels,[old.NT,*panel])].copy()
        note['subset_rule']='All original panel targets with 400 independent real cells; fixed before scores'
    probability=mapping['thinning_probability']
    note.update(native_control_median_umi=mapping['native_control_median_umi'],thinning_probability=probability,
                control_rng_blocks=mapping['control_rng_blocks'])
    lab=real.obs.target_gene.astype(str).to_numpy();ctrl=real[lab==old.NT].copy()
    counts=pd.Series(lab[lab!=old.NT]).value_counts()
    if 'source_row' in real.obs:
        assert real.obs.source_row.nunique()==real.n_obs,'duplicated reference source row'
    note.update(n_targets=len(panel),genes=real.n_vars,cells_per_target_min=int(counts.min()),
                cells_per_target_median=float(counts.median()),targets_ge400=int((counts>=400).sum()),
                control_median_umi=float(np.median(np.asarray(ctrl.X.sum(1)).ravel())),real_sha256=old.scoring.matrix_id(real))
    print('TASK_READY',dest,note,flush=True)
    return real,ctrl,panel,note


def anchors(dest,real,folder,note):
    path=folder/'anchors.json'
    if path.exists():
        saved=json.loads(path.read_text());assert saved['real_sha256']==note['real_sha256']
        if saved['status']!='VALID':raise ValueError(saved['error'])
        return saved
    from baseline_protocol import mean_response
    t0=time.monotonic();cfg=replace(old.scoring.config(cache_real=None),num_threads=int(os.environ.get('PROXY_DE_THREADS','1')))
    previous_path=P/'artifacts/scores'/dest/'anchors.json'
    previous=json.loads(previous_path.read_text()) if previous_path.exists() else {}
    same_real=previous.get('real_sha256')==note['real_sha256']
    if same_real:
        base=previous['baseline'];desc=previous['baseline_construction']
        print('REUSE_BASELINE_EXACT_MATRIX',previous_path,flush=True)
    else:
        mean,desc=mean_response(real)
        tidy,_=old.bounded_score.score_iter(real,old.bounded_score.baseline_iter(real,mean,cells=400),replace(cfg,allow_fractional_counts=True))
        base=old.aggregate(tidy)
        tidy.write_parquet(folder/'baseline_tidy.parquet')
    base_seconds=time.monotonic()-t0
    old.runtime.save(folder/'baseline_raw.json',dict(raw=base,real_sha256=note['real_sha256'],seconds=base_seconds))
    print('BASELINE',dest,base,base_seconds,flush=True)
    try:
        if same_real and previous.get('n_splits')==3:
            splits,frame=incremental_anchor(real,cfg,previous,previous_path,folder)
        else:
            splits,frame=compute_replicate_anchor(real,config=cfg,n_splits=5,base_seed=0)
    except Exception as error:
        old.runtime.save(path,dict(status='INVALID',error=str(error),baseline=base,real_sha256=note['real_sha256'],
                                  n_splits=5,task=note,baseline_seconds=base_seconds,
                                  anchor_seconds=time.monotonic()-t0-base_seconds))
        raise
    splits.write_parquet(folder/'anchor_splits.parquet');frame.write_parquet(folder/'anchor_agg.parquet')
    rep=dict(zip(frame['metric'].to_list(),frame['replicate'].to_list()))
    from cell_eval2.catalog import CATALOG
    headroom={m:float(base[m]-rep[m] if CATALOG[m].scoring.direction=='lower' else rep[m]-base[m]) for m in old.scoring.SIX}
    saved=dict(baseline=base,replicate=rep,real_sha256=note['real_sha256'],n_splits=5,
               derived_seeds=list(derived_seeds()),baseline_construction=desc,baseline_seconds=base_seconds,
               anchor_seconds=time.monotonic()-t0-base_seconds,task=note,competition_rule_digest=competition_digest(),headroom=headroom)
    if (folder/'anchor_cached_prefix.json').exists():
        saved['cached_prefix']=json.loads((folder/'anchor_cached_prefix.json').read_text())
    try:
        scaled,overall=old.scoring.scale(base,base,rep)
        assert all(abs(v)<1e-12 for v in scaled.values()) and abs(overall)<1e-12
        scaled,overall=old.scoring.scale(rep,base,rep)
        assert all(abs(v-1)<1e-12 for v in scaled.values()) and abs(overall-1)<1e-12
        saved['status']='VALID'
    except Exception as error:
        saved.update(status='INVALID',error=str(error))
    old.runtime.save(path,saved)
    print('ANCHORS',dest,saved['status'],rep,flush=True)
    if saved['status']!='VALID':raise ValueError(saved['error'])
    return saved


def incremental_anchor(real,cfg,previous,previous_path,folder):
    """Extend the identical old 3-split mean with the missing official seeds 3/4.

    No fictitious per-split rows: only new splits are saved. The old prefix's
    individual spread is unavailable, so merged SD/min/max remain null.
    """
    from fractions import Fraction
    import cell_eval2.anchor as official_anchor
    seeds=list(derived_seeds());original=official_anchor._derive_seeds
    print('REUSE_ANCHOR_PREFIX',previous_path,'new seeds',seeds[3:],flush=True)
    official_anchor._derive_seeds=lambda base,n: seeds[3:]
    try:
        splits,frame=compute_replicate_anchor(real,config=cfg,n_splits=2,base_seed=0)
    finally:
        official_anchor._derive_seeds=original
    splits=splits.with_columns((pl.col('split_index')+3).alias('split_index'))
    means=[]
    for metric in frame['metric'].to_list():
        values=splits.filter(pl.col('metric')==metric)['value'].to_list()
        assert len(values)==2
        exact=(3*Fraction(previous['replicate'][metric])+sum(Fraction(v) for v in values))/5
        means.append(float(exact))
    frame=frame.with_columns(pl.Series('replicate',means))
    for name in ('replicate_sd','replicate_min','replicate_max','n_perturbations_min','n_perturbations_max'):
        frame=frame.with_columns(pl.lit(None,dtype=frame.schema[name]).alias(name))
    record=dict(n_cached=3,cached_seed_prefix=seeds[:3],new_seeds=seeds[3:],
                source=str(previous_path.relative_to(P)),source_sha256=old.runtime.sha(previous_path),
                cached_mean=previous['replicate'],
                qualification='Five-split mean extended from the identical archived three-split mean. Individual prefix split values/SD are unavailable, so archival rounding cannot be independently reconstructed. Combination is exact on the stored inputs and rounded once; no statistical re-fit or synthetic split rows.')
    old.runtime.save(folder/'anchor_cached_prefix.json',record)
    return splits,frame


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--dest',required=True,choices=old.FAMILY)
    parser.add_argument('--recipe',choices=('v1','v2b','v3b','v4a','v4b','v5a','v5b'))
    parser.add_argument('--recipes',help='comma-separated arms, one task load, sequential predictions')
    parser.add_argument('--prediction',type=Path,help='arbitrary candidate counts on this public proxy destination')
    parser.add_argument('--anchors-only',action='store_true')
    args=parser.parse_args();provenance=old.activate_scorer();install_bounded_hash()
    from execution_v2 import install_guard,install_mwu_cache,install_sorted_mwu,install_stat_cache
    install_guard(old.runtime)
    sorted_fixture=install_sorted_mwu();mwu_counts=install_mwu_cache(old.runtime.guard)
    stat_counts=install_stat_cache(old.bounded_score)
    provenance['sorted_mwu_backend']=sorted_fixture['status']
    provenance['execution_wrapper_sha256']=old.runtime.sha(P/'execution_v2.py')
    provenance['mapped_task_sha256']=old.runtime.sha(P/'mapped_task.py')
    folder=P/'artifacts/scores_v2'/args.dest;folder.mkdir(parents=True,exist_ok=True)
    # Same adapter, but a separate regenerable cache avoids touching original numerical evidence.
    old.runtime.ROOT=P/'v2_cache';(old.runtime.ROOT/'artifacts').mkdir(parents=True,exist_ok=True)
    real,ctrl,panel,note=task(args.dest)
    old.runtime.guard()
    with (folder/'anchors.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        old.runtime.save(folder/'task.json',note)
        anc=anchors(args.dest,real,folder,note)
    if args.anchors_only:
        print('MWU_CACHE',mwu_counts,flush=True);return
    if args.prediction:
        assert not args.recipe and not args.recipes
        score_prediction(args.dest,args.prediction,real,panel,anc,note,provenance)
        return
    recipes=args.recipes.split(',') if args.recipes else [args.recipe]
    assert all(r in ('v1','v2b','v3b','v4a','v4b','v5a','v5b') for r in recipes)
    for recipe in recipes:
        score_recipe(args.dest,recipe,real,ctrl,panel,anc,folder,note,provenance)
        gc.collect()
    print('MWU_CACHE',mwu_counts,'STAT_CACHE',stat_counts,flush=True)


def score_recipe(dest,recipe,real,ctrl,panel,anc,folder,note,provenance):
    previous_path=P/'artifacts/scores'/dest/(recipe+'.json')
    previous_anchor=P/'artifacts/scores'/dest/'anchors.json'
    if (previous_path.exists() and previous_anchor.exists() and
        json.loads(previous_anchor.read_text()).get('real_sha256')==note['real_sha256']):
        record=json.loads(previous_path.read_text())
        scaled,overall=old.scoring.scale(record['raw'],anc['baseline'],anc['replicate'])
        record.update(scaled=scaled,overall=overall,cached_raw_source=str(previous_path.relative_to(P)),
                      cached_raw_sha256=old.runtime.sha(previous_path),runtime_is_prior_run=True,
                      reuse='Exact reference matrix identity; unchanged prediction recipe; only replicate scale rebuilt from 3 to 5 splits')
        print('REUSED_RAW',dest,recipe,overall,flush=True)
    else:
        record=old.score(dest,recipe,real,ctrl,panel,anc,folder,workers=1)
    record.update(task_sha256=note['real_sha256'],scorer=provenance,anchor_splits=5,
                  anchor_sha256=old.runtime.sha(folder/'anchors.json'))
    old.runtime.save(folder/(recipe+'.json'),record)


def score_prediction(dest,path,real,panel,anc,note,provenance):
    start=time.monotonic();prediction=old.ad.read_h5ad(path,backed='r')
    assert np.array_equal(prediction.var_names.to_numpy(str),real.var_names.to_numpy(str))
    labels=prediction.obs.target_gene.astype(str).to_numpy()
    assert set(labels)-{old.NT}==set(panel),'candidate must reproduce the frozen public proxy panel'
    for target in panel:assert int((labels==target).sum())==400,'exactly 400 predicted cells required'
    gen=((t,sparse.csr_matrix(prediction.X[np.flatnonzero(labels==t)])) for t in panel)
    cfg=replace(old.scoring.config(),num_threads=int(os.environ.get('PROXY_DE_THREADS','1')))
    tidy,digest=old.bounded_score.score_iter(real,gen,cfg)
    prediction.file.close();raw=old.aggregate(tidy);scaled,overall=old.scoring.scale(raw,anc['baseline'],anc['replicate'])
    out=P/'artifacts/candidate_scores'/path.stem/dest;out.mkdir(parents=True,exist_ok=True)
    record=dict(destination=dest,candidate=path.stem,raw=raw,scaled=scaled,overall=overall,seconds=time.monotonic()-start,
                task_sha256=note['real_sha256'],anchor_sha256=old.runtime.sha(P/f'artifacts/scores_v2/{dest}/anchors.json'),
                scorer=provenance,prediction_sha256=digest,note='Public proxy destination; not hidden contexts A/B/C or an official score')
    tidy.write_parquet(out/'tidy.parquet');old.runtime.save(out/'score.json',record);print(json.dumps(record,indent=2))


if __name__=='__main__':
    logging.basicConfig(level=logging.INFO,format='%(message)s')
    main()
