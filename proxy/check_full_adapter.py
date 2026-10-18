"""Unmodified official scorer/anchors versus bounded execution and caching, bit for bit."""
import os
for k in ('NUMBA_NUM_THREADS','OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','POLARS_MAX_THREADS'):
    os.environ.setdefault(k,'1')
import hashlib,json
import tempfile
from pathlib import Path
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
import polars as pl
from dataclasses import replace
from cell_eval2 import EvalConfig,compute_metrics,aggregate_metrics
from cell_eval2.anchor import compute_replicate_anchor
from execution_v2 import install_sorted_mwu,install_mwu_cache,install_stat_cache,numeric_adapter_sha256


def main():
    rng=np.random.default_rng(18);genes=[f'G{i}' for i in range(160)];labels=[];real=[];pred=[]
    for i,t in enumerate(['non-targeting','G0','G1','G2']):
        mu=np.full(160,8.);n=120 if i==0 else 60
        if i:mu[20*i:20*(i+1)]=24;mu[:10]=3
        real.append(rng.poisson(mu,(n,160)));pred.append(rng.poisson(mu*.95+.3,(n,160)));labels.extend([t]*n)
    def matrix(blocks):return ad.AnnData(sparse.csr_matrix(np.vstack(blocks)),obs=pd.DataFrame({'target_gene':labels},index=[str(i) for i in range(len(labels))]),var=pd.DataFrame(index=genes))
    r=matrix(real);p=matrix(pred);cfg=EvalConfig.from_preset('vcc2026')
    cfg=replace(cfg,pert_col='target_gene',device='cpu',num_threads=1,gather_threads=1,de=replace(cfg.de,backend='pdex'))
    original=compute_metrics(p,r,config=cfg).sort(['metric','perturbation']);osplits,oanchor=compute_replicate_anchor(r,config=cfg,n_splits=2)
    import run_v2 as runner
    old=runner.old
    from baseline_protocol import mean_response
    profile,_=mean_response(r)
    blocks=[];blabel=[]
    for target in ['G0','G1','G2']:
        blocks.append(sparse.csr_matrix(np.tile(profile.astype(np.float32),(400,1))));blabel.extend([target]*400)
    blocks.append(r.X[np.array(labels)=='non-targeting']);blabel.extend(['non-targeting']*120)
    baseline=ad.AnnData(sparse.vstack(blocks,format='csr').astype(np.float32),obs=pd.DataFrame({'target_gene':blabel},index=[str(i) for i in range(len(blabel))]),var=r.var.copy())
    bcfg=replace(cfg,allow_fractional_counts=True)
    obase=compute_metrics(baseline,r,config=bcfg).sort(['metric','perturbation'])
    old.activate_scorer();runner.install_bounded_hash()
    fixture=install_sorted_mwu();mwu=install_mwu_cache(lambda:None);stat=install_stat_cache(old.bounded_score)
    (Path(__file__).resolve().parent/'v2_cache').mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parent/'v2_cache',prefix='full_fixture_') as scratch:
        old.runtime.ROOT=Path(scratch);(Path(scratch)/'artifacts').mkdir()
        changed=old.bounded_score.resident(p,r,config=cfg).sort(['metric','perturbation'])
        splits,anchor=compute_replicate_anchor(r,config=cfg,n_splits=2)
        newbase=old.bounded_score.resident(baseline,r,config=bcfg).sort(['metric','perturbation'])
    assert original.equals(changed),'full tidy frame differs'
    assert osplits.equals(splits),'anchor splits differ'
    assert oanchor.equals(anchor),'anchor aggregate differs'
    assert obase.equals(newbase),'400-cell tiled baseline differs'
    root=Path(__file__).resolve().parent
    (root/'artifacts').mkdir(parents=True,exist_ok=True)
    result=dict(status='PASS bitwise official full scorer and two-split anchors',tidy_rows=original.height,
                split_rows=splits.height,anchor_rows=anchor.height,sorted_mwu_fixture=fixture,
                baseline_rows=newbase.height,mwu_cache=mwu,stat_cache=stat,
                numeric_adapter_sha256=numeric_adapter_sha256(),
                wrapper_sha256=hashlib.sha256((root/'execution_v2.py').read_bytes()).hexdigest())
    (root/'artifacts/full_adapter_fixture.json').write_text(json.dumps(result,indent=2));print(result['status'],flush=True)


if __name__=='__main__':main()
