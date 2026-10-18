"""Official cell-eval2 metrics and scaling; no metric formula is reimplemented.
The separately verified CPU DE memory adapter is activated by benchmark.py.
"""
import runtime as rt
import hashlib,json
from dataclasses import replace
from pathlib import Path
import numpy as np,polars as pl
import cell_eval2
from cell_eval2 import EvalConfig,compute_metrics,aggregate_metrics
import cell_eval2.run as run
import cell_eval2.anchor as anchor
import cell_eval2.lfc_nmae_ref as lfc_ref
from cell_eval2.score import _replicate_entries,_reference_column
from cell_eval2.scoring import score_one
SIX=['pds_cosine','expr_mse_unbiased_capped_norm','de_wilcoxon_lfc_nmae','de_wilcoxon_direction_fidelity_yield_raw','de_wilcoxon_direction_reach_raw','de_wilcoxon_sig_jaccard']
ORIGINAL_DE=run._compute_de_side

def config(**kwargs):
 c=EvalConfig.from_preset('vcc2026');return replace(c,pert_col='target_gene',device='cpu',num_threads=2,gather_threads=1,de=replace(c.de,backend='pdex'),**kwargs)

def matrix_id(a):
 h=hashlib.sha256();x=a.X.tocsr()
 for value in [x.data,x.indices,x.indptr]:h.update(value.tobytes())
 h.update('\0'.join(a.obs.target_gene.astype(str)).encode());h.update('\0'.join(a.var_names).encode());return h.hexdigest()

def aggregate(pred,real,cfg=None):
 result=compute_metrics(pred,real,config=cfg or config());agg=aggregate_metrics(result);return dict(zip(agg['metric'].to_list(),agg['mean'].to_list()))

def scale(raw,baseline,replicate):
 entries=_replicate_entries(baseline,pl.DataFrame({'metric':list(replicate),'replicate':list(replicate.values())}));column=_reference_column(SIX+['avg_score'],[raw[m] for m in SIX],SIX,entries,column='from_replicate',label='RSI measured local reference');values=dict(zip(SIX,column.to_list()[:-1]));assert all(np.isfinite(list(values.values())));return values,float(column[-1])

def provenance():
 package=Path(cell_eval2.__file__).parent
 from importlib.metadata import version,PackageNotFoundError
 versions={}
 for name in ['cell-eval2','pdex','numba-mwu','numpy','scipy','anndata','polars','pandas','h5py','numba']:
  try:versions[name]=version(name)
  except PackageNotFoundError:versions[name]=None
 return dict(version=cell_eval2.__version__,packages=versions,official_files={str(p.relative_to(package)):rt.sha(p) for p in package.rglob('*.py')},config=config().to_dict(),execution='Official compute_metrics/aggregate_metrics/replicate anchor/scaling; Original official CPU DE in target-complete pieces, preserving row order; passes bitwise fixture')
