"""Official CPU DE called on target-complete pieces; original row order preserved.
There are no metric or DE formulas here. Activation requires bitwise fixture equality.
The official full-panel metric dispatcher, jackknife cap, ranking, and aggregation stay whole.
"""
import runtime as rt
import hashlib,json
from pathlib import Path
import numpy as np,polars as pl
import cell_eval2.run as run,cell_eval2.anchor as anchor,cell_eval2.lfc_nmae_ref as ref
ORIGINAL=run._compute_de_side

def compute(adata,*,cfg,fp,store,side,reference_adata=None):
 assert cfg.de.backend=='pdex' and cfg.de.fdr_scope=='per_pert' and cfg.target_sum==1e6 and reference_adata is None
 from official_scoring import matrix_id
 content=matrix_id(adata);params=cfg.to_dict()
 # Matrix already contains its exact reference controls. These fields are operational here.
 normalized=dict(params)
 for k in ['cache_real','cache_pred','outdir']:normalized[k]=None
 normalized['num_threads']=1;normalized['gather_threads']=1;normalized['control_source']='pred'
 normalized['allow_fractional_counts']=bool(cfg.allow_fractional_counts and side=='pred')
 key=hashlib.sha256((content+json.dumps(normalized,sort_keys=True,default=str)+':official-de-order-preserving-v4').encode()).hexdigest();path=rt.ROOT/'artifacts/de_cache'/f'{key}.parquet';path.parent.mkdir(exist_ok=True)
 # Reuse exact old-config entries from the first attempt; no numerical recalculation.
 if not path.exists():
  for nt in [1,2]:
   old=dict(params);old['num_threads']=nt;oldkey=hashlib.sha256((content+json.dumps(old,sort_keys=True,default=str)+':official-de-order-preserving-v3').encode()).hexdigest();previous=path.with_name(oldkey+'.parquet')
   if previous.exists():
    frame=pl.read_parquet(previous);frame.write_parquet(path);return frame
 if path.exists():return pl.read_parquet(path)
 labels=adata.obs[cfg.pert_col].astype(str).to_numpy();targets=sorted(set(labels)-{cfg.control});frames=[]
 for lo in range(0,len(targets),8):
  part=adata[np.isin(labels,[cfg.control]+targets[lo:lo+8])].copy()
  # Whole-half NMAE anchors reuse the same exact target-complete DE pieces.
  frame=compute(part,cfg=cfg,fp=None,store=None,side=side) if len(targets)>8 else ORIGINAL(part,cfg=cfg,fp=None,store=None,side=side)
  frames.append(frame);del part
 # Do not alphabetize features: NMAE's floating-point sums retain official gene order.
 frame=pl.concat(frames);frame.write_parquet(path);return frame

def install():
 receipt=rt.ROOT/'artifacts/bounded_fixture.json'
 assert receipt.exists() and rt.json.loads(receipt.read_text())['status']=='PASS bitwise equality'
 assert rt.sha(__file__)==rt.json.loads(receipt.read_text())['adapter_sha256'],'Adapter changed since fixture'
 from official_scoring import provenance
 assert provenance()['official_files']==rt.json.loads(receipt.read_text())['cell_eval2']['official_files'],'Official scorer changed since bitwise fixture'
 install_unchecked()
def install_unchecked():
 run._compute_de_side=compute;anchor._compute_de_side=compute;ref._compute_de_side=compute
