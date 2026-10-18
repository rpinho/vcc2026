"""Cheap exact official expression-anchor preflight before costly DE calibration."""
import runtime as rt
import sys,time
from pathlib import Path
import official_scoring as s
from cell_eval2.anchor import compute_replicate_anchor
from cell_eval2.catalog import resolve_metrics,CATALOG
from benchmark import setup
name=sys.argv[1];m,ctrl,real=setup(name);del ctrl;names=[n for n in resolve_metrics(s.config().metrics,version='v2')[0] if CATALOG[n].kind=='anndata'];start=time.time();cfg=s.config(metrics=names)
try:
 splits,frame=compute_replicate_anchor(real,config=cfg,n_splits=5);frame.write_csv(f'artifacts/{name}_expression_preflight.csv');rt.save(f'artifacts/{name}_expression_preflight.json',dict(status='PASS expression anchor only; DE anchor not checked',seconds=time.time()-start,source_sha256=s.matrix_id(real),anchor=frame.to_dicts(),official_version=__import__('cell_eval2').__version__))
except ValueError as e:
 rt.save(f'artifacts/{name}_expression_preflight.json',dict(status='INVALID expression anchor: full six-metric benchmark impossible at this reference size',error=str(e),seconds=time.time()-start,source_sha256=s.matrix_id(real)));print('INVALID',name,str(e),flush=True)
