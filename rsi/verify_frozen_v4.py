"""Replay all canonical H1 v4a counts after source-isolation/seed-plumbing edits."""
import runtime as rt
from pathlib import Path
import hashlib
from scipy import sparse
from benchmark import setup
from recipes import predict
m,ctrl,real=setup('H1_search');del real;h=hashlib.sha256();count=0
for t,x in predict(ctrl,m['targets'],m['context'],'v4a',42):
 x=sparse.csr_matrix(x);x.eliminate_zeros()
 for value in [x.data,x.indices,x.indptr]:h.update(value.tobytes())
 h.update(t.encode());count+=1
 if count%25==0:print('Canonical replay',count,flush=True)
old=rt.json.loads(Path('artifacts/scores/H1_search/v4a_s42.json').read_text());assert h.hexdigest()==old['prediction_sha256'];rt.save('artifacts/frozen_v4_fixture.json',dict(status='PASS all canonical prediction counts bitwise unchanged',targets=count,prediction_sha256=h.hexdigest(),recipe_sha256=rt.sha('recipes.py'),scope='Old H1 v4a seed42 count stream versus current isolated-source and explicit alternate-seed plumbing; no metric recomputation'));print('PASS canonical v4a replay',count)
