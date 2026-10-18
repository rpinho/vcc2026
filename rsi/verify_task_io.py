"""Independent reference checks for bounded in-place CSR compaction."""
import runtime as rt
from pathlib import Path
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from task_io import split
from official_scoring import matrix_id
rng=np.random.default_rng(415);cases=[];file=Path('artifacts/tmp/io_fixture.h5ad')
for mode in ['mixed','prefix','middle','suffix']:
 x=sparse.random(71,23,density=.3,format='csr',random_state=rng);x.data=rng.integers(1,20,len(x.data)).astype(np.int32);mask=np.zeros(71,bool)
 if mode=='mixed':mask[rng.choice(71,19,replace=False)]=True
 elif mode=='prefix':mask[:19]=True
 elif mode=='middle':mask[21:40]=True
 else:mask[-19:]=True
 a=ad.AnnData(x,obs=pd.DataFrame({'target_gene':['x']*71,'evaluation_role':np.where(mask,'input_control','perturbation')},index=[str(i) for i in range(71)]));expected=[a[mask].copy(),a[~mask].copy()];a.write_h5ad(file);actual=split(file)
 for left,right in zip(expected,actual):
  assert matrix_id(left)==matrix_id(right)
 cases.append(mode)
file.unlink();ctrl,real=split('artifacts/tasks/H1_search.h5ad');assert matrix_id(real)==rt.json.loads(Path('artifacts/scores/H1_search/anchors.json').read_text())['real_sha256'];a=ad.read_h5ad('artifacts/tasks/H1_search.h5ad',backed='r');expected=a[a.obs.evaluation_role.astype(str)=='input_control'].to_memory();assert matrix_id(ctrl)==matrix_id(expected);a.file.close();rt.save('artifacts/task_io_fixture.json',dict(status='PASS bitwise CSR and label identities',cases=cases+['actual H1 reference versus previously computed anchor identity','actual H1 controls versus ordinary backed slicing'],implementation_sha256=rt.sha('task_io.py')));print('PASS bounded task I/O')
