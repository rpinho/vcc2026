"""Check the metadata-only baseline wrapper against the previous full-copy wrapper."""
import runtime as rt
from pathlib import Path
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from cell_eval2.baseline import generic_response_profile
from baseline_protocol import mean_response
from official_scoring import matrix_id
def check(real):
 before=matrix_id(real);old=real.copy();labels=old.obs.target_gene.astype(str).to_numpy();column='target_gene'
 if 'guide_id' in old.obs:column='baseline_construct';old.obs[column]=np.where(labels=='non-targeting','non-targeting',old.obs.guide_id.astype(str).to_numpy())
 expected=generic_response_profile(old,pert_col=column,control='non-targeting',exclude_target_gene=False);actual,description=mean_response(real);np.testing.assert_array_equal(expected.values,actual);assert matrix_id(real)==before and description['n_constructs']==expected.n_perturbations
rng=np.random.default_rng(641);labels=np.repeat(['non-targeting','a','b','c'],20);real=ad.AnnData(sparse.csr_matrix(rng.poisson(3,(80,41))),obs=pd.DataFrame({'target_gene':labels,'guide_id':[str(t)+':'+str(i%2) for i,t in enumerate(labels)]},index=[str(i) for i in range(80)]));check(real)
a=ad.read_h5ad('artifacts/tasks/H1_search.h5ad',backed='r');labels=a.obs.target_gene.astype(str).to_numpy();chosen=np.r_[np.flatnonzero(labels=='non-targeting')[:40],*[np.flatnonzero(labels==t)[:30] for t in sorted(set(labels)-{'non-targeting'})[:4]]];real=a[chosen].to_memory();a.file.close();check(real);selection=pd.read_parquet('artifacts/tasks/H1_search_selection.parquet');guides=dict(zip(selection.source_row,selection.guide_id.astype(str)));real.obs['guide_id']=[guides[i] for i in real.obs.source_row];check(real)
rt.save('artifacts/baseline_profile_fixture.json',dict(status='PASS bitwise profile equality and unchanged reference',builder_sha256=rt.sha('baseline_protocol.py'),compatible_legacy_builder_sha256=['a1cb00180772c3e41e1ec1067c2a1ef6c0cde7c7ee3c99e36588c105d80c4754'],cases=['synthetic multi-guide constructs','actual H1 subset, multiple target and control guides'],scope='Only full-copy versus shared-readonly matrix construction changed; same official generic_response_profile call and metadata'));print('PASS baseline profile memory refactor')
