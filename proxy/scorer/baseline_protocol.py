"""Arc-documented mean-response comparator; control-mean is sensitivity only.
Profile construction and every scored metric are official cell-eval2 functions.
The public H1 harness README disagrees with Arc; that ambiguity is recorded, not hidden.
"""
import runtime as rt
import numpy as np,anndata as ad
from cell_eval2.baseline import generic_response_profile

def mean_response(real):
 a=ad.AnnData(real.X,obs=real.obs.copy(),var=real.var.copy());labels=a.obs.target_gene.astype(str).to_numpy();column='target_gene'
 if 'guide_id' in a.obs:
  column='baseline_construct';a.obs[column]=np.where(labels=='non-targeting','non-targeting',a.obs.guide_id.astype(str).to_numpy())
 profile=generic_response_profile(a,pert_col=column,control='non-targeting',exclude_target_gene=False)
 return profile.values,dict(profile='Official generic_response_profile; equal construct means when guide IDs exist, otherwise equal target means',n_constructs=profile.n_perturbations,exclude_target_gene=False,emit='float32 identical fractional cells, 400 per target',qualification='Uses already selected reference cells. HCT116 guide calls pass source filter; guide-level knockdown-filter parameters for the hidden official baseline are not published. Exact official bundle reproduction is unverified.',builder_sha256=rt.sha(__file__))
