"""Exact v5b emitter, bounded process worker; no changed model logic."""
import runtime, sys
sys.path.append(str(runtime.ROOT.parent/'model_v5'))
import model
import numpy as np
from scipy import sparse
def initialize(q,d):
 global Q,D
 Q,D=q,d
def emit(job):
 b,seed=job
 runtime.guard()
 return sparse.csr_matrix(model.emit(Q,None,b,D,seed,'anchor_single_mean').astype(np.int32))
