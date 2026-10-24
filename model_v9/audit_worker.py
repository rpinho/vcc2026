"""Independent emitted-count checks; no imports from the emitter."""
import runtime,sys,hashlib
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.stats import spearmanr
from cell_eval2.moments import jackknife_correction
sys.path.append(str(Path(__file__).resolve().parents[1]/'validation'))
import gates_core
sys.path.append(str(Path(__file__).resolve().parents[1]/'model_v7'))
from verify_counts import moment_certificate
from de import Reference
REFS={}
def audit(x,original,c,confidence,own,selected,baseline=None):
 runtime.guard()
 if c not in REFS:REFS[c]=Reference(c)
 ref=REFS[c];T=x.sum(0,dtype=np.int64);d=x.sum(1,dtype=np.int64);od=original.sum(1,dtype=np.int64)
 assert np.array_equal(T,original.sum(0,dtype=np.int64))
 assert np.array_equal(d,od)
 assert np.array_equal(x[:,~ref.keep],original[:,~ref.keep])
 assert np.array_equal(x[:,own],original[:,own])
 mean0=(original/(od[:,None]/1e6)).mean(0) if baseline is None else baseline['mean0'];mean=(x/(d[:,None]/1e6)).mean(0)
 z,p,q=ref.test(x);sig=q<.05;sig[own]=False
 jk=float(jackknife_correction(sparse.csr_matrix(x),np.zeros(400,np.intp),1,50000.,chunk=512)[0])
 jk0=float(jackknife_correction(sparse.csr_matrix(original),np.zeros(400,np.intp),1,50000.,chunk=512)[0]) if baseline is None else float(baseline['jk0'])
 independent=gates_core.jackknife_cp(x,T,T.sum(),chunk=512)
 assert np.isclose(jk,independent,rtol=2e-8,atol=1e-7)
 cert=moment_certificate(x,d,mean0/1e6,T/T.sum())
 selected_mask=np.zeros(len(T),bool);selected_mask[selected]=True
 lfc=np.log2((mean+1e-9)/(ref.mean+1e-9));lfc0=np.log2((mean0+1e-9)/(ref.mean+1e-9))
 rho=float(spearmanr(confidence[sig],-np.log10(np.maximum(p[sig],1e-300))).statistic) if sig.sum()>2 else 0.
 rec=dict(n_calls=int(sig.sum()),selected_calls=len(selected),unintended_calls=int(np.sum(sig&~selected_mask)),confidence_rho=rho,jk=jk,v5b_jk=jk0,jk_ratio=jk/jk0,mean_CPM_L1=float(abs(mean-mean0).sum()/1e6),LFC_change_median=float(np.median(abs(lfc[ref.keep]-lfc0[ref.keep]))),LFC_change_max=float(np.max(abs(lfc[ref.keep]-lfc0[ref.keep]))),mean_certificate=cert,totals_exact=True,depths_exact=True,nonDE_exact=True,own_target_exact=True,changed_entries=int(np.count_nonzero(x!=original)),counts_sha256=hashlib.sha256(x.astype('<i4').tobytes(order='C')).hexdigest())
 rec['called_direction_agreement_with_v5b']=float(np.mean(np.sign(lfc[sig])==np.sign(lfc0[sig]))) if sig.any() else None
 if baseline is None:_,p0,q0=ref.test(original)
 else:p0,q0=baseline['p0'],baseline['q0']
 return rec,dict(p=p,q=q,lfc=lfc,lfc0=lfc0,mean=mean,mean0=mean0,p0=p0,q0=q0,jk0=np.array(jk0))
