import runtime_local as j
j.guard(launch=True)
assert j.os.getpriority(j.os.PRIO_PROCESS,0)>=10
import numpy as np, anndata as ad
from scipy import sparse
from cell_eval2.moments import jackknife_correction

def moments(x):
    x=sparse.csr_matrix(x,dtype=np.float64)
    sums=np.asarray(x.sum(0)).ravel()
    mean=np.log1p(50000*sums/sums.sum())
    jk=jackknife_correction(x,np.zeros(x.shape[0],np.intp),1,50000.,chunk=32)[0]
    return mean,float(jk)

if __name__=='__main__':
    path=j.VCC/'rsi/artifacts/tasks/H1_rich_search.h5ad'
    a=ad.read_h5ad(path,backed='r')
    labels=a.obs.target_gene.astype(str).to_numpy()
    inputs=a.obs.evaluation_role.astype(str).to_numpy()=='input_control'
    targets=np.array(sorted(set(labels[~inputs])-{'non-targeting'}))
    expected=j.json.loads((j.VCC/'h2h2_1026/results/official_B0.json').read_text())['targets']
    assert targets.tolist()==sorted(expected) and len(targets)==79
    order=np.r_[targets,['non-targeting']]
    means=[];jk=[];counts=[]
    for i,t in enumerate(order):
        j.guard()
        ix=np.flatnonzero((~inputs)&(labels==t)); x=a.X[ix]
        m,c=moments(x);means.append(m);jk.append(c);counts.append(len(ix))
        if i%10==0:j.status(f'Reference moments {i+1}/80; {t}, n={len(ix)}')
    ctrl=sparse.csr_matrix(a.X[np.flatnonzero(inputs)])
    control_mean,control_jk=moments(ctrl)
    # Prespecified Monte Carlo average, with replacement: population null at n_p.
    rng=np.random.default_rng(20261026)
    resamples={}
    for n in sorted(set(counts[:-1])):
        draws=[]
        for k in range(16):
            j.guard()
            _,c=moments(ctrl[rng.choice(ctrl.shape[0],int(n),replace=True)])
            draws.append(c)
        resamples[str(n)]=draws
    np.savez_compressed(j.ROOT/'reference.npz',targets=order,genes=a.var_names.to_numpy(str),means=means,
        jk=jk,counts=counts,input_mean=control_mean,input_jk=control_jk,
        null_jk=np.array([np.mean(resamples[str(n)]) for n in counts[:-1]]))
    j.save('reference.json',dict(task=str(path),sha256=j.sha(path),targets=targets.tolist(),
        null_source='input controls only; bootstrap with replacement at each real n_p',
        bootstrap_seed=20261026,bootstrap_replicates=16,resamples=resamples,
        input_n=ctrl.shape[0],free_gib=j.guard()))
    a.file.close();j.status('Reference and label-free control-resample corrections saved.')
