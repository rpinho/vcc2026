"""Response-blind control similarity and task feasibility; bounded backed reads."""
import os
for key in ('OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','OMP_NUM_THREADS','POLARS_MAX_THREADS'):
    os.environ[key] = '1'
import hashlib
import json
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import spearmanr

P = Path(__file__).resolve().parent
V = P.parent


def summary(path, role=None):
    a = ad.read_h5ad(path, backed='r')
    col = 'target_gene' if 'target_gene' in a.obs else 'gene'
    labels = a.obs[col].astype(str)
    controls = labels.str.contains('non-targeting', regex=False).to_numpy()
    if role:
        controls = a.obs.evaluation_role.astype(str).isin(['input_control','reference_control']).to_numpy()
    rows = np.flatnonzero(controls)
    # Uniform deterministic sample; identical size prevents a cell-count similarity advantage.
    if len(rows) > 3200:
        rows = np.sort(np.random.default_rng(20261018).choice(rows,3200,replace=False))
    genes = a.var['gene_name'].astype(str).to_numpy() if 'gene_name' in a.var else a.var_names.to_numpy(str)
    totals = np.zeros(len(genes)); depth = []; expressed = np.zeros(len(genes))
    for lo in range(0,len(rows),400):
        x = sparse.csr_matrix(a.X[rows[lo:lo+400]], dtype=float)
        d = np.asarray(x.sum(1)).ravel(); good = d > 0; x = x[good]; d=d[good]
        norm = sparse.diags(1e6/d) @ x
        totals += np.asarray(norm.sum(0)).ravel(); expressed += np.asarray((x>0).sum(0)).ravel(); depth.extend(d.tolist())
    counts = labels[~controls].value_counts()
    ids = {}
    for key in ('ntc_id','guide_id','sgID_AB','batch','gem_group'):
        if key in a.obs:
            ids[key] = a.obs.loc[controls,key].astype(str).value_counts().to_dict()
    s = pd.Series(totals/max(len(depth),1),index=genes).groupby(level=0).sum()
    info = dict(path=str(path.relative_to(V)),shape=list(a.shape),controls=int(controls.sum()),sampled=len(depth),
                median_umi=float(np.median(depth)),umi_quantiles=np.quantile(depth,[.1,.5,.9]).tolist(),
                n_targets=int(len(counts)),targets_ge400=int((counts>=400).sum()),targets_ge800=int((counts>=800).sum()),
                control_ids=ids)
    a.file.close()
    return s,info


def main():
    profiles={}; metadata={}
    paths={**{c:V/f'data/context_{c}.h5ad' for c in 'ABC'},
           **{d:P/f'tasks/{d}.h5ad' for d in ('HEK293T','HCT116','K562','CD4')},
           'H1':V/'rsi/artifacts/tasks/H1_rich_search.h5ad',
           **{d:V/f'public_data/{d}.h5ad' for d in ('HepG2','Jurkat','RPE1')}}
    for name,path in paths.items():
        profiles[name],metadata[name]=summary(path,role=True if name=='H1' else None)
        print(name, {k:v for k,v in metadata[name].items() if k!='control_ids'},flush=True)
    common=sorted(set.intersection(*(set(s.index) for s in profiles.values())))
    # All comparisons use the same measured genes and normalized fractions on that axis.
    vectors={k:np.log1p(s.reindex(common).to_numpy()/s.reindex(common).sum()*1e6) for k,s in profiles.items()}
    candidates=[k for k in profiles if k not in 'ABC']
    feasible=('HEK293T','HCT116','K562','CD4','H1')
    similarity={}; selected={}
    for c in 'ABC':
        similarity[c]={d:dict(pearson=float(np.corrcoef(vectors[c],vectors[d])[0,1]),
                              spearman=float(spearmanr(vectors[c],vectors[d]).statistic),
                              cosine=float(vectors[c]@vectors[d]/np.linalg.norm(vectors[c])/np.linalg.norm(vectors[d]))) for d in candidates}
        selected[c]=max(feasible,key=lambda d:similarity[c][d]['pearson'])
    weights={d:list(selected.values()).count(d)/3 for d in sorted(set(selected.values()))}
    out=dict(rule='maximum Pearson log1p mean cell CPM on common measured genes; equal contexts, no score data',
             common_genes=len(common),common_gene_sha256=hashlib.sha256('\0'.join(common).encode()).hexdigest(),
             metadata=metadata,similarity=similarity,selected=selected,destination_weights=weights)
    (P/'artifacts/control_audit.json').write_text(json.dumps(out,indent=2))
    print('SELECTED',selected,weights,flush=True)


if __name__=='__main__':
    try:
        os.nice(10)
    except OSError:
        pass
    main()
