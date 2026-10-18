"""Diagnostic only: how many real DE genes per target (depth-20k H1), how shared across targets."""
import runtime, numpy as np, anndata as ad, json
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import evaluate_v5 as E
ctrl, real, anchors = E.load_fold('depth20k')
labels = real.obs.target_gene.astype(str).to_numpy(); targets = anchors['targets']
X = real.X.tocsr().astype(np.float64); d = np.asarray(X.sum(1)).ravel()
N = (X.multiply(1e6 / d[:, None])).tocsc()
nt = labels == 'non-targeting'
keep = np.asarray((X[nt] > 0).mean(0)).ravel() > 0.01
genes = np.flatnonzero(keep)
C = N[nt][:, genes].toarray()
sig = {}
for t in targets:
    T = N[labels == t][:, genes].toarray()
    p = mannwhitneyu(T, C, axis=0).pvalue
    q = multipletests(np.nan_to_num(p, nan=1), method='fdr_bh')[1]
    s = np.sign(T.mean(0) - C.mean(0))
    sig[t] = {int(g): int(sg) for g, sg, qq in zip(genes, s, q) if qq < 0.05}
    print(t, len(sig[t]), flush=True)
allg = {}
for t in sig:
    for g, s in sig[t].items():
        allg.setdefault(g, []).append(s)
shared = {g: v for g, v in allg.items() if len(v) >= 6}
cons = np.mean([abs(np.mean(v)) for v in shared.values()]) if shared else None
print('genes DE in >=6/11 targets', len(shared), 'mean |sign agreement|', cons)
