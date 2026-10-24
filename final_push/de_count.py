"""How many reference DE calls do panel-like targets get under the official rules (400 cells, ~20k UMI median,
Wilcoxon rank-sum two-sided normal approximation with tie correction, BH within perturbation, 5-CPM control filter)?
Uses only the RSI SEARCH tasks (never the sealed HOLDOUT): H1_rich_search (79 x 400 cells) and
CD4_depth192_search (79 x 192 cells, a lower-power caveat). Writes artifacts/de_count_<task>.csv."""
import sys, h5py, numpy as np, pandas as pd
from pathlib import Path
from scipy import sparse
from scipy.stats import mannwhitneyu, false_discovery_control

HERE = Path(__file__).resolve().parent; VCC = HERE.parent
TASKS = VCC / 'rsi/artifacts/tasks'


def load(task):
    with h5py.File(TASKS / f'{task}.h5ad') as h:
        X = h['X']; shape = tuple(h['X'].attrs['shape'])
        M = sparse.csr_matrix((X['data'][...], X['indices'][...], X['indptr'][...]), shape=shape)
        obs = h['obs']
        role = obs['evaluation_role']['categories'][...].astype(str)[obs['evaluation_role']['codes'][...]]
        tg = obs['target_gene']['categories'][...].astype(str)[obs['target_gene']['codes'][...]]
        vi = h['var']['_index']
        if isinstance(vi, h5py.Group):      # nullable-string encoding: {mask, values}
            assert not np.any(vi['mask'][...]); genes = vi['values'][...].astype(str)
        else:
            genes = vi[...].astype(str)
    return M, role, tg, genes


def thin(M, target_median, seed=0):
    d = np.asarray(M.sum(1)).ravel(); med = np.median(d)
    if med <= target_median:
        return M, float(med), 1.0
    f = target_median / med; rng = np.random.default_rng(seed); M = M.tocoo()
    data = rng.binomial(M.data.astype(np.int64), f)
    M = sparse.csr_matrix((data, (M.row, M.col)), shape=M.shape); M.eliminate_zeros()
    return M, float(np.median(np.asarray(M.sum(1)).ravel())), f


def main(task, target_median=20000):
    M, role, tg, genes = load(task)
    M, med, f = thin(M, target_median)
    print(task, 'median depth after thinning', med, 'fraction', round(f, 3), flush=True)
    ctrl = M[role == 'reference_control']
    depth = np.asarray(ctrl.sum(1)).ravel(); C = sparse.diags(1e6 / depth) @ ctrl
    cmean = np.asarray(C.mean(0)).ravel(); keep = cmean > 5.0
    Cd = C[:, keep].toarray(); kept = genes[keep]
    print('tested genes', int(keep.sum()), 'controls', ctrl.shape[0], flush=True)
    rows = []
    for t in np.unique(tg[role == 'perturbation']):
        P = M[(role == 'perturbation') & (tg == t)]; dp = np.asarray(P.sum(1)).ravel()
        Pd = (sparse.diags(1e6 / dp) @ P)[:, keep].toarray()
        p = mannwhitneyu(Pd, Cd, alternative='two-sided', method='asymptotic', use_continuity=False, axis=0).pvalue
        own = kept == t; p = p[~own]
        q = false_discovery_control(p); sig = q < 0.05
        lfc = np.log2((Pd.mean(0) + 1e-9) / (Cd.mean(0) + 1e-9))[~own]
        rows.append(dict(task=task, target=t, cells=int(P.shape[0]), n_real=int(sig.sum()), n_real_finite_lfc=int((sig & np.isfinite(lfc)).sum()),
                         n_q10=int((q < 0.10).sum()), sum_abs_lfc_sig=float(np.abs(lfc[sig]).sum())))
        print(rows[-1], flush=True)
    df = pd.DataFrame(rows); df.to_csv(HERE / 'artifacts' / f'de_count_{task}.csv', index=False)
    n = df.n_real.to_numpy()
    print(task, 'targets', len(n), 'n_real==0:', int((n == 0).sum()), '<10:', int((n < 10).sum()), '<50:', int((n < 50).sum()),
          'median', float(np.median(n)), 'p90', float(np.percentile(n, 90)), flush=True)


if __name__ == '__main__':
    for task in sys.argv[1:] or ['H1_rich_search', 'CD4_depth192_search']:
        main(task)
