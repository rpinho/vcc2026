"""What predicts a target's true response magnitude (||Delta_p||, own gene excluded) on held-out contexts?
Features are label-free (sources + destination controls). Leave-one-context-out ridge on log norms."""
import sys, numpy as np, pandas as pd
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import banks as B, recipes as R, metrics as M
from run_selector import context_panel, predict
from scipy.stats import spearmanr

all_banks = B.discover(); rows = []
for ctx in ['HCT116','HEK293T','CD4_Rest','CD4_Stim8hr','CD4_Stim48hr','H1','K562_GWPS','HepG2','Jurkat','RPE1','K562']:
    bank = all_banks[ctx]; genes = bank.genes; targets = context_panel(bank, all_banks, 300, 0)
    truth = np.nan_to_num(bank.rows(targets)); ctrl = np.nanmean(bank.rows(targets, 'control'), 0); ctrl /= ctrl.sum()
    entries = R.centre(R.gather(all_banks, targets, bank.family), targets, genes, 'source')
    pred, nsrc, nfam, cov = R.pool(entries, targets, genes, 'equal')
    agree = R.agreement(entries, pred, targets, genes)
    m = M.compute(pred, truth, targets, genes, bank.measured)
    gi = {g: i for i, g in enumerate(genes)}
    own = np.array([gi.get(t, -1) for t in targets])
    # per-source norms (max and mean over sources), destination control expression of the target
    src_norms = np.full((len(targets), len(entries)), np.nan)
    for k, (b, r, has) in enumerate(entries):
        rr = np.where(has[:, None] & b.measured[None, :], r, 0.0).astype(np.float64)
        for i in np.flatnonzero(has):
            if own[i] >= 0: rr[i, own[i]] = 0
        src_norms[has, k] = np.linalg.norm(rr[has], axis=1)
    dot = np.einsum('ij,ij->i', np.nan_to_num(pred.astype(np.float64)), truth.astype(np.float64))
    oracle = dot / np.maximum((pred.astype(np.float64) ** 2).sum(1), 1e-12)
    for i, t in enumerate(targets):
        rows.append(dict(context=ctx, target=t, true_norm=m['true_norm'][i], pred_norm=m['pred_norm'][i], oracle_amp=oracle[i], cos=m['per_target_cos'][i],
                         nsrc=nsrc[i], agree=agree[i], src_norm_max=np.nanmax(src_norms[i]), src_norm_mean=np.nanmean(src_norms[i]),
                         ctrl_cpm=ctrl[own[i]] * 1e6 if own[i] >= 0 else np.nan, n_cells=bank.n_cells[t]))
df = pd.DataFrame(rows); df.to_csv('artifacts/magnitude_features.csv', index=False)
df['log_true'] = np.log(df.true_norm + 1e-6); df['log_pred'] = np.log(df.pred_norm + 1e-6); df['log_srcmax'] = np.log(df.src_norm_max + 1e-6)
df['log_ctrl'] = np.log1p(df.ctrl_cpm.fillna(0)); df['agree_f'] = df.agree.fillna(0.3)
print('Spearman with true_norm per context:')
for c, d in df.groupby('context'):
    print(f"{c:12s} pred_norm {spearmanr(d.pred_norm, d.true_norm).correlation:.3f} src_max {spearmanr(d.src_norm_max, d.true_norm).correlation:.3f} "
          f"nsrc {spearmanr(d.nsrc, d.true_norm).correlation:.3f} agree {spearmanr(d.agree_f, d.true_norm).correlation:.3f} ctrl_cpm {spearmanr(d.log_ctrl, d.true_norm).correlation:.3f} "
          f"| oracle_amp vs pred_norm {spearmanr(d.pred_norm, d.oracle_amp).correlation:.3f} vs agree {spearmanr(d.agree_f, d.oracle_amp).correlation:.3f} vs ctrl {spearmanr(d.log_ctrl, d.oracle_amp).correlation:.3f}")
# leave-one-context-out ridge for log(true_norm) and for the oracle amplitude factor (clipped >= 0)
feats = ['log_pred', 'log_srcmax', 'nsrc', 'agree_f', 'log_ctrl']
from numpy.linalg import solve
def ridge(X, y, lam=1.0):
    Xb = np.c_[X, np.ones(len(X))]; A = Xb.T @ Xb + lam * np.eye(Xb.shape[1]); A[-1, -1] -= lam
    return solve(A, Xb.T @ y)
out = []
for c in df.context.unique():
    tr = df[df.context != c]; te = df[df.context == c]
    w = ridge(tr[feats].to_numpy(float), tr.log_true.to_numpy(float))
    yhat = np.c_[te[feats].to_numpy(float), np.ones(len(te))] @ w
    # rho_bar gain if we rescale each target's prediction by exp(yhat)/pred_norm (predicted norm) vs no rescale: a_p ∝ predicted d_p * cos
    out.append(dict(context=c, spearman_pred=spearmanr(yhat, te.log_true).correlation, spearman_base=spearmanr(te.log_pred, te.log_true).correlation, coef=np.round(w, 3).tolist()))
print(pd.DataFrame(out).to_string())
