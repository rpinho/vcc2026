"""Target-strength predictor: which knockdowns do anything at all in THIS context?

Truth (per fold, for validation only): noise-corrected centred response energy of the held-out screen's own row.
A bank row pooled over n cells has E[energy] = signal + c/n; c is estimated as the 30th percentile of (energy x n)
over that screen's rows (the coordinator's crude floor, kept identical so the numbers are comparable).

Features (all label-free, computable for an unseen context from public sources + that context's controls):
  sim_strength   similarity-weighted mean source signal strength (weight = billm compatibility
                 exp(-mean |log1p(5e4 p_src_ctrl) - log1p(5e4 p_dest_ctrl)|) over measured genes)
  max_strength   max source signal strength
  n_src          how many source screens measure the target
  agreement      mean cosine of each source's centred row to the pooled direction
  depletion      min over sources of n_cells(target) / median n_cells(screen)  (essential genes lose cells)
  ctrl_cpm       the target gene's own expression in the destination's controls (CRISPRi of an unexpressed gene)
  pred_norm      norm of the pooled predicted response

Validation is leave-one-FAMILY-out: coefficients are fitted on the other families' observed strengths and applied
to the held-out family, so nothing from the evaluated context is used to fit.
Writes artifacts/strength_features.csv and artifacts/strength_lofo.json.
"""
import sys, json
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
import banks as B, recipes as R
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
CONTEXTS = ['HCT116', 'HEK293T', 'CD4_Rest', 'CD4_Stim8hr', 'CD4_Stim48hr', 'H1', 'K562_GWPS', 'HepG2', 'Jurkat', 'RPE1', 'K562']
FEATS = ['log_sim_strength', 'log_max_strength', 'n_src', 'agreement', 'log_depletion', 'log_ctrl_cpm', 'log_pred_norm']


def screen_strength(bank, floor_pct=30):
    """Noise-corrected per-target signal strength for every row of a screen, in its own source-wide centred space."""
    targets = bank.targets
    r = bank.rows(targets); mean = bank.source_mean()['mean']
    gi = {g: i for i, g in enumerate(bank.genes)}
    M = bank.measured
    out = {}
    e = np.zeros(len(targets)); n = np.array([bank.n_cells[t] for t in targets], float)
    for i, t in enumerate(targets):
        v = np.where(M, r[i] - mean, 0.0).astype(np.float64)
        j = gi.get(t)
        if j is not None:
            v[j] = 0.0
        e[i] = float(v @ v) / max(M.sum(), 1)
    c = float(np.percentile(e * n, floor_pct))
    signal = np.maximum(e - c / np.maximum(n, 1), 0.0)
    for i, t in enumerate(targets):
        out[t] = dict(energy=e[i], signal=signal[i], n_cells=n[i])
    return out, c, e, signal, n


def main():
    all_banks = B.discover(include_shards=False)
    strengths = {}; medcells = {}
    for name, bank in all_banks.items():
        s, c, e, sig, n = screen_strength(bank)
        strengths[name] = s; medcells[name] = float(np.median(n))
        print(f'{name:12s} floor c={c:.4g} median n={medcells[name]:.0f} silent(signal==0) {float((sig==0).mean()):.3f} '
              f'top10%share {float(np.sort(sig)[::-1][:max(1,len(sig)//10)].sum()/max(sig.sum(),1e-30)):.3f}', flush=True)
    rows = []
    for ctx in CONTEXTS:
        bank = all_banks[ctx]; genes = bank.genes
        from run_selector import context_panel
        targets = context_panel(bank, all_banks, 300, 0)
        dest_ctrl = np.nanmean(bank.rows(targets, 'control'), 0).astype(np.float64); dest_ctrl /= dest_ctrl.sum()
        dst = np.log1p(50000 * dest_ctrl)
        entries = R.centre(R.gather(all_banks, targets, bank.family), targets, genes, 'source')
        pred, nsrc, nfam, cov = R.pool(entries, targets, genes, 'equal')
        agree = R.agreement(entries, pred, targets, genes)
        gi = {g: i for i, g in enumerate(genes)}
        sim = {}
        for b, r, has in entries:
            sc = b.source_mean()['control']
            m = b.measured
            sim[b.source] = float(np.exp(-np.mean(np.abs(dst[m] - np.log1p(50000 * sc[m])))))
        pn = np.linalg.norm(np.nan_to_num(pred), axis=1)
        for i, t in enumerate(targets):
            ws = []; vs = []; dep = []
            for b, r, has in entries:
                if not has[i]:
                    continue
                st = strengths[b.source].get(t)
                if st is None:
                    continue
                ws.append(sim[b.source]); vs.append(st['signal']); dep.append(st['n_cells'] / max(medcells[b.source], 1))
            if not vs:
                continue
            ws = np.array(ws); vs = np.array(vs)
            j = gi.get(t)
            rows.append(dict(context=ctx, family=bank.family, target=t,
                             sim_strength=float((ws * vs).sum() / ws.sum()), max_strength=float(vs.max()),
                             n_src=int(len(vs)), agreement=float(agree[i]) if np.isfinite(agree[i]) else np.nan,
                             depletion=float(np.min(dep)), ctrl_cpm=float(dest_ctrl[j] * 1e6) if j is not None else np.nan,
                             pred_norm=float(pn[i]),
                             true_signal=float(strengths[ctx][t]['signal']), true_energy=float(strengths[ctx][t]['energy'])))
        print('features', ctx, len(rows), flush=True)
    df = pd.DataFrame(rows)
    df['log_sim_strength'] = np.log(df.sim_strength + 1e-9); df['log_max_strength'] = np.log(df.max_strength + 1e-9)
    df['log_depletion'] = np.log(df.depletion.clip(1e-3)); df['log_ctrl_cpm'] = np.log1p(df.ctrl_cpm.fillna(0))
    df['log_pred_norm'] = np.log(df.pred_norm + 1e-9); df['agreement'] = df.agreement.fillna(df.agreement.median())
    df['log_true'] = np.log(df.true_signal + 1e-9); df['silent'] = (df.true_signal <= 0).astype(int)
    df.to_csv(HERE / 'artifacts/strength_features.csv', index=False)

    def ridge(X, y, lam=1.0):
        Xb = np.c_[X, np.ones(len(X))]; A = Xb.T @ Xb + lam * np.eye(Xb.shape[1]); A[-1, -1] -= lam
        return np.linalg.solve(A, Xb.T @ y)

    out = []
    for fam in sorted(df.family.unique()):
        tr = df[df.family != fam]; te = df[df.family == fam]
        w = ridge(tr[FEATS].to_numpy(float), tr.log_true.to_numpy(float))
        yhat = np.c_[te[FEATS].to_numpy(float), np.ones(len(te))] @ w
        rho = spearmanr(yhat, te.true_signal).correlation
        rho_base = spearmanr(te.log_max_strength, te.true_signal).correlation
        rho_pred = spearmanr(te.log_pred_norm, te.true_signal).correlation
        # does the predictor separate silent from non-silent?
        s = te.silent.to_numpy()
        auc = np.nan
        if 0 < s.sum() < len(s):
            r = pd.Series(yhat).rank().to_numpy()
            auc = float((r[s == 0].mean() - r[s == 1].mean()) / len(s) + 0.5)
        # tier accuracy: predicted bottom third vs actual bottom third of true signal
        k = max(1, len(te) // 3)
        pred_bottom = set(np.argsort(yhat)[:k]); true_bottom = set(np.argsort(te.true_signal.to_numpy())[:k])
        pred_top = set(np.argsort(yhat)[-k:]); true_top = set(np.argsort(te.true_signal.to_numpy())[-k:])
        out.append(dict(family=fam, n=len(te), spearman_predictor=float(rho), spearman_max_strength=float(rho_base),
                        spearman_pred_norm=float(rho_pred), silent_auc=auc, silent_frac=float(s.mean()),
                        bottom_third_recall=len(pred_bottom & true_bottom) / k, top_third_recall=len(pred_top & true_top) / k,
                        coef=dict(zip(FEATS + ['intercept'], np.round(w, 4).tolist()))))
        print(out[-1], flush=True)
    w_all = ridge(df[FEATS].to_numpy(float), df.log_true.to_numpy(float))
    (HERE / 'artifacts/strength_lofo.json').write_text(json.dumps(dict(
        features=FEATS, lofo=out, pooled_coef=dict(zip(FEATS + ['intercept'], np.round(w_all, 4).tolist())),
        note='LOFO coefficients validate the predictor; the pooled coefficients are what a blind final would use'), indent=2))
    d = pd.DataFrame(out)
    print('\nMACRO spearman: predictor %.3f | max_strength %.3f | pred_norm %.3f | silent AUC %.3f | bottom-third recall %.3f | top-third recall %.3f'
          % (d.spearman_predictor.mean(), d.spearman_max_strength.mean(), d.spearman_pred_norm.mean(),
             d.silent_auc.mean(skipna=True), d.bottom_third_recall.mean(), d.top_third_recall.mean()))


if __name__ == '__main__':
    main()
