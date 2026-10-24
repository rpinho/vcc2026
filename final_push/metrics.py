"""Selector metrics in log1p(CP50K) group-sum bulk space (the space of pds_cosine and expr_mse).

rho_bar: weighted Frobenius cosine between predicted and measured responses over the panel, own target gene excluded,
official 1/|G_p| row weights (cell-eval2 delta._row_weights). At break-even credit (B = claim) the official MSE is
u = 1 - 2 * rho_bar * |A|/|Delta| (consult_1021/FABLE.md section 1), so rho_bar orders direction/magnitude changes.
pds: the official pds_cosine formula (all panel target genes removed, tie-aware rank).
"""
import numpy as np
from scipy.stats import spearmanr


def own_cols(targets, genes):
    gi = {g: i for i, g in enumerate(genes)}
    return np.array([gi.get(t, -1) for t in targets])


def compute(pred, true, targets, genes, measured):
    T, G = pred.shape
    own = own_cols(targets, genes)
    M = np.broadcast_to(np.asarray(measured, bool)[None, :], (T, G)).copy()
    for i, c in enumerate(own):
        if c >= 0:
            M[i, c] = False
    Gp = M.sum(1).astype(np.float64); w = np.divide(1.0, Gp, out=np.zeros(T), where=Gp > 0)
    P = np.where(M, pred, 0.0).astype(np.float64); Y = np.where(M, true, 0.0).astype(np.float64)
    dot = (P * Y).sum(1); pp = (P * P).sum(1); yy = (Y * Y).sum(1)
    wpp = max((w * pp).sum(), 1e-300); wyy = max((w * yy).sum(), 1e-300)
    rho_bar = float((w * dot).sum() / np.sqrt(wpp * wyy))
    cos_p = np.divide(dot, np.sqrt(pp * yy), out=np.zeros(T), where=(pp > 0) & (yy > 0))
    opt = np.divide(dot ** 2, pp, out=np.zeros(T), where=pp > 0)       # oracle per-target amplitude
    rho_opt = float(np.sqrt(max((w * opt).sum(), 0) / wyy))
    mag = float(spearmanr(np.sqrt(pp), np.sqrt(yy)).correlation) if T > 3 else float('nan')
    keep = np.asarray(measured, bool).copy()
    for c in own:
        if c >= 0:
            keep[c] = False
    Pk = pred[:, keep].astype(np.float64); Yk = true[:, keep].astype(np.float64)
    pn = np.linalg.norm(Pk, axis=1); yn = np.linalg.norm(Yk, axis=1)
    D = 1 - (Pk @ Yk.T) / np.maximum(pn[:, None] * yn[None, :], 1e-300)
    D[(pn == 0)[:, None] | (yn == 0)[None, :]] = 1.0
    pds = np.empty(T)
    for i in range(T):
        d = D[i]; dii = d[i]
        pds[i] = 1 - ((d < dii).sum() + 0.5 * ((d == dii).sum() - 1)) / max(T - 1, 1)
    e = w * yy; share = float(np.sort(e)[::-1][:max(1, T // 10)].sum() / max(e.sum(), 1e-300))
    # shared component of the truth: energy explained by assigning the panel-mean true response to every target
    ybar = (w[:, None] * Y).sum(0) / max(w.sum(), 1e-300); Yc = np.where(M, Y - ybar[None, :], 0.0)
    mean_share = float(1 - (w * (Yc * Yc).sum(1)).sum() / wyy)
    # same for the prediction
    pbar = (w[:, None] * P).sum(0) / max(w.sum(), 1e-300); Pc = np.where(M, P - pbar[None, :], 0.0)
    pred_mean_share = float(1 - (w * (Pc * Pc).sum(1)).sum() / wpp)
    # rho_bar of the CENTRED matrices (what is left after both shared components are removed)
    rho_centred = float((w * (Pc * Yc).sum(1)).sum() / np.sqrt(max((w * (Pc * Pc).sum(1)).sum(), 1e-300) * max((w * (Yc * Yc).sum(1)).sum(), 1e-300)))
    return dict(rho_bar=rho_bar, rho_bar_oracle_amp=rho_opt, cos_mean=float(cos_p.mean()), cos_median=float(np.median(cos_p)),
                pds=float(pds.mean()), magnitude_spearman=mag, pred_over_true_norm=float(np.sqrt(wpp / wyy)),
                top_decile_energy_share=share, n_targets=int(T), truth_mean_share=mean_share, pred_mean_share=pred_mean_share, rho_bar_centred=rho_centred,
                est_scaled_mse_fullcredit=float(np.clip((1.41 * rho_bar - 0.01) / 0.95, 0, 1)),
                est_scaled_pds=float((pds.mean() - 0.5) / 0.45), per_target_cos=cos_p, per_target_pds=pds,
                pred_norm=np.sqrt(pp), true_norm=np.sqrt(yy))
