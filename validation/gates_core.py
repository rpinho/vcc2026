
"""Independent Layer-0 gates G1-G4 for VCC2026 candidates.

Mirrors the vendored cell_eval2 0.18.0 kernels by hand (not imported) so the numbers are an
INDEPENDENT check of model_v8's own gate implementation:
  - bulk comparator            : log1p(bulk_target_sum * P_pg / S_p)          [norm/prep]
  - jackknife C_p              : moments.jackknife_correction / _loo_bulk
  - across-target budget B     : metrics.delta._across_pert_budget
  - row weights w_p            : metrics.delta._row_weights
  - reference-only CPM gate    : de_compute._ref_cpm_from_cells  (control mean CPM > 5)
  - DE                         : tied-asymptotic Mann-Whitney, 400 pred vs 18,400 real controls,
                                 BH per target (fdr_scope=per_pert), target gene excluded
"""
import numpy as np, h5py, json, hashlib, time, sys, os
from scipy.stats import norm, spearmanr

TS = 50_000.0          # bulk_target_sum
CPM_T = 1e6            # de target_sum
ALPHA = 0.05
GATE = 5.0             # filter_gene_min_cpm_cell


def group_bulk(P, S):
    """log1p(TS * P_g / S) -- the expression comparator. P:[P,G] float64, S:[P]."""
    out = np.zeros_like(P, dtype=np.float64)
    ok = S > 0
    out[ok] = np.log1p(TS * P[ok] / S[ok, None])
    return out


def jackknife_cp(Y, P_g, S_p, chunk=256):
    """Delete-1 jackknife C_p = (n-1)/n * sum_g sum_i (v_ig - vbar_g)^2, v from _loo_bulk."""
    n, G = Y.shape
    if n < 2:
        return 0.0
    lib = Y.sum(axis=1, dtype=np.float64)
    r = S_p - lib
    s1 = np.zeros(G); s2 = np.zeros(G)
    for a in range(0, n, chunk):
        b = min(a + chunk, n)
        num = P_g[None, :] - Y[a:b]
        rr = r[a:b, None]
        v = np.where(rr > 0, np.log1p(TS * num / np.where(rr > 0, rr, 1.0)), 0.0)
        s1 += v.sum(axis=0); s2 += (v * v).sum(axis=0)
    ss = s2 - s1 * s1 / n                      # sum_i (v - vbar)^2 per gene
    return float((n - 1) / n * np.maximum(ss, 0.0).sum())


def row_weights(n_rows, n_genes, cols):
    nr = np.full(n_rows, float(n_genes))
    for k in cols:
        nr[k] = n_genes - 1
    return np.divide(1.0, nr, out=np.zeros_like(nr), where=nr > 0)


def across_pert_budget(pred_sel, cols, weights):
    """Verbatim port of metrics.delta._across_pert_budget."""
    n_rows, n_genes = pred_sel.shape
    w = np.asarray(weights, dtype=np.float64)
    dev = pred_sel - pred_sel.mean(axis=0)
    s1 = np.einsum("ij,i->j", dev, w)
    s2 = np.einsum("ij,ij,i->j", dev, dev, w)
    wsum = np.full(n_genes, float(w.sum())); cnt = np.full(n_genes, float(n_rows))
    for k, col in cols.items():
        v = dev[k, col]
        s1[col] -= w[k] * v; s2[col] -= w[k] * v * v
        wsum[col] -= w[k];   cnt[col]  -= 1.0
    ok = (cnt >= 2.0) & (wsum > 0.0)
    if not ok.any():
        return 0.0
    ss = s2[ok] - s1[ok] ** 2 / wsum[ok]
    return float(np.maximum(ss, 0.0).sum())


def bh(p):
    """Benjamini-Hochberg adjusted p, ascending-sorted step-up with cumulative min."""
    m = p.size
    o = np.argsort(p, kind="stable")
    q = p[o] * m / np.arange(1, m + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(m); out[o] = np.minimum(q, 1.0)
    return out


def ctrl_tie_base(ec_vals, ec_offs, ec_zeros):
    """sum_v (c_v^3 - c_v) over the CONTROL's own distinct values, per universe gene."""
    nU = ec_offs.size - 1
    base = np.zeros(nU)
    for j in range(nU):
        a, b = ec_offs[j], ec_offs[j + 1]
        t = [float(ec_zeros[j])] if ec_zeros[j] > 0 else []
        if b > a:
            cv = ec_vals[a:b]
            runs = np.diff(np.flatnonzero(np.concatenate(([True], cv[1:] != cv[:-1], [True])))).astype(float)
            t.extend(runs.tolist())
        tt = np.asarray(t)
        base[j] = float((tt ** 3 - tt).sum())
    return base


def mwu_vs_controls(pred_cpm, ec_vals, ec_offs, ec_zeros, nC, tie_base):
    """Tied-asymptotic two-sided Mann-Whitney of each gene's predicted CPMs against the
    context's control CPMs. Exact U with midranks; sigma uses the COMBINED sample's full tie
    structure (control tie base, adjusted for the blocks predicted values join).
    Validated against scipy.stats.mannwhitneyu(method='asymptotic') to rtol 1e-12.
    """
    nU = ec_offs.size - 1
    n1 = pred_cpm.shape[0]; n2 = nC; N = n1 + n2
    z = np.zeros(nU); pv = np.ones(nU)
    for j in range(nU):
        a, b = ec_offs[j], ec_offs[j + 1]
        cv = ec_vals[a:b]; nz = ec_zeros[j]
        q = pred_cpm[:, j]
        uq, cnt = np.unique(q, return_counts=True)
        lo = np.empty(uq.size); eq = np.empty(uq.size)
        pos = uq > 0
        li = np.searchsorted(cv, uq[pos], "left"); ri = np.searchsorted(cv, uq[pos], "right")
        lo[pos] = nz + li; eq[pos] = ri - li
        lo[~pos] = 0.0;    eq[~pos] = nz
        U = float((cnt * (lo + 0.5 * eq)).sum())
        merged = eq + cnt
        T = tie_base[j] + float(((merged ** 3 - merged) - (eq ** 3 - eq)).sum())
        var = n1 * n2 / 12.0 * ((N + 1) - T / (N * (N - 1.0)))
        if var <= 0:
            continue
        zz = (U - n1 * n2 / 2.0) / np.sqrt(var)
        z[j] = zz; pv[j] = 2.0 * norm.sf(abs(zz))
    return z, pv
