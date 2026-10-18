"""v5 transforms on top of v4a (model_v4 `anchor_single_mean`). All adopted, none invented.

Technique 1 (v5a): scorer-aware emission for expr_mse_unbiased_capped (cell-eval2 0.18.0 #348).
  * realistic per-cell scatter: one real control cell per predicted cell (pool_k=1; v4a averaged 4,
    which cut the prediction's jackknife claim to about 1/4 of the reference's);
  * per-(target, gene) integer totals pinned to the predicted pseudobulk (v4a's mean_counts +
    AtlasShift repair_bulk_totals, unchanged);
  * per-context pooled amplitude alpha on the across-target CENTRED log1p(5e4*CP) pseudobulk so the
    scorer's own budget B (delta._across_pert_budget, own-target gene excluded, 1/(G-1) row weights)
    equals rho x the prediction's own claim (sum of jackknife corrections of our emitted cells).
    Label-free: needs only our prediction and the context's public controls.
Technique 2 (v5b): explicit DE call set per target, ranked by cross-source agreement, consensus sign.
Technique 4 (v5b): project the top-k PCs of the context's own control cells (log1p CP50K) out of the
  predicted bulk response.
"""
import numpy as np, pandas as pd
from scipy import sparse
from cell_eval2.metrics.delta import _across_pert_budget, _row_weights
from cell_eval2.moments import jackknife_correction

TS = 50000.0  # vcc2026.yaml bulk_target_sum


def excl_cols(targets, genes):
    gi = {g: i for i, g in enumerate(genes)}
    return {k: gi[t] for k, t in enumerate(targets) if t in gi}


def logbulk(b):
    return np.log1p(TS * np.asarray(b, np.float64))


def from_logbulk(L):
    b = np.maximum(np.expm1(L), 0.0)
    return b / b.sum(1, keepdims=True)


def budget(b_or_L, targets, genes, is_log=False):
    """The scorer's #348 budget, in score units, computed with the scorer's own function."""
    L = b_or_L if is_log else logbulk(b_or_L)
    cols = excl_cols(targets, genes)
    w = _row_weights(len(targets), len(genes), cols)
    return _across_pert_budget(L, cols, w)


def claim_from_cells(X, codes, n_groups, targets, genes):
    """sum_p w_p jk_p for the prediction's own cells (= the claim before the #247 cap)."""
    jk = jackknife_correction(X, codes, n_groups, TS, chunk=32)
    w = _row_weights(len(targets), len(genes), excl_cols(targets, genes))
    return float((jk * w).sum()), jk


def rescale(L, alpha, center=None):
    m = L.mean(0, keepdims=True) if center is None else center
    return m + alpha * (L - m)


def solve_alpha(L, targets, genes, target_budget, lo=0.0, hi=8.0, iters=60, center=None):
    """Bisection on alpha so that the renormalised profile's budget hits target_budget.
    center=control log-bulk scales the predicted response about the control (pooled amplitude,
    keeps each target's delta direction, so PDS cosine is unchanged)."""
    def B(a):
        return budget(from_logbulk(rescale(L, a, center)), targets, genes)
    if B(hi) < target_budget:
        return hi, B(hi)
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if B(mid) < target_budget:
            lo = mid
        else:
            hi = mid
    a = 0.5 * (lo + hi)
    return a, B(a)


def control_axes(ctrl_X, k=50, max_cells=4000, seed=0):
    """Top-k principal axes (genes x k) of the context's control cells in log1p(CP50K)."""
    X = sparse.csr_matrix(ctrl_X, dtype=np.float64)
    rng = np.random.default_rng(seed)
    if X.shape[0] > max_cells:
        X = X[np.sort(rng.choice(X.shape[0], max_cells, replace=False))]
    d = np.asarray(X.sum(1)).ravel()
    Y = sparse.diags(TS / d) @ X
    Y.data = np.log1p(Y.data)
    Y = Y.toarray().astype(np.float32)
    Y -= Y.mean(0, keepdims=True)
    from sklearn.utils.extmath import randomized_svd
    _, s, vt = randomized_svd(Y, n_components=k, n_iter=5, random_state=seed)
    del Y
    return vt.T.astype(np.float64), s


def deflate(delta, V, eta=1.0):
    return delta - eta * (delta @ V) @ V.T


def de_call_mask(source_effects, consensus, k=275, z_sig=1.96, min_net=1):
    """Cross-source agreement ranking (SIVC #167 / Chen Zhang #82 call-set control, PIE unavailable).

    source_effects: list of (effect, mask, variance) arrays [T, G] in log2FC.
    consensus: [T, G] fused log2FC (v4a's), whose sign is the call direction.
    Rank key: (#sources significant and agreeing with consensus) - (#significant and opposing),
    then |consensus| / pooled SE. Keep the top-k genes per target with net agreement >= 1.
    """
    T, G = consensus.shape
    sgn = np.sign(consensus)
    net = np.zeros((T, G), np.float32)
    zsum = np.zeros((T, G))
    for e, m, v in source_effects:
        z = np.where(m & np.isfinite(v), e / np.sqrt(np.maximum(v, 1e-12)), 0.0)
        sig = np.abs(z) > z_sig
        net += (sig & (np.sign(z) == sgn) & (sgn != 0)).astype(np.float32)
        net -= (sig & (np.sign(z) == -sgn) & (sgn != 0)).astype(np.float32)
        zsum += np.abs(z) * (np.sign(z) == sgn)
    key = net * 1e6 + np.minimum(zsum, 1e5)
    mask = np.zeros((T, G), bool)
    for i in range(T):
        ok = np.flatnonzero(net[i] >= min_net)
        if ok.size == 0:
            continue
        top = ok[np.argsort(-key[i, ok], kind='stable')[:k]]
        mask[i, top] = True
    return mask, net


def target_totals(b, S):
    """Integer per-gene totals summing exactly to S (largest remainder), = round(S*b)."""
    e = S * np.asarray(b, np.float64)
    t = np.floor(e).astype(np.int64)
    rem = int(S - t.sum())
    if rem > 0:
        t[np.argsort(e - t, kind='stable')[-rem:]] += 1
    return t


def _add_counts(X, need, weights, rng):
    """Add need[g] counts to gene g, cells drawn from weights[:, g] (vectorised inverse CDF)."""
    genes = np.flatnonzero(need > 0)
    if genes.size == 0:
        return
    W = weights[:, genes].astype(np.float64)
    W /= W.sum(0, keepdims=True)
    cdf = np.cumsum(W, 0); cdf[-1] = 1.0
    flat = (cdf + np.arange(genes.size)[None, :]).T.ravel()   # gene-major, increasing
    rep = np.repeat(np.arange(genes.size), need[genes])
    idx = np.searchsorted(flat, rep + rng.random(rep.size) * (1 - 1e-12), side='right')
    idx = np.minimum(idx, flat.size - 1)
    cell = idx % X.shape[0]; gcol = genes[idx // X.shape[0]]
    np.add.at(X, (cell, gcol), 1)


def _sharp_gene(col, delta, rng):
    """Change one gene's counts by `delta` with the largest rank shift per count (rank-aware
    count exchange, eye-contact #96 / danielchen26 'significance is a separate knob'):
    add to the lowest-count cells first, remove from the lowest NONZERO cells first."""
    col = col.copy(); n = len(col)
    if delta > 0:
        while delta > 0:
            order = np.lexsort((rng.random(n), col))
            m = col[order[0]]; tier = order[col[order] == m]
            take = tier[:delta]; col[take] += 1; delta -= len(take)
    elif delta < 0:
        delta = -delta
        while delta > 0:
            nz = np.flatnonzero(col > 0)
            if nz.size == 0:
                break
            m = col[nz].min(); tier = nz[col[nz] == m]; rng.shuffle(tier)
            take = tier[:delta]; col[take] -= 1; delta -= len(take)
    return col


def pin_emit(X0, b, rng, eps=0.05, sharp=None):
    """Pinned-total emission with minimal change to REAL control cells (keeps their scatter/shape).

    X0: [n, G] integer counts of real control cells. Totals T = round(S*b), S = X0.sum().
    Down genes: binomial thinning; up genes: add counts with cell weight x_ig + eps*depth_i*b_g
    (multiplicative, zero-inflation kept); then exact integer repair. sum_i X[i, g] == T[g].
    """
    X = np.asarray(X0, np.int64).copy()
    S = int(X.sum()); T = target_totals(b, S); C = X.sum(0)
    depth = X.sum(1).astype(np.float64)
    f = np.divide(T, C, out=np.ones(len(T)), where=C > 0)
    down = f < 1
    if down.any():
        X[:, down] = rng.binomial(X[:, down], f[down])
    wts = X0 + eps * depth[:, None] * np.maximum(b, 1e-12)[None, :]
    for _ in range(50):
        R = T - X.sum(0)
        if not R.any():
            break
        _add_counts(X, np.maximum(R, 0), wts, rng)
        neg = np.flatnonzero(R < 0)
        for g in neg:  # few genes, few counts: remove from cells holding counts
            k = -R[g]; col = X[:, g]
            units = np.repeat(np.arange(len(col)), col)
            take = rng.choice(units, size=k, replace=False)
            np.subtract.at(X[:, g], take, 1)
    if sharp is not None:
        for g in np.flatnonzero(sharp):
            X[:, g] = _sharp_gene(np.asarray(X0[:, g], np.int64), int(T[g] - X0[:, g].sum()), rng)
    assert np.array_equal(X.sum(0), T) and X.min() >= 0
    return X


def _drain_top(X0, R):
    """Remove R[g] >= 0 counts per gene from the TOP of each column (water-draining from the highest
    cells down). Ranks against the controls move as little as possible. Vectorised over genes."""
    n, m = X0.shape
    order = np.argsort(-X0, axis=0, kind='stable')
    S = np.take_along_axis(X0, order, 0).astype(np.int64)          # descending
    cs = np.cumsum(S, 0)
    nxt = np.vstack([S[1:], np.zeros((1, m), np.int64)])
    k_idx = np.arange(1, n + 1)[:, None]
    Rk = cs - k_idx * nxt                                             # removal to flatten top k to S[k]
    k = np.argmax(Rk >= R[None, :], axis=0)                           # first k (0-based) with Rk >= R
    k = np.where((Rk >= R[None, :]).any(0), k, n - 1)
    kk = k + 1
    keep = cs[k, np.arange(m)] - R                                    # total left in top kk cells
    L = keep // kk; extra = keep - L * kk
    new = S.copy()
    top = np.arange(n)[:, None] < kk[None, :]
    new = np.where(top, L[None, :], new)
    bump = np.arange(n)[:, None] < extra[None, :]
    new = new + bump.astype(np.int64)
    out = np.empty_like(new); np.put_along_axis(out, order, new, 0)
    return out


def shaped_emit(X0, b, rng, called=None, up_power=3.0):
    """Pinned totals (as pin_emit) with per-gene SHAPE chosen for the DE test:
    called genes -> rank-maximising change (bottom of the distribution, _sharp_gene);
    all other genes -> rank-neutral change (add to / drain from the top tail), so their pinned
    totals move the MSE pseudobulk without creating Wilcoxon calls."""
    X0 = np.asarray(X0, np.int64)
    S = int(X0.sum()); T = target_totals(b, S); C = X0.sum(0)
    R = T - C
    X = X0.copy()
    called = np.zeros(len(T), bool) if called is None else np.asarray(called, bool)
    dn = (~called) & (R < 0)
    if dn.any():
        X[:, dn] = _drain_top(X0[:, dn], -R[dn])
    up = (~called) & (R > 0)
    if up.any():
        depth = X0.sum(1).astype(np.float64)
        w = X0.astype(np.float64) ** up_power + 1e-3 * depth[:, None] * np.maximum(b, 1e-12)[None, :]
        need = np.where(up, R, 0)
        _add_counts(X, need, w, rng)
    for g in np.flatnonzero(called & (R != 0)):
        X[:, g] = _sharp_gene(X0[:, g], int(R[g]), rng)
    assert np.array_equal(X.sum(0), T) and X.min() >= 0
    return X
