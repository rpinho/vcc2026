"""Label-free direction/magnitude recipes evaluated by the selector. Every knob is a keyword so one component can be
switched on/off. Inputs: source banks (other families only), the destination context's control profile, the panel."""
import numpy as np

V4_WEIGHTS = {'K562_GWPS': 2.0, 'H1': 2.0, 'HCT116': 1.0, 'HEK293T': 1.0, 'CD4': 0.5,
              'K562': 1.0, 'HepG2': 1.0, 'Jurkat': 1.0, 'RPE1': 1.0}


def source_weight(bank, scheme):
    if scheme == 'equal':
        return 1.0
    if scheme == 'v4':
        return V4_WEIGHTS.get(bank.family if bank.family != 'K562' else bank.source, 1.0)
    raise ValueError(scheme)


def gather(banks, targets, exclude_family, min_cells=30):
    """Per-source centred responses for the panel: list of (bank, resp[T,G] with NaN where absent, has[T])."""
    out = []
    for b in banks.values():
        if b.family == exclude_family:
            continue
        has = b.has(targets) & np.array([b.n_cells.get(t, 0) >= min_cells for t in targets])
        if not has.any():
            continue
        r = b.rows(targets); r[~has] = np.nan
        out.append((b, r, has))
    return out


def centre(entries, targets, genes, mode='source'):
    """Centre each source: 'source' = source-wide mean over ALL its rows (full-screen centring);
    'panel' = mean over the panel rows it holds (what a 437-target cut does); 'none'."""
    gi = {g: i for i, g in enumerate(genes)}
    res = []
    for b, r, has in entries:
        r = r.copy()
        if mode == 'source':
            m = b.source_mean()['mean']
        elif mode == 'panel':
            rr = r[has].astype(np.float64); w = np.ones_like(rr)
            for k, t in enumerate(np.asarray(targets)[has]):
                j = gi.get(t)
                if j is not None:
                    w[k, j] = 0
            m = np.divide((rr * w).sum(0), w.sum(0), out=np.zeros(len(genes)), where=w.sum(0) > 0)
        elif mode == 'none':
            m = 0.0
        else:
            raise ValueError(mode)
        r[has] = r[has] - m
        res.append((b, r, has))
    return res


def pool(entries, targets, genes, weights='v4', family_balance=True):
    """Weighted mean of centred source responses. family_balance: a family with several screens (CD4 conditions,
    K562 GWPS + essential) shares one family weight (v4 behaviour for CD4). Returns pred, n_sources, n_families, agreement."""
    T, G = len(targets), len(genes)
    num = np.zeros((T, G)); den = np.zeros((T, G)); nsrc = np.zeros(T, int)
    fam_count = {}
    for b, r, has in entries:
        fam_count.setdefault(b.family, 0); fam_count[b.family] += 1
    fams = np.zeros((T, len(fam_count)), bool); fam_ix = {f: i for i, f in enumerate(fam_count)}
    stack = []
    for b, r, has in entries:
        w = source_weight(b, weights) / (fam_count[b.family] if family_balance else 1.0)
        m = (has[:, None] & b.measured[None, :])
        rr = np.where(m, r, 0.0)
        num += w * rr; den += w * m; nsrc += has; fams[has, fam_ix[b.family]] = True
        stack.append((rr, m))
    pred = np.divide(num, den, out=np.zeros_like(num), where=den > 0)
    return pred.astype(np.float32), nsrc, fams.sum(1), (den > 0)


def agreement(entries, pred, targets, genes):
    """Mean cosine between each source's centred row and the pooled prediction (own target excluded); NaN if < 2 sources."""
    T = len(targets); gi = {g: i for i, g in enumerate(genes)}
    out = np.full(T, np.nan); cnt = np.zeros(T)
    for b, r, has in entries:
        m = b.measured.copy()
        for i in np.flatnonzero(has):
            mm = m.copy(); j = gi.get(targets[i])
            if j is not None:
                mm[j] = False
            a = r[i][mm].astype(np.float64); p = pred[i][mm].astype(np.float64)
            na, npp = np.linalg.norm(a), np.linalg.norm(p)
            if na > 0 and npp > 0:
                c = float(a @ p / (na * npp)); out[i] = (0 if np.isnan(out[i]) else out[i]) + c; cnt[i] += 1
    out = np.where(cnt > 0, out / np.maximum(cnt, 1), np.nan)
    out[cnt < 2] = np.nan
    return out


def compat_mask(pred, dest_control, genes, tau_cpm=1.0):
    """Expression-compatibility masking: zero the response on genes not expressed in the destination controls
    (control CPM below tau)."""
    expressed = dest_control * 1e6 >= tau_cpm
    return pred * expressed[None, :]


def silence_gate(pred, targets, genes, dest_control, cpm_floor=1.0, cpm_full=10.0, floor=0.0):
    """Per-target silence gate: a CRISPRi target the destination barely expresses cannot respond.
    Row factor ramps from 0 at control CPM <= floor to 1 at >= full (log scale); targets absent from the axis keep 1."""
    gi = {g: i for i, g in enumerate(genes)}
    f = np.ones(len(targets))
    for i, t in enumerate(targets):
        j = gi.get(t)
        if j is not None:
            cpm = dest_control[j] * 1e6
            f[i] = float(np.clip((np.log(max(cpm, 1e-9)) - np.log(cpm_floor)) / (np.log(cpm_full) - np.log(cpm_floor)), 0, 1))
    f = floor + (1 - floor) * f
    return pred * f[:, None], f


def confidence_scale(pred, nsrc, agree, mode='none'):
    """Per-target amplitude allocation proportional to confidence (MSE-optimal a_p ∝ rho_p * d_p)."""
    if mode == 'none':
        return pred, np.ones(len(nsrc))
    if mode == 'nsrc':
        f = nsrc / (nsrc + 1.0)
    elif mode == 'agree':
        a = np.where(np.isnan(agree), 0.5, agree); f = np.clip(a, 0, 1)
    elif mode == 'agree_sqrt':
        a = np.where(np.isnan(agree), 0.5, agree); f = np.sqrt(np.clip(a, 0, 1))
    else:
        raise ValueError(mode)
    return pred * f[:, None], f


def sign_gate(entries, pred, targets, genes, min_net=1, z_tau=None, tau=0.05):
    """v5b-style cross-source call gate without variances: a source 'votes' on a gene when |centred response| > tau
    (log1p CP50K units); keep genes with net agreement >= min_net in the consensus direction, zero the rest."""
    T, G = pred.shape; sgn = np.sign(pred); net = np.zeros((T, G), np.float32)
    for b, r, has in entries:
        z = np.where(has[:, None] & b.measured[None, :], r, 0.0)
        sig = np.abs(z) > tau
        net += (sig & (np.sign(z) == sgn) & (sgn != 0)); net -= (sig & (np.sign(z) == -sgn) & (sgn != 0))
    mask = net >= min_net
    return pred * mask, mask


def deflate(pred, V, eta=1.0):
    """Remove the destination's own control principal axes (genes x k) from the response."""
    if V is None:
        return pred
    return pred - eta * (pred @ V) @ V.T


def blend_learned(pred, learned, coefficient=0.3):
    """Bill M mix: learned unit direction mixed into the pooled unit direction before amplitude (model_billm/network.blend)."""
    norm = np.linalg.norm(pred, axis=-1, keepdims=True); ln = np.linalg.norm(learned, axis=-1, keepdims=True)
    unit = pred / np.maximum(norm, 1e-12); lu = np.where(ln > 1e-12, learned / np.maximum(ln, 1e-12), unit)
    mixed = (1 - coefficient) * unit + coefficient * lu
    return mixed / np.maximum(np.linalg.norm(mixed, axis=-1, keepdims=True), 1e-12) * norm


# ----------------------------------------------------------------------------------------------------------
# Two-tier recipe: "not all genes are equal" (2026-10-21). Spend amplitude and calls on the targets and genes
# that actually respond; give everything else the baseline (a zero centred response).
# ----------------------------------------------------------------------------------------------------------
_SIGNAL = {}


def source_signal_table(path=None):
    """Cached noise-corrected per-(source, target) signal strengths (final_push/source_signal.py)."""
    import pandas as pd
    from pathlib import Path
    if not _SIGNAL:
        p = Path(path) if path else Path(__file__).resolve().parent / 'artifacts/source_signal.csv'
        d = pd.read_csv(p)
        for r in d.itertuples():
            _SIGNAL[(r.source, r.target)] = (float(r.signal), float(r.depletion))
    return _SIGNAL


def strength_score(entries, targets, sim=None):
    """Label-free predicted target strength for an unseen context: max over available source screens of that
    screen's noise-corrected signal for the target (optionally control-similarity weighted). Returns log scores;
    targets no source measures get -inf (treated as unknown, NOT as silent)."""
    tab = source_signal_table()
    out = np.full(len(targets), -np.inf)
    for b, r, has in entries:
        w = 1.0 if sim is None else float(sim.get(b.source, 1.0))
        for i in np.flatnonzero(has):
            v = tab.get((b.source, targets[i]))
            if v is None:
                continue
            s = v[0] * w
            if s > 0:
                out[i] = max(out[i], np.log(s))
    return out


def tiers(score, q_silent=0.33, q_strong=0.67):
    """Quantile tiers within the panel (so the rule needs no absolute calibration across contexts).
    Unknown-strength targets (-inf) are put in the WEAK tier, never silenced."""
    known = np.isfinite(score)
    t = np.full(len(score), 1, int)            # 1 = weak
    if known.sum() >= 10:
        lo, hi = np.quantile(score[known], [q_silent, q_strong])
        t[known & (score <= lo)] = 0           # 0 = silent
        t[known & (score >= hi)] = 2           # 2 = strong
    return t


def apply_tiers(pred, tier, weak=0.6, strong=1.0, silent=0.0):
    f = np.where(tier == 0, silent, np.where(tier == 2, strong, weak)).astype(np.float64)
    return pred * f[:, None], f


def responsive_gene_mask(entries, genes, quantile=0.5, min_sources=2):
    """Genes that respond across many targets and many screens: for each source, the per-gene mean squared centred
    response over its rows; a gene is responsive when it is above the within-source median in >= min_sources
    screens. Label-free (sources only) and panel-independent."""
    votes = np.zeros(len(genes), int); n = 0
    for b, r, has in entries:
        if not has.any():
            continue
        v = np.nanmean(np.where(has[:, None], r, np.nan) ** 2, 0)
        v = np.where(np.isfinite(v) & b.measured, v, 0.0)
        pos = v > 0
        if pos.sum() < 100:
            continue
        thr = np.quantile(v[pos], quantile); votes += (v >= thr) & pos; n += 1
    keep = votes >= min(min_sources, max(n, 1))
    return keep, dict(sources=n, genes_kept=int(keep.sum()), quantile=quantile, min_sources=min_sources)
