"""Label-free, panel-agnostic dual-moment emitter with independently switchable knobs.

Heads (independent by construction, coupled only through the depth distribution):
  * BULK head  -> the group's count-weighted pseudobulk, which expr_mse and pds_cosine read.
  * CPM head   -> the per-cell arithmetic CPM mean, which lfc_nmae and the Wilcoxon table read.
AtlasShift's `dual_moment_counts` fits both moments on one integer matrix within the depth-feasible bounds
(pooled composition minus arithmetic composition = Cov(L, x)/mean(L), so the two heads can differ only as far as
the realised depth spread allows; the feasible projection is applied BEFORE the budget is solved, as in v7).

Knobs (each independently switchable; `Knobs.off()` reproduces the v7a behaviour):
  claim_rho        B / capped-claim target (v7a used .95). HARD GATE: the realised ratio must land in [0.9, 1.05].
  claim_jackknife  'control400' (jackknife of real 400-cell control draws at the context's own depth) and/or 'probe'
                   (our own emitted probe groups); the claim is min(probe, control400) per v7.
  amp_alloc        per-target amplitude allocation: 'uniform' (v5/v7) or 'magnitude_conf'
                   (a_p proportional to predicted magnitude x confidence, MSE-optimal allocation).
  silence          per-target silence gate from the context's OWN controls (target barely expressed -> amplitude 0).
  cpm_amplitude    per-cell CPM head amplitude on sign-confident genes (v7 used .9), 0 elsewhere.
  call_budget      per-target DE call budget ('n_real' = predicted responder count, or a fixed int, or None).
  pvalue_order     p-value goal order = confidence order (reach walks our order over the real DE pool).

Nothing here reads hidden labels: the claim, the silence gate and the compatibility mask come from the released
controls of the context being predicted; the direction and confidence come from public source screens.
"""
import dataclasses, json
import numpy as np


@dataclasses.dataclass
class Knobs:
    claim_rho: float = 0.97
    claim_source: str = 'min_probe_control400'   # 'probe' | 'control400' | 'min_probe_control400'
    amp_alloc: str = 'magnitude_conf'            # 'uniform' | 'magnitude_conf'
    amp_alloc_power: float = 1.0
    silence: tuple | None = (1.0, 10.0)          # (cpm_floor, cpm_full) or None
    compat_cpm: float | None = 1.0               # zero the response on genes under this control CPM, or None
    scatter: str = 'selective'                   # 'pool4' (v5/v7a) | 'pool1' | 'selective' (v7b: pool1 off-support, pool4 on-support)
    cpm_amplitude: float = 0.9
    call_budget: str | int | None = 'n_real'
    pvalue_order: bool = True
    # hard safety gate (the slim_v5b lesson: an amplitude knob must be gated on B/claim, never on call counts)
    ratio_gate: tuple = (0.9, 1.05)

    @classmethod
    def off(cls):
        return cls(claim_rho=0.95, amp_alloc='uniform', silence=None, compat_cpm=None, cpm_amplitude=0.9,
                   scatter='pool4', call_budget=None, pvalue_order=False)

    def asdict(self):
        return dataclasses.asdict(self)


def allocation(delta, confidence, knobs):
    """Per-target amplitude factors f_p with mean 1 (so the context's total budget is unchanged by the allocation;
    the overall level is set afterwards by the B = rho x claim solve).

    MSE-optimal allocation is a_p proportional to rho_p * d_p: we cannot know rho_p, so confidence stands in for it
    and the predicted magnitude d_p for the magnitude. Returns (factors, info)."""
    T = len(delta)
    if knobs.amp_alloc == 'uniform':
        return np.ones(T), dict(mode='uniform')
    norm = np.linalg.norm(delta, axis=1)
    c = np.asarray(confidence, np.float64)
    c = np.where(np.isfinite(c), c, np.nanmedian(c[np.isfinite(c)]) if np.isfinite(c).any() else 0.5)
    c = np.clip(c, 0, None)
    if c.max() > 0:
        c = c / c.max()
    f = (np.maximum(c, 1e-3) ** knobs.amp_alloc_power)
    # amplitude multiplies an already-magnitude-carrying delta, so only the CONFIDENCE part is new here;
    # the magnitude part is delta's own norm. Normalise to mean 1 to keep the budget solve well-posed.
    f = f / max(f.mean(), 1e-12)
    return f, dict(mode='magnitude_conf', power=knobs.amp_alloc_power, factor_min=float(f.min()),
                   factor_max=float(f.max()), pred_norm_mean=float(norm.mean()))


def silence_factors(targets, genes, control_probability, knobs):
    """Row factor from the context's own controls: a CRISPRi target the context barely expresses cannot respond."""
    if not knobs.silence:
        return np.ones(len(targets)), dict(mode='off')
    lo, hi = knobs.silence
    gi = {g: i for i, g in enumerate(genes)}
    f = np.ones(len(targets)); cpms = np.full(len(targets), np.nan)
    for i, t in enumerate(targets):
        j = gi.get(t)
        if j is None:
            continue
        cpm = float(control_probability[j] * 1e6); cpms[i] = cpm
        f[i] = float(np.clip((np.log(max(cpm, 1e-9)) - np.log(lo)) / (np.log(hi) - np.log(lo)), 0, 1))
    return f, dict(mode='control_cpm_ramp', cpm_floor=lo, cpm_full=hi, silenced=int((f == 0).sum()),
                   factor_mean=float(f.mean()), target_cpm_median=float(np.nanmedian(cpms)))


def compat_mask(delta, control_probability, knobs):
    """Expression compatibility: a gene the context does not express cannot carry a response."""
    if knobs.compat_cpm is None:
        return delta, dict(mode='off')
    keep = control_probability * 1e6 >= knobs.compat_cpm
    return delta * keep[None, :], dict(mode='control_cpm_floor', cpm=knobs.compat_cpm, genes_kept=int(keep.sum()))


def ratio_gate(B, claim, knobs, label=''):
    """HARD gate. slim_v5b shipped at B/claim 6.5 because its gate had been demoted to a diagnostic; this raises."""
    ratio = float(B) / max(float(claim), 1e-300)
    lo, hi = knobs.ratio_gate
    ok = lo <= ratio <= hi
    rec = dict(label=label, B=float(B), claim=float(claim), ratio=ratio, gate=[lo, hi], pass_=bool(ok))
    if not ok:
        raise RuntimeError(f'B/claim gate FAILED {label}: {ratio:.4f} outside [{lo}, {hi}] ' + json.dumps(rec))
    return rec
