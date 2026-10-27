# DECISION — 2026-10-25 23:15Z — fix the eligibility gate + per-target call budget

**Authorisation: Astra (`openai/gpt-6-astra`), consult `decisions/consult_1025_gate/astra.md`.**
Kimi and Fable had not returned at write time; their answers are appended verbatim when they land, and if
they contradict, the *narrower* proposal runs — I am not voting.

## What changed the diagnosis

`de_wilcoxon_direction_fidelity_yield_raw` was reverse-engineered exactly
(`diagnostics/FIDELITY_DOMINANCE.md`; verified to <1e-6 against the job's own `tidy_CD4_incumbent.csv` on
all 79 CD4 targets). It is:

```
fidelity_yield_raw(t) = k / max(n_pred, n_conf)  ==  precision x min(1, coverage)
```

and it is **99.4% of `score_avg`** (-31.53 of -5.2871 on CD4).

Per-target decomposition on CD4 (`diag_fid_1025/`, 77 targets):

- **32 of 77 targets call nothing** -> fidelity forced to exactly 0. They hold **9,899 of 19,658** `n_conf` genes.
- median `n_pred` = 1; mean `n_conf` = 255 (range 1 ... 3218).
- precision where we *do* call = **0.585** on the 16 healthy targets, ~= baseline anchor **0.5866**.

We are not out-signed. We are out-**called**, by a renderer gate:

```python
layout.py:145   eligible = keep & (confidence >= 1e6) & (direction == bulkdir)
confidence = net*1e6 + zs,  net = #(donor bank significant & agreeing) - #(disagreeing)
```

`zs` is O(1e4) at most, so `confidence >= 1e6` is exactly **`net >= 1`**: a gene must be called significant
*and directionally agreed* by a net majority of **cross-cell-line** donor banks — on a fold where per-gene
donor->target cosine is 0.013.

## Astra's verdict (verbatim)

> VERDICT: Run one **soft-evidence eligibility + cross-fitted target-budget renderer** experiment; leave the
> MLP, expression amplitudes, and donor weights unchanged.
>
> EXPERIMENT: Admit every structurally retained, non-on-target gene with a finite, nonzero frozen predicted
> effect. Replace donor-significance and donor-direction agreement requirements with ranking evidence—not
> vetoes. Select the top `B_t` genes, retaining their frozen predicted signs.
>
> PROTOCOL: Fit a shallow CPU count regressor for `B_t`, using public-training-family reference DE counts as
> labels and inference-available target features … Calibrate … on inner leave-family-out folds using the exact
> fidelity metric … freeze everything and render the untouched outer family; never access its treated cells or
> DE table.
>
> KILL TEST: Score the actual rendered cells, not a hypothetical selected-gene list. Reject if outer-fold mean
> fidelity fails to improve over 0.1790, or if rescued eligibility does not produce actual significant calls.
> Report previously-zero targets separately, plus precision, realized `n_pred`, and budget-versus-realized-call
> discrepancy.

## Why Astra's caveat is the correct one

The measured headroom (+0.432 leak-free at constant depth, +0.520 per-target oracle) is a **hypothetical gene
list**. Selecting `B_t` genes is not the same as *producing* `B_t` significant genes: the renderer's count
edits must actually push each selected gene's statistic past `p_adj < 0.05`. Adjudicability of the incumbent
ranking is 1.000, so the evaluation side is clean — but realisation is unverified and is exactly what must be
scored end-to-end. Arm A below measures realisation separately from budget choice.

## What runs — one job, two arms

- **Arm A — GATE.** Structural eligibility only (`keep`, non-on-target, finite nonzero frozen effect); rank by
  the existing `confidence`; budget = the current scalar (~940). Isolates *realisation*. Leak-free by
  construction: no reference statistic is consulted.
- **Arm B — GATE+BT.** Same widened gate; `B_t` predicted per target by a small CPU regressor fitted on the
  seven *training* families' reference DE counts, features from inference-time quantities only, selected by
  inner leave-family-out against the exact fidelity metric.

Both compared to the frozen incumbent on the official `cell_eval2` vcc2026 path, on **CD4 and H1** (the only
folds with outer-excluded checkpoints). All 8 folds carry the cheap map-level diagnostics.

## Pre-registered kill tests

1. Arm A must lift outer-fold mean `direction_fidelity` above **0.1790** *on rendered, officially scored cells*.
2. Realised `n_pred` must rise for the previously-silent targets. If eligibility widens but cells do not become
   significant, the gate is **not** the binding constraint and this whole diagnosis is wrong — say so plainly.
3. Arm B must beat Arm A on >=6 of 8 folds with a paired 95% CI excluding zero; otherwise it is capacity, not
   signal, and is rejected (the same standard that killed the B-shuffled map arm).
4. Neither arm may push `pds_cosine` below 0.48 or raise `expr_mse` more than 10%.

No submission. Frozen entry `5NeeoWn96lfm23hgExcm` untouched throughout.

---

## Kimi's verdict (returned 23:19Z, `decisions/consult_1025_gate/kimi.md`) — agrees on mechanism

> VERDICT: The donor-vote eligibility floor (net >= 1) is the disease; widening it to a calibrated
> model-confidence threshold removes the 42% zero-call collapse, and a per-target budget regression trained on
> the 8 public families captures the remaining oracle headroom.
>
> PROTOCOL: (1) leave-one-family-out on the public families, sweep a threshold tau on the model z-score
> (dropping the net*1e6 term) with a fixed 940 budget; select tau* maximising mean fidelity. (2) fit
> log(n_conf) ~ log(m), m = per-target count of genes with model |z| > tau*/2. (3) eligible = keep &
> (direction == bulkdir) & (z >= tau*), render exactly min(s_hat_n_conf, n_eligible) top genes by z.
> (4) Run the official scorer.
>
> KILL TEST: If tau* = 0 (direction-only gate) on the 8 public families with a fixed 940 budget fails to beat
> the replicate anchor (0.5996), model sign noise — not the gate — is the blocker and the experiment dies.

**Two advisors, same two chokepoints, independently.** Astra: soft-evidence eligibility + cross-fitted
budget, `confidence` as ranking. Kimi: replace the `net*1e6` term with a calibrated z threshold, regression
budget on model dispersion. These are the same change; Kimi's is the more specific parameterisation (it names
tau explicitly instead of "rank by confidence").

### One discrepancy, resolved conservatively as instructed
Kimi writes "trained on the 8 public families" and "the 300 panel targets". **For a CD4 or H1 outer fold that
is a leak**: the held-out family's reference DE counts cannot enter the budget regression. Astra's protocol
correctly says "the seven training families". Where the two disagree, the narrower one runs — per my own note
at the top of this file — so the budget model is fitted on the **7 training families for each outer fold**.
Kimi's tau sweep is retained because it is leak-free under leave-one-family-out (7 train, 1 validate, rotate).

### Consequence for the running job
`vcc_gate_1025` was launched with Astra's formulation, which is the narrower one, so **no change is needed**.
Kimi's tau-sweep is a refinement of the *ranking* inside Arm A and will be tried as a sub-arm only if Arm A's
realisation rate is high (i.e. only if widening eligibility actually produces significant calls). If
realisation is low, the tau sweep is moot and the diagnosis is wrong.

Fable had not returned at 23:20Z; the handler chain had it queued behind a Kimi retry loop. Still no vote from me.
