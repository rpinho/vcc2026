# 2026-10-25 21:00Z — Unanimous: run amplitude calibration, scored on the official composite

## What preceded this
Mechanism 1 (reliability-weighted donor pooling) was preregistered, run, independently verified and
**rejected** at 20:42Z: macro MSE +0.206% (worse), paired 95% CI [-1.470%, +1.343%], 5/8 families.
See `source_weight_1025/REPORT.md`. Astra's and Kimi's earlier mechanism-1 recommendation is dead.

## The new evidence put to the advisors
`source_weight_1025/diagnose.py`, `diagnose_scale.py` (honest; scalars fit on the seven training families
only, applied to the held-out family):
- Uniform-pooling baseline family-macro MSE **0.015061**.
- One global nonnegative scalar: **0.013189 = -12.43% relative, 7/8 families improve**.
- Fitted scale stable across all eight folds: 0.2547, 0.2715, 0.3140, 0.3150, 0.3211, 0.3240, 0.3448, 0.3568.
- Mean `|pred|/|true|` = 0.70; per-gene cosine only 0.013–0.218.
- Ceilings: per-row oracle scale only reaches -17.57%; predicting ZERO reaches -8.20%.

## Votes (verbatim files in `decisions/consult_1025_amplitude/`)
- **Astra** (`astra.md`): YES. Apply one learned nonnegative scalar **after** the incumbent's final
  norm-rescaling; freeze MLP, 0.3 blend, pooling. Fit on leakage-free predictions from the seven training
  families. Evaluate on the **official composite and every DE component** with paired family-level
  uncertainty. Promote only if the composite gain has a positive lower 95% bound **and** preregistered DE
  noninferiority margins hold. Key correction: *fit the scalar to the actual incumbent output, not the
  pooled-response scalar* — the learned residual changes the optimum.
- **Kimi** (`kimi.md`): YES. Learn a scalar initialized ~0.3 instead of forcing the rescale-to-pooled-norm.
  Kill if official MSE improvement <5% or any DE metric degrades >10%.
- **Fable** (`fable.md`): YES, but as a *cheap gated check*, not the headline bet. Sharp caveat:
  **zero-response already gets -8.20%, so most of the -12.43% is shrinkage compensating for poor
  direction.** "Direction is the real gap to 0.472." Sweep s in {0.2,...,1.0}, choose by official
  composite. Kill if composite gain < +0.005, CI includes 0, selected s >= 0.85, or reach/Jaccard losses
  offset the MSE gain. Next big experiment either way: improve per-gene response **direction**.

**Unanimous: yes.** All three independently require scoring on the official composite with DE components
reported — not on log-MSE alone. That is the binding constraint this experiment must honour.

## Decision
Run `calib_1025` as a bounded, preregistered, gated check:

1. Freeze the incumbent (`model_billmB` MLP + 0.3 blend + donor pooling + renderer). Do not retrain it.
2. For each of 8 held-out-family folds, fit ONE nonnegative scalar `s` on the **actual incumbent's**
   post-rescale predictions over the seven training families, then apply that `s` to the held-out family.
3. Score with the real `cell_eval2` vcc2026 profile (via `.venv`), reporting **all** components:
   `expr_mse_unbiased_capped_norm`, `pds_cosine`, `de_wilcoxon_lfc_nmae`,
   `de_wilcoxon_direction_fidelity_yield_raw`, `de_wilcoxon_direction_reach_raw`,
   `de_wilcoxon_sig_jaccard`, plus the composite.
4. Paired family-level bootstrap CI on the composite change.

Promotion requires the **intersection** of the advisors' gates:
- composite gain >= +0.005 (Fable) with a positive lower 95% bound (Astra, Fable), and
- selected `s` < 0.85 (Fable), and
- no DE component degrading beyond the stated margin, i.e. not >10% relative (Kimi) and not offsetting
  the MSE gain (Fable).
Kimi's 5%-of-MSE figure is reported separately.

## Honest framing to carry into the report
Fable's correction matters and must not be lost: a large part of the -12.43% is shrinkage, not new signal.
If the official composite does not improve, the correct conclusion is **"amplitude is a symptom, direction
is the disease"**, and the next experiment is per-gene direction quality — which is also what the low
cosines (0.01–0.22) independently suggest.

## Constraints
Frozen entry `5NeeoWn96lfm23hgExcm` and `preserve_1025/` untouched. No submission. Write only in
`calib_1025/`. CPU only, `nice -n 10`.
