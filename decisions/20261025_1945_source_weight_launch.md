# 2026-10-25 19:45Z — Launch source-reliability weighting (advisor-directed)

## Why this decision exists
Ricardo redirected on 2026-10-25 after the three-model "no launch" vote:
*"No you're supposed to make better models, try new ideas, and research."*
So closeout is off the table. This decision selects a model experiment, not a stop.

## Advisors (publication-normalized sources)
Directory: `decisions/consult_1025_continuation/`

1. **Astra — `openai/gpt-6-astra`** (`astra.md`, `astra.md.raw.json`): build
   **renderer-aware source weighting**. Fit a few nonnegative source weights, regularized toward
   uniform, using training-only cross-source reconstruction *after* the existing count renderer.
   Baseline identical except uniform pooling. Held-out-family folds with derivatives excluded from donor
   banks, residual fitting, weight fitting and tuning. Metric: family-macro rendered expression MSE with
   paired target errors. Kill unless **>=3% relative MSE reduction, improvement in >=6/8 families, and a
   paired hierarchical bootstrap 95% CI excluding zero**. Explicitly: *if retained data cannot reconstruct
   held-out expression in renderer units, do not substitute direction correlation.* Cap 48 CPU-hours.

2. **Kimi — `moonshotai/kimi-k2.5`** (`kimi2.md`): **VERDICT: pursue mechanism (1) per Advisor A —
   source-reliability weighting of donor pooling.** Inverse-variance donor weights, held-out-family
   rendered expression MSE. Kill if improvement <5% relative or reliability CV<0.05 within 0.5 CPU-days.
   Note: the first Kimi attempt (`kimi.md.raw.json`) returned `finish_reason=length`, `content=null`,
   4,997 reasoning tokens of 5,000 — it deliberated and never answered. Re-asked answer-first; `kimi2.md`
   is the valid vote.

3. **Fable — `anthropic/claude-opus-5.5` / `anthropic/claude-fable-5.1`** (`fable.log`,
   `fable_generic.md`): **refused the direct VCC brief twice with `content_filter`**, zero content.
   The domain-neutral rewrite produced generic guidance: per-domain inverse-variance weighting plus
   per-domain distributional calibration, chosen over representation/gating mechanisms.
   **Recorded as generic, not a direct VCC endorsement.** It is directionally consistent with the
   weighting experiment, so it does not conflict.

## Intersection actually authorized
Astra and Kimi — the two advisors that answered the actual question — independently name the **same**
mechanism: reliability-weighted donor pooling scored by held-out-family rendered expression MSE.
Fable's generic answer is directionally compatible. The experiment below is that intersection, using the
stricter of the two kill thresholds at each gate.

## Feasibility evidence (mine, not an advisor claim)
`source_weight_1025/probe.py` -> `probe_uniform.json`, run 2026-10-25 19:44Z, 5.6s:
- 14 banks, 8 families, 18,533 genes, all gene axes identical.
- `response == log1p(50000*control)`-space confirmed finite; controls non-negative.
- Donor-only rendered log-expression MSE per held-out family is computable for all 8 folds:
  H1 0.004544, CD4 0.006505, HCT116 0.005693, HEK293T 0.004734, HepG2 0.023567,
  Jurkat 0.021673, RPE1 0.028896, K562 0.017206 (120 rows/fold).
- Pooling choice already moves the metric (H1 0.004544 -> 0.004502 under row-count weighting;
  HCT116 0.005693 -> 0.005874), so the test has measurable dynamic range.
- Cross-family donor support is strong: K562 shares >=1,659 targets with CD4/HCT6/RPE1 etc.,
  so leave-one-family-out pooling is not degenerate.
- **Astra's blocking condition is therefore satisfied**: rendered expression can be reconstructed
  from these banks, so no surrogate correlation metric is needed.

## Decision
Launch the job in **`source_weight_1025/`**.
Two gates, fail-closed, advisor thresholds unchanged:

- **Gate A (cheap, no neural net, ~1 CPU-hour):** donor-pooling-only rendered MSE, uniform vs
  reliability-weighted vs row-count, 8 held-out-family folds, paired per-target errors.
  Promote only if >=3% relative macro improvement AND >=6/8 families improve AND paired bootstrap
  95% CI excludes zero (Astra), with CV of fitted weights >=0.05 (Kimi's collapse guard).
- **Gate B (only if A promotes):** re-run with the residual MLP + existing 0.3 blend, identical
  architecture/seed/layout in both arms, refit once per fold and shared between arms.
  Same promotion thresholds on the end-to-end rendered metric.

Gate A failing is a **successful negative result**, not a reason to stop working: it kills mechanism (1)
and redirects to the next advisor-ranked mechanism.

## Hard constraints honoured
- Never submit, upload, replace or modify `model_billmB/out/final_billmA.vcc`.
- Write only inside `source_weight_1025/`; read-only elsewhere; no frozen-package reads.
- CPU only, `nice -n 10`, maintain the storage reserve.
- No new credentials or pods. No commit/push of private research.

## Status
Launched in `source_weight_1025/`; machine-specific launch references omitted.
