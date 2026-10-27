**Rank A first, B second. Run A next: a two-arm, budget-matched H1 comparison.** A has measured improvements in response direction, PDS and MSE; B still requires discovering whether label-free bins separate beneficial amplitude changes from NMAE damage. A therefore has the better expected official gain per hour, although neither has a calibrated numerical forecast.

Two findings qualify A’s case:

- **The recorded CD4 check is not demonstrably the claimed one-family no-op.** cross_context.py (source reference: `mse_direction_1024/cross_context.py:26`; not included) constructs `Ensemble('CD4')` and explicitly excludes CD4-family banks. The one-family cancellation argument is correct for pure inter-family reweighting, but it does not describe this recorded comparison. Its negative result remains weak adverse evidence.
- The old `source_scorer_best` also had requested B/claim **0.8673**, below its recorded 0.90 gate. Its emitted capped-budget check passed. Those are distinct checks; the old gain cannot simply be called deployment-ready. Experiment report (source reference: `mse_direction_1024/REPORT.md:43`; not included)

**1. Why A wins the next slot**

A offers an already specified candidate with a positive measured direction signal. B offers a plausible mechanism but additional selection degrees of freedom: expression versus predicted-effect bins, cutoffs and multipliers. That costs more experiments and creates more SEARCH overfitting opportunities.

For perspective only, applying the supplied layout transfer factors to A’s three DE-component deltas changes its +0.0167 local overall to approximately **+0.0133**, holding expression-component deltas unchanged. That is a sensitivity calculation, **not an official forecast**: layout-derived slopes are not validated for source-direction changes.

The next A comparison costs approximately two H1 arm-hours, plus preparation. Its purpose is to establish whether the gain survives a consistent, label-free budget rule.

**2. Freeze this experiment**

Exactly two arms, with no selection based on the pending factorial:

| Arm | Frozen definition |
|---|---|
| A0 | Incumbent response direction, passed through the common budget calibration below. |
| A1 | The existing `family_quality` direction-replacement recipe: measured-gene family means; the existing capped square-root cell-count quality weights; unchanged missing-source fallback, own-target treatment and initial per-target norm matching. Then the same budget calibration as A0. |

Use the existing label-free B=claim solver and identical probe-selection rule in both arms. Freeze source files, recipient-family exclusion, genes, targets, control split, scorer TS=50,000, renderer, seeds, call-selection policy and all other settings. No new fusion exponent, amplitude grid or combination with an unmask winner.

These are **budget-matched prospective arms**. Do not substitute the archived incumbent score for A0, or assume A1 retains the old +0.0167. Record requested and emitted budget checks separately.

**Gate:** orient all component deltas so positive means better, using recovered official anchors rather than harness anchors. Define

\[
C=(\Delta P+\Delta M+\Delta N+\Delta F+\Delta R+\Delta J)/6
\]

and a prespecified sensitivity score

\[
T=(\Delta P+\Delta M+\Delta N+0.1\Delta F+1.5\Delta R+\Delta J/6)/6.
\]

Advance from H1 only if budget/emission checks pass and **both C and T are positive**. Report all six deltas and paired target-bootstrap intervals. Require no standalone +0.030 MSE hurdle.

For confirmation, require positive C and T on the second recipient and paired 95% lower bounds above zero there, recomputing aggregate metrics within each resample. H1 remains discovery evidence; these intervals quantify target-sampling uncertainty, not official transfer. Passing means “worth retaining for final selection,” not “proven official improvement.”

**Second recipient: HCT116 is structurally legitimate**, provided the actual eligible donor pool contains multiple families and excludes HCT116 throughout prediction/training. But it is not currently ready: the local HCT116 SEARCH truth matrices are absent, and `model_billm/checkpoints` contains only H1, CD4 and production checkpoints. The existing availability audit (source reference: `context_identity_1025/REPORT.md:94`; not included) documents the missing truth. Count preparation in the cost; never use a production model containing recipient-family information as its comparator. HEK293T has local truth but likewise lacks the matched incumbent checkpoint.

**3. The traps**

**A:** confusing correct scaling with correct transfer, and attributing a whole direction replacement’s gain specifically to quality weighting. Its selection among several H1 candidates also makes the reported bootstrap interval optimistic as evidence of generalization. The CD4 pool mismatch prevents discarding contrary evidence on structural grounds alone.

**B:** assuming high-expression/high-predicted-|lfc| genes and reference-significant genes are separable. They can overlap heavily. Your global-amplitude result already shows MSE/PDS gains overwhelmed by NMAE loss. Fixed predicted calls do not freeze the reference’s NMAE evaluation set; B=claim recalibration can also undo the intended bin allocation. Judge emitted profiles and all six components, not pre-calibration amplitude curves.

Read-only inspection only; no jobs, downloads or submissions. Safety entry unchanged.
