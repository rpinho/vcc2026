# Amendment RENDER_U — uniform-confidence layout ablation at the frozen 885-call budget (registered before scoring)

Authority: decisions/consult_1026_amp round 5 (2026-10-26 17:54 UTC).
- Q1 majority (Fable, Astra): Astra's uniform-confidence ablation; Opus's 2x budget is the minority and runs only
  if this ablation shows the layout matters.
- Q3 composition: one diagnostic arm, then stop unless it yields an actionable hypothesis (Fable: stop; Opus: one
  experiment with a fixed stop rule; Astra: continue bounded, stop earlier if nothing is actionable).
- DIAGNOSTIC ONLY, no promotion (Fable). Any later candidate needs Astra's confirmation on untouched targets plus
  Opus's Q2 list.

Arm U885 (results/predictions_U885.npz = B0 response verbatim; sha256 in results/render_u_arms.json):
- Confidence = B0's on ineligible genes; on the eligible genes (confidence >= 1e6, i.e. the same eligible set as B0)
  it is 1e6 + u, u ~ U[0,1) from np.random.default_rng(20261027). This is a preregistered, confidence-independent
  tie-break: B0's ranking is replaced by a seeded permutation, so which 885 genes are called and their goal order
  change, and nothing else does. Budget 885, seed base 20260910, anchors and depth are all identical to B0.
- Pre-check (Astra): 78/79 targets have > 885 eligible genes; the median overlap of the top-885 identities with B0
  is 0.217 (min 0.098). Post-check: counts_sha256 must differ from render_B0.json on those targets.
- Runner score_render_u.py = score_run.py confidence_all/run with textual substitutions only: confidence cache
  path -> scratch/confidences_RENDER_U.npz; render-json budget recorded from base['renderer']['budget'] instead
  of the literal 885 (Opus provenance fix).

Preregistered reading (Fable), paired composite delta d = U885 - B0 (20k boot, seed 20261026):
- |d| <= 0.002: the confidence allocation is inert. Close the renderer line; B0 is final for research.
- d < -0.002: confidence-guided allocation contributes. That explains B0 but is no entry; put Opus's 2x budget
  back to the advisors as the only remaining test.
- d > +0.002: an anomaly, not a candidate. Back to the advisors for a confirmation design on untouched targets.
PDS and raw MSE must be bit-identical to B0 (the layout preserves column sums); if they are not, the run is invalid.
