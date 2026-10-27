Round 9 of consult_1026_amp. You are a science advisor to a Virtual Cell Challenge 2026 team (deadline 11-05 23:59
UTC; frozen entry B0 = billmA 0.16438 validation, rank ~463; top-10 median 0.323; final phase blind; nothing is
submitted without Ricardo). Ricardo is not convinced progress is impossible and asked for public research. Read
competitors/IDEAS-2026-10-26.md, signaudit_1026/REPORT.md and calib3_1026/REPORT.md. <= 300 words.

## New from public research (cell-eval2 docs/brief, board history)
- Official fidelity F = k / max(n_pred, n_real); random signs score 0.5. Our official raw fidelity is ~0.49
  (v8a 0.4933, v9a 0.4923, SLIM 0.4981): below chance. All top-60 teams sit at 0.517-0.564. Board fit: scaled
  fid ~ 3.34 x (raw - 0.513), so reaching ~0.53 is worth ~+0.02 overall. The validation median n_real is 46; we
  call ~885 per target, so officially F ~ sign precision of our calls.
- Reach depends only on the order of calls (smallest p first) over reference-significant genes. Arc's docs value
  a habitual-direction head at about +0.006 overall. FlowState raised reach .246 -> .308 with everything else
  flat (+0.011).
- Skadi's "mixs75": call genes by expected reference detectability at a fixed budget: fid +.011, reach +.013,
  jac +.008 scaled, +0.0057 overall.
- Jaccard is about |intersection|/n_pred when n_pred >> n_real. X-team leads jaccard.
- MSE: the across-perturbation budget B binds for us (emitted B/claim ~0.5 under the v9 layout, known since
  10-22). Top-60 raw normalised MSE is 0.46-0.61; ours is 0.91.

## Our local sign audit (B0, analysis only)
Called-gene precision is 0.567, falling with confidence rank on weak+mid targets: 0.563 (top 50) -> 0.534
(200-300) -> 0.520 (750+). Weak+mid fidelity is ~0.55, which equals precision; strong targets are 0.24 because they
are under-called. The renderer flips 11.6% of signs, but the flipped calls are MORE often right, so the renderer is
not the cause. B0's intended sign is right on only 57% of truly DE genes. There is no down-call skew.
Implication: on the weak/mid-heavy official panel, our ~885 calls per target pay a rank-decaying precision. Fewer,
better-ordered calls should raise fidelity and jaccard officially, which the strong-heavy local harness cannot see
(it is 33% strong, against an estimated 3.5% officially).

## Rule
1. Pick the first 1-3 experiments (each ~11 min official local arm) from: (a) smaller flat budgets (e.g. 100/250)
   or a per-target budget from inference-time predicted n_real; (b) reach-head ordering (positions 1-20 go to
   sign-certain, likely-significant genes); (c) Skadi-style detectability call selection; (d) other. Name each
   arm exactly.
2. Decision instrument: the local harness is mis-weighted. Is it acceptable to read results on the 53 weak+mid
   targets, plus a composite reweighted to the official stratum mix, with calib3's component ratios? If not, what?
3. With 10 days left, can any such change honestly become a submission recommendation to Ricardo? Under what bar?
