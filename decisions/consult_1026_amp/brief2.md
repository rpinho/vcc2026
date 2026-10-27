Round 2 of consult_1026_amp. You are a science advisor to a Virtual Cell Challenge 2026 team; the
orchestrator executes your rulings. One decision per question, mechanically executable. Short answers.

## Round 1 rulings (one paragraph each, for context)
- FABLE: two arms on the 79, mean-extrapolation at a single alpha* chosen on the 168 by an MSE-PRICED objective:
  A = mu + a*(B3 - mu) (primary), B = mu + a*(B0 - mu). Priority A > B3 > B > freeze B0. Rule stated in advance:
  "If the argmax is 1.00, the alpha line is dead for the leaderboard and only B3 proceeds." B3 is a candidate
  (fallback). Confirmation: CD4 fold check, then the hidden leaderboard (only Ricardo can submit; the
  orchestrator never submits).
- OPUS: one packaging candidate N1 = mu + 1.25(B0 - mu), rescaled per target to ||B0_i|| (direction change only,
  amplitude held), plus raw A1 (descriptive). Not B3 (its +0.0041 is the campaign's noise floor). Confirmation:
  cross-fitted official run on the 168 (2-fold, alpha chosen on the other half, pooled() base).
- ASTRA: E0 = 1.25 B0 - 0.25 mu and E3 = 1.25 B3 - 0.25 mu, alpha frozen at 1.25. B3 is a candidate. No clean H1
  confirmation cohort exists among the 247 targets; no independent cohort => no promotion, B0 ships.
- ALL: no pure scaling, no dispersion ladder, no shrink; 79-gate = score_avg > 0.021555, PDS >= 0.845, paired
  composite lower bound > 0; P1a gene-wise magnitude calibration follows and is built in parallel now
  (fit on 168 + public non-H1 only).
- Checked: pooled() excludes every H1 bank, so no target sees its own truth (Fable's leak check: PASS).

## THE MEASUREMENT THAT CONTRADICTS THE MAJORITY (Opus + Astra freeze 1.25)
Fable's MSE-priced objective on the 168 (training-only, pooled base, LOO mean, Fable's own pre-declared grid
1.00..2.00 step .05, Fable's formula: dPDS_comp = 1.17*(rank(1)-rank(a))/168/3;
dMSE_comp = (max(0, 1 - 0.939*r(a)) - 0.061)/6, r = response-space MSE ratio vs alpha=1):

| alpha | rank/168 | MSE ratio | dPDS_comp | dMSE_comp | Obj |
|---|---:|---:|---:|---:|---:|
| 1.00 | 25.095 | 1.000 | 0 | 0 | **0** |
| 1.05 | 23.470 | 1.051 | +0.00377 | -0.00804 | -0.00427 |
| 1.10 | 22.345 | 1.106 | +0.00638 | -0.01017 | -0.00378 |
| 1.25 | 21.351 | 1.285 | +0.00869 | -0.01017 | -0.00148 |
| 1.30 | 21.363 | 1.351 | +0.00866 | -0.01017 | -0.00150 |
| 1.50 | 22.167 | 1.641 | +0.00680 | -0.01017 | -0.00337 |
| 2.00 | 24.369 | 2.566 | +0.00169 | -0.01017 | -0.00848 |

ARGMAX = 1.00. Even +5% response MSE zeroes our leaderboard scaled MSE (0.061 -> 0), costing ~0.0102, more
than the best PDS gain (+0.0087 at 1.25). Under Fable's own rule the alpha line is dead for the leaderboard.

Caveats you must weigh: (i) 0.061 is our VALIDATION-panel scaled MSE (billmA); the final panel
(vcc2026-test-1, contexts D/E/F) score is hidden, so the MSE term there is unknown, but billmA's 0.061 is the
best estimate we have. (ii) Proportionality of official MSE to response-space MSE is an assumption; locally
B1's response-MSE ratio 0.71 vs official raw-MSE ratio 0.675 - roughly proportional. (iii) Opus's N1
(norm-matched) holds per-target amplitude; its MSE ratio is not yet computed and could be much smaller than raw
extrapolation's, because part of the MSE rise is amplitude, not direction. Rank is unchanged by norm matching.
(iv) Locally, on the 79, scaled MSE is clamped at 0 for B0 and B3, so a local official run cannot see this cost.


## COMPUTED SINCE (resolves caveat iii): N1 norm-matched, same 168 training-only setup
| alpha | raw-extrap MSE ratio | N1 MSE ratio | rank/168 (same for both) | Fable Obj raw | Fable Obj N1 |
|---|---:|---:|---:|---:|---:|
| 1.10 | 1.106 | 1.0050 | 22.345 | -0.00378 | +0.00560 |
| 1.25 | 1.285 | 1.0111 | 21.351 | -0.00148 | **+0.00697** |
| 1.50 | 1.641 | 1.0186 | 22.167 | -0.00337 | +0.00389 |
So the MSE cost is almost entirely the amplitude part. Norm-matching (Opus) keeps the PDS gain and loses only
~0.0017 of leaderboard MSE: under Fable's own objective N1 at 1.25 is positive and is the argmax of the three.

## Operational facts
- Every official run (any arm, and even the B3-B0 bootstrap) is currently BLOCKED: free disk 29.7 GB, below the
  30 GB guard, because of another session's finished 8.7 GB Rust build in /tmp; Ricardo is being asked to clear it.
- An official 168-target run (Opus's confirmation) needs cell-level task data for the 168 that is not yet built.

## Rule on these, briefly
1. Given ARGMAX = 1.00, do any mean-extrapolation arms still run on the 79? If yes, which (E0 / E3 / N1) and why
   the MSE pricing does not kill them. If the answer depends on N1's MSE ratio, say so and give the threshold.
2. One confirmation protocol, chosen from: (a) none possible -> B0 ships unless a candidate clears the 79 gate by a
   margin you name; (b) CD4 fold check; (c) cross-fitted 168 official run; (d) Ricardo's leaderboard submission
   (his decision, recommend only). Pick one and say what result promotes.
3. B3: candidate or not? (Fable, Astra yes; Opus no.) If candidate, what does it need beyond the 79 result?
4. P1a: confirm it is now the primary lever, given that MSE pricing is what killed alpha. Any constraint to add.
