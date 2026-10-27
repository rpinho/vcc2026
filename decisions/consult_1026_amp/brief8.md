Round 8 of consult_1026_amp. You are a science advisor to a Virtual Cell Challenge 2026 team (deadline 11-05 23:59
UTC; frozen entry B0 = billmA, 0.16438 validation; final-phase scores hidden; nothing is submitted without Ricardo).
<= 250 words, one decision per question. Read calib3_1026/REPORT.md (the calibration) before answering.

## Correction to my R7 brief (your unanimous R7 ruling relied on it)
I said the local harness "reversed the sign" of the only official density comparison. That is true for the overall
composite only. On v5b minus v9a, all four DE component signs agree with the official board. Locally, fidelity is
inflated 9-12x and jaccard 7-9x, reach is understated (0.6-0.75x), and MSE is invisible. Native lost officially
through reach (-0.104 scaled). D1770 does not lose reach.

## Calibration result (Ricardo asked for it; analysis only, no new scoring)
- D1770 projects to about 0.000 to +0.004 official, central about +0.002, MSE included. That is at or below the
  +0.002 to +0.003 the same method assigns to pure seed noise (B0s1/B0s2).
- D1770's local reach gain equals a seed replicate's, so it is noise.
- Its fidelity and jaccard gains come from the 26 strong targets (n_real >= 1000), and the 41 mid targets lose
  fidelity.
- A stratum model fitted to the official v5b/v9a deltas implies the official panel is about 31% weak, 65% mid and
  3.5% strong (locally 15/52/33%). Break-even is about 4-5% strong.
- Verdict: indeterminate, about 60% positive, no size quotable.
- U885, A125 and A150 project negative.

Proposed test (not run): D1770 on the CD4 recipient (layout_cd4_1025 harness, frozen CD4 incumbent, budget 885 ->
1770, about 20 min). Its stratum mix (29/63/8%) is close to the fitted official mix. Rule: convert the raw
fidelity/jaccard deltas with the official slopes; <= 0 -> D1770 closed as officially null; clearly above the seed
floor -> "plausible" only.

## Rule
1. Does the correction change your R7 ruling (B0 final; D1770 not a candidate)? Say plainly.
2. Run the CD4 D1770 test? Yes / no. If yes, is the decision rule above right, or what instead?
3. If CD4 says "plausible", what (if anything) would make D1770 worth presenting to Ricardo before 11-05?
