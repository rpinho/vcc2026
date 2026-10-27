Round 7 of consult_1026_amp. You are a science advisor to a Virtual Cell Challenge 2026 team (deadline 11-05 23:59
UTC; frozen entry B0 = 0.16438 validation, built with the 885-call v9 layout; final-phase scores hidden, validation
board frozen since 10-22; last entry counts; NOTHING is submitted without Ricardo). <= 250 words, one decision per
question. You may read h2h2_1026/results/d_gate.json, VCC_OVERNIGHT_LOG.md (search "21:45" and "fidelity hole"),
layout_cd4_1025/REPORT.md, FINAL_PUSH.md (search "monotone").

## Your R6 experiment: D1770 = B0 response + B0 confidence, layout call budget 1770 instead of 885. PASSED the gate.
score_avg 0.056045 vs B0 0.021555; paired composite d = +0.0345 [+0.0219, +0.0477]. Valid: pseudobulk and PDS
bit-identical, scaled MSE 0 for both, raw MSE ratio 0.994. Calls per target 378/1655/1770 (min/median/max).
Fidelity 0.4902 vs 0.4487 (scaled +0.150, most of the gain), jaccard 0.0592 vs 0.0430, reach 0.1343 vs 0.1266.
Per target, fidelity is better on 39 and worse on 38 (median 0); the top-5 targets give 38% of the summed gain.
Official fidelity = k / max(n_pred, n_real), with k = direction agreement.

## Conflicting evidence already on record (the official board, not local)
- Official validation: v5b (native dense layout, ~4,475 calls) = 0.1528 vs v9a (same biology/amplitude, 885-call
  v9 layout) = 0.16135. Denser calls LOST officially by 0.0085, while this same local harness preferred dense.
  The local harness REVERSED the sign of that official comparison.
- A 10-25 local model of 1200->3500 calls was monotone and pointed at that dense endpoint, so no sweep was run
  ("would be fitting H1").
- The CD4 second recipient replicated the local dense-layout gain, so it shares the harness bias and is not an
  independent confirmation of transfer.
- H1 local median n_real is 107 (34% < 10); 885 already over-calls most H1 targets.

## Rule
1. Does D1770's local +0.0345 carry evidential weight for the hidden test, given that the harness mis-signed the
   only official density comparison we have? Yes / no, and why.
2. Is any honest confirmation possible before 11-05 with no live board? Name it, or say none exists.
3. Recommendation to Ricardo (he decides submissions): keep B0 final, or present D1770 to him as a candidate? If
   the latter, with what caveat in one sentence.
