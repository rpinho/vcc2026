Round 6 of consult_1026_amp. You are a science advisor to a Virtual Cell Challenge 2026 team (deadline 11-05 23:59
UTC; frozen entry B0; final-phase scores hidden; we never submit without Ricardo). <= 200 words, one decision per
question. You may read h2h2_1026/AMENDMENT_RENDER_U.md and h2h2_1026/results/u_gate.json.

## Result of your R5 majority experiment (uniform-confidence ablation U885, B0 response, 885 calls, seeded shuffle)
score_avg 0.016127 vs B0 0.021555; paired composite d = -0.00543 [-0.00982, -0.00056] (20k target bootstrap).
Fidelity -0.0050, jaccard -0.0037, reach -0.0034, nmae ~0. PDS bit-identical. Rendered pseudobulk bit-identical.
Counts changed on 79/79 targets; median called-set overlap with B0 0.54.

## Correction to my R5 brief
I said official MSE is deterministic given the prediction. It is not: expr_mse_unbiased subtracts a per-cell variance
term, so the cell layout moves it slightly even with identical pseudobulk (U885 raw 1.40920 vs 1.41072, -0.11%; the
seed replicate B0s1 moved it by -2.6e-5). Scaled MSE is 0 for both (clamped), so the composite delta is all DE.
The amendment's validity clause ("PDS and raw MSE bit-identical") therefore fails literally, on my wrong premise.
The same clause sits in Opus's R5 promotion list.

## Rule
1. Accept U885 as valid under a corrected check (pseudobulk and PDS bit-identical; MSE moves only through the
   unbiased variance term), or treat it as void? If accepted, the preregistered branch is d < -0.002: confidence
   allocation contributes; Opus's 2x budget (1770 calls) is the remaining test.
2. Run 2x budget (one ~11 min arm, diagnostic, Opus's R5 stop rule: point < +0.002 or LB <= 0 -> close renderer
   research, B0 final), or stop now?
3. Replace "MSE bit-identical" in the promotion list with what?
