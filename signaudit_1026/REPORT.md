# signaudit_1026 — sign audit of B0's DE calls on the h2h2 H1 harness (analysis only)
Inputs: captured official DE tables h2h2_1026/artifacts/de_pred_B0.parquet and de_real.parquet, render_B0.json and
predictions_B0.npz. Scripts: audit.py, decomp.py, rank.py. Reproduced fidelity 0.44866 (official local 0.44872).

- Called genes (p_adj<0.05 in the prediction): 80,044, median 1001 per target (885 forced by the layout plus
  incidental ones); median n_real 322 locally.
- Sign precision on called genes is 0.567 overall: 0.677 on genes the reference also calls, 0.548 on genes it doesn't.
- By stratum (n_real <10 / 10-999 / >=1000): 9 / 44 / 26 targets; fidelity 0.547 / 0.550 / 0.243;
  precision 0.547 / 0.550 / 0.601. Strong targets lose fidelity through under-calling (n_real 2848 vs n_pred 1016).
- Forced layout calls have precision 0.554; incidental (natively significant) calls 0.645.
- Precision falls with confidence rank on weak+mid targets: top 50 0.563, 50-100 0.553, 100-200 0.544,
  200-300 0.534, 300-450 0.535, 450-600 0.525, 600-750 0.528, 750+ 0.520. Strong targets: 0.63 -> 0.59.
- The renderer flips 11.6% of called signs against B0's intended sign, but flipped calls agree with truth more
  (0.602) than kept ones (0.563). The renderer is not the cause. B0's intended sign is right on only 57.2% of
  reference-significant genes: the head itself is weak.
- No down-call skew: up calls agree 0.586, down calls 0.549. Among reference-null genes, 50.006% have negative
  real LFC. The sweep's null-skew idea (competitors/IDEAS-2026-10-26.md #2) is not supported locally.
- Official raw fidelity is ~0.49 for v8a/v9a/SLIM. The official panel is weak/mid-heavy (calib3: ~96%; val median
  n_real 46), so there fidelity is roughly precision, and official precision is roughly 0.49 (below chance).
