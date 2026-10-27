**R6 rulings**

**Q1: Accept U885 as valid under the corrected check.** The clause existed to prove the change touched only the cell layout, and that holds. Rendered pseudobulk is bit-identical, PDS is bit-identical, and `expr_distance_unbiased` is identical (0.0020131 in both). Raw MSE moved only through the per-cell variance term. The scaled-MSE term is 0 in the composite for both arms, so d is all DE. This fixes a wrong premise; it doesn't rescue a bad result, because d doesn't depend on the corrected term. The preregistered branch applies: d = −0.0054 < −0.002, so confidence allocation contributes. Note in the record that the CI's upper end (−0.00056) is inside ±0.002, so the effect is real but its size is uncertain.

**Q2: Run the 2x arm (1,770 calls) once, now.** It is the preregistered next step and costs about 11 minutes. Opus's stop rule is binding and set in advance: point < +0.002 or lower bound ≤ 0 closes renderer research and makes B0 final. A pass is still diagnostic only. Promotion still needs Astra's confirmation on untouched targets before 11-05. No further arms after this.

**Q3: Replace "MSE bit-identical" with:**
- pseudobulk sha256 and PDS bit-identical to B0;
- scaled-MSE contribution to the composite exactly 0;
- raw `expr_mse_unbiased` reported and not worse than B0 by more than 0.1%.

That last check guards the clamped term, which could hide a worse raw MSE in hidden scoring.
