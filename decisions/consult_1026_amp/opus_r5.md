**Q1: Run Fable's 2× budget (1,770 calls) on the B0 response. Skip the uniform-confidence ablation for now.**
`layout.emit` asserts that per-gene count sums and per-cell depth are preserved (`score_run.py:71-72`). So any layout change leaves PDS and raw MSE exactly where they are and moves only the cell-level DE metrics (reach, fidelity, jaccard). That makes budget the cleanest lever left: Opus's MSE gate is met automatically, and only the composite has to clear it. Uniform confidence keeps the same 885 calls, so it only explains why B0 works. Its expected effect is smaller and it is unlikely to produce an entry.
Fix this before launch: line 91 hardcodes `budget=885` in `render_{arm}.json`. Record `base['renderer']['budget']` instead, or the provenance will be wrong.

**Q2: Yes, it's legitimate.** The submission is just counts, and the renderer is our own pipeline. The risk isn't the rules; it's overfitting to 79 targets. A renderer change replaces B0 only if all of these hold:
- Paired composite lower bound > 0 on H1 SEARCH.
- The gain repeats under two more render seeds, each beating B0 by more than 0.0008 (the seed noise floor).
- PDS and MSE are bit-identical to B0 (check the asserts and the hashes).
- The gain is spread across targets, not driven by a few large reach swings.
- The submission file passes full format validation.
- Ricardo signs off.

**Q3: Continue, but only on this one experiment, with a stop rule fixed in advance.** Budget is the one lever with a mechanism that can't damage PDS or MSE, and an arm costs about 11 minutes. That doesn't justify opening a tenth lever.
Stop rule: if the 2× composite point estimate is below +0.002 or its lower bound is ≤ 0, close renderer research, make B0 final, and stop. If it passes, run 1.5× and 3× to check the trend is monotone, then go to the Q2 evidence. B0 stays submitted throughout.
