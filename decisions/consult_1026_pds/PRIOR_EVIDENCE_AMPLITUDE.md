# Prior evidence the advisors' Monday recommendation has ALREADY been partly run

Written 08:20 EDT before launching anything. All three advisors recommended "amplitude-only
calibration of the shipped model" as Monday's first experiment. `mse_1025` already ran a version of it.

## What was already measured (mse_1025, CD4 SEARCH fold, 79 targets, official pipeline, integer v9)
| arm | raw expr MSE | raw direction fidelity |
|---|---:|---:|
| incumbent (integer v9) | 1.28615 | 0.176704 |
| **source-fitted amplitude scalar 0.343363** | **1.11041** | **0.125058** |
| control mean | 1.10871 | (n/a) |

So a source-fitted scalar of **0.343** was tried, frozen before scoring, and produced:
- **MSE −13.66%** (1.28615 → 1.11041), nearly reaching the control-mean arm (1.10871, still 0.15% worse)
- **fidelity −0.0516 raw** (0.17670 → 0.12506), a **29% relative drop**

`mse_1025/REPORT.md` verdict: *"Fidelity falls from 0.176703749 to 0.125058135. Do not promote it as an
expression solution."*

## Why this matters for Monday
1. **The scalar direction is not untested.** What is untested is the *composite*: `mse_1025` recorded
   only raw MSE and fidelity for the scalar arm, not PDS, jaccard, nmae, or reach. The advisors' claim
   that "PDS preservation must be measured" is exactly right, and it is the one thing we lack.
2. **The fidelity cost is real and large**, not hypothetical. Fable predicted a scalar >1 might *raise*
   fidelity by fixing under-calling. On CD4 the opposite happened at alpha=0.343: fidelity fell 29%.
   That is evidence *against* the "B0 is over-confident, shrink it" story — at least on CD4.
3. **The break-even arithmetic is therefore tighter than Round 2 assumed.** If alpha ≈ 0.343 costs
   fidelity, and the H1 evidence says fidelity is the single largest recoverable block
   (−0.0799 pts), then a scalar that costs fidelity is spending from the very block we want to grow.

## What remains genuinely open, and is worth Monday
- **PDS was never measured for a scalar arm.** If PDS is truly scale-invariant through the renderer
  (Astra's caveat: it may not be), the trade is MSE-gain vs fidelity-loss only.
- **The fitted value may be wrong.** 0.343 was "source-fitted" on some earlier job's protocol. The
  advisors want alpha fitted by LOO on the 168 — a different, cleaner estimate. It could be nearer 1.0.
- **H1 vs CD4 may differ.** The 0.10 response-space optimum I measured was on H1.

## Concrete recommendation to put back to the advisors
Do not run a blind alpha sweep. Run **one arm per candidate alpha** through the official pipeline on H1
with all six components: alpha fitted-by-LOO-on-168 (the advisors' proposal) **plus** alpha=1 (B0
itself, already have it) **plus** the CD4-derived 0.343 as a replication check. That is 2 new scoring
runs, ~30 min, and it answers "does the scalar hold PDS and does it gain composite" directly rather
than by inference from a partial CD4 table.

## The uncomfortable possibility this raises
If amplitude calibration also loses (scalar costs fidelity, PDS is unchanged, net negative), then every
single lever we have tried — DE gating, recipient conditioning, manifold direction, amplitude
calibration — has failed, and the honest conclusion is that B0 is at a local optimum for this model
family and the remaining upside requires a genuinely different model, which 10 days and a 16 GB CPU box
cannot deliver. That is worth knowing explicitly rather than discovering after another night.
