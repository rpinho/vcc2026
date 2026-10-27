VCC consultation, 2026-10-25 00:05 ET. Read-only, answer in <=700 words, plain prose, no preamble.

## Just measured (accept as given; do not re-derive)
Expression-leg 2x2 factorial, paired, same run, official scorer, independently verified, 79 targets each
on two recipients. Official-anchored scaled deltas (positive = better); official overall = mean of six.

H1_rich_search (incumbent sum=0.514):
- no_mask (remove DE mask, keep PC deflation): PDS -0.0688, MSE -0.0081, NMAE -0.0069, fid +0.0047,
  reach -0.0258, jac -0.0012 -> sum -0.1060 (overall -0.0177). rho 95% CI [-0.084, 0.000].
- no_pc (keep mask, remove 50-control-PC deflation): PDS +0.0054, MSE +0.0035, NMAE -0.0008, fid +0.0067,
  reach -0.0051, jac -0.0001 -> sum +0.0096 (overall +0.0016). rho 95% CI [+0.004, +0.028], MSE CI excludes zero (better).
- unmasked (remove both): sum -0.0654 (overall -0.0109).
CD4_depth192_search (incumbent sum=0.199):
- no_mask: PDS +0.0144, MSE +0.0024, NMAE +0.0015, fid +0.0180, reach -0.0234, jac -0.0009 -> sum +0.0121 (overall +0.0020).
So the mask REMOVAL has OPPOSITE SIGNS on the two recipients (-0.106 H1 vs +0.012 CD4), while the PC
deflation removal is small-positive on H1 and still running on CD4.

Also measured: alpha collapses 0.444->0.261 (H1) and 1.546->0.795 (CD4) whenever the mask is removed, at
constant B/claim ~0.43/0.64. So removing the mask changes shape, and the B=claim solve rescales it down.
And PDS removes all 300 PANEL target genes from the comparison, so the mask's large PDS effect (0.069 H1)
operates entirely on NON-target genes.

## Standing facts
Official best billmA 0.16438: pds .706, mse .061, nmae .116, fid -.066, reach .186, jac -.016.
MSE is 43% of the gap to the top-10 median and is the component the local harness measures worst.
Final panel D/E/F is BLIND; only the last entry counts; 2 validation submissions/day; deadline Nov 5.
Safety entry 5NeeoWn96lfm23hgExcm is published with a hidden score. Never submit without Ricardo.
Local-to-official transfer: fidelity ~10x attenuated, Jaccard ~6x, reach ~1.5x amplified, MSE uncalibrated.

## Closed (do not propose)
Global amplitude/rho, pool_k, bulk multipliers, B/claim as a lever, learned-blend coefficient,
source breadth/family fusion (one H1-positive unpromoted attempt), per-target budget shrinkage,
density sweep, global dense layout for a package, re-solving alpha after the blend (measured to COST
-0.041 scaled MSE), public MOLM-13/U-87MG acquisition.

## Questions
1. The mask's sign reverses between recipients. Is there anything worth building on that, or is it a
   definitive close of the expression-leg programme? Answer yes/no first, then the reason.
2. Given the remaining ~7.5 hours of local CPU tonight and ~10 days to deadline, name the single highest
   expected-value next experiment and the exact arm set + promotion gate. It must be something a local
   harness can measure honestly, and it must not be on the closed list.
3. If your answer is "nothing local is worth it", say what SHOULD be done instead with the time (e.g.
   package verification, writeup, a specific public-data recon), and name the minimum bar for it.
4. What is the most likely way we fool ourselves tonight?
