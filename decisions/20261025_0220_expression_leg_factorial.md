# DECISION 2026-10-25 02:20 ET — test the expression-leg mask/deflation as a 2x2 factorial

## What was decided
Launch one paired emitted job, `expr_unmask_1025`, on `H1_rich_search` and `CD4_depth192_search`, with four arms:

| arm | 50-PC `L-Lc` projection | mask applied to profile |
|---|---|---|
| incumbent | on | on |
| no_pc | off | on |
| no_mask | on | off |
| unmasked | off | off |

Everything else is frozen: response head, learned blend coefficient 0.3, per-arm label-free
amplitude solve to `B = claim`, integer totals, 885-call v9 layout, confidence ordering, scorer,
seeds. The mask is still computed from the **incumbent** deflated direction and used as the call-set
generator in every arm. An incumbent replay must reproduce the archived prediction SHA in the same
run before any variant is believed.

## Why (measured basis)
1. MSE is the largest open component. 289 leaderboard entries with raw PDS within 0.03 of ours have
   median implied `rho=sqrt(1-u)` **0.363**, max 0.621, versus our **0.154**
   (`decisions/20261025_0115_mse_headroom_is_unlocked.md`). Same discrimination, 2.4x the usable
   direction: the difference is how the response is constructed.
2. The shipped profile applies two strong, precisely located suppressions
   (`rsi/v5b_loco_baseline.py:107-118`): `L = Lc + deflate(L-Lc,V)` with 50 control PCs, then
   `L = Lc + np.where(mask,L-Lc,0)` with the 18,533-gene cross-source mask. Either can remove real
   target-specific response, and the mask zeroes most of the DE universe.
3. Prior work closed **global amplitude/rho, pool_k, bulk multipliers, B/claim, PC projection alone,
   learned-blend coefficient, source breadth, per-target shrinkage, and global dense layout**. The
   PC screen varied the projection while leaving the mask in place; the layout work varied call
   placement, not the expression leg. The joint and mask-only arms are genuinely untested.
4. `emitter_screen_1025` screened blend 0.0 vs 0.3 at the correct coordinate but its arms rebuild
   direction from the ensemble; it is **not like-for-like** and is not being promoted.

## Pre-registered reading (frozen before results)
A variant is a candidate for the next stage only if, **on both recipients in the same run**:
- raw `expr_mse_unbiased_capped_norm` improves versus the incumbent replay;
- the paired target-bootstrap 95% interval for the MSE/rho delta excludes zero;
- raw PDS falls by at most 0.010; and
- each of fidelity, Jaccard, reach, and NMAE falls by at most 0.020 raw.

Otherwise it is not promoted and the lever closes. This is a mechanism test, not an official-score
forecast; local fidelity attenuates ~10x, Jaccard ~6x, reach amplifies ~1.5x, and MSE is the least
transferable component. **No package change, no submission, safety entry unchanged.**

## Falsification / what a null means
If all three variants are null or worse, the expression mask and control-PC projection are not the
missing direction signal, and the next work must be new evidence or a different model class rather
than another transform of the same banks. If `unmasked` helps but neither single-factor arm does,
report the interaction and do not attribute it to either mechanism alone.

## Status
Job directory: `expr_unmask_1025/`. Nothing submitted; safety entry unchanged.
