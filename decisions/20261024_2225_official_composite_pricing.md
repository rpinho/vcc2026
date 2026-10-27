# DECISION 2026-10-24 22:25 ET — the official composite is the unweighted mean of six components; use it to price levers

## Verification
`OFFICIAL_SCORES.md` gives per-component values for 7 entries. Reconstructing `overall` as
`(pds + mse + nmae + fid + reach + jac) / 6` reproduces every one within rounding:

| entry | reported | mean of six | diff |
|---|---:|---:|---:|
| v4a | 0.14380 | 0.14383 | +0.00003 |
| v5a | 0.14990 | 0.15000 | +0.00010 |
| v5b | 0.15280 | 0.15283 | +0.00003 |
| v8a | 0.14700 | 0.14683 | -0.00017 |
| v8b | 0.14120 | 0.14117 | -0.00003 |
| slim_v5b | 0.08833 | 0.08833 | +0.00000 |
| billmA | 0.16438 | 0.16450 | +0.00012 |

Max residual 1.7e-4, consistent with the 4-5 decimals published. Treat the rule as established:
**one scaled unit of any single component buys exactly 1/6 = 0.1667 of overall.** No component is
weighted more than another, so "which metric is worth the most" is purely "which has the most room
to move", and the denominator is the same for all six.

## What our current best (billmA 0.16438) is actually made of
| component | value | contributes | share of overall |
|---|---:|---:|---:|
| **pds** | 0.706 | **+0.1177** | **71.6%** |
| reach | 0.186 | +0.0310 | 18.8% |
| nmae | 0.116 | +0.0193 | 11.8% |
| mse | 0.061 | +0.0102 | 6.2% |
| jac | -0.016 | -0.0027 | -1.6% |
| fid | -0.066 | -0.0110 | -6.7% |

**PDS alone is 71.6% of the score.** Any lever that moves PDS moves overall 5x harder than the same
number of scaled units anywhere else — but PDS has been pinned at 0.666-0.708 across every entry we
have, and Sidechain (our reference reproduction at 0.1657 overall) sits at PDS 0.678. So PDS is not
currently movable, and the remaining fight is over the other five.

## Re-pricing the layout lever we already ran (v8a vs v5b, official, same totals)
| component | v5b | v8a | delta | overall effect |
|---|---:|---:|---:|---:|
| reach | 0.077 | 0.162 | **+0.085** | **+0.0142** |
| mse | 0.082 | 0.016 | **-0.066** | **-0.0110** |
| fid | -0.031 | -0.064 | -0.033 | -0.0055 |
| jac | -0.000 | -0.015 | -0.015 | -0.0025 |
| nmae | 0.109 | 0.102 | -0.007 | -0.0012 |
| **overall** | 0.1528 | 0.1470 | **-0.0058** | (sum = -0.0060) |

This is the single most important number in the project right now. The confidence-ordered call
layout **did exactly what it was designed to do** — reach doubled, +0.085 scaled, the largest single
component gain recorded — and it still **lost overall**, because MSE fell 0.066 and fidelity fell
0.033 on the same change. Reach is the *cheapest* component we own; MSE is the *most fragile*.

**Rule that follows:** any future layout/DE-call arm must be scored on MSE and fidelity *in the same
run as the reach gain*, and must beat `+0.0142` (the reach ceiling it can win) in costs. A layout
proposal that cannot show MSE preservation cannot win, no matter how large its reach gain looks.

## Where the real room is
- **fidelity** at -0.066 is the only *negative* component we carry. Getting it to 0 buys +0.011 overall
  — more than the entire reach lever ever bought — and getting it to +0.066 buys +0.022.
- **mse** at 0.061 vs Sidechain 0.151 is +0.090 of headroom = +0.015 overall, and it is where Sidechain
  beats us while having *identical* PDS (0.678).
- **reach** at 0.186 vs our own earlier 0.077 shows +0.109 of headroom on paper, but the only measured
  route to it costs more MSE than it gains (above).
- **jac** is negative and small; low priority but never let a layout arm make it worse without pricing it.

## What this forbids
Do not report a layout result as a win on its reach number alone. Our own history (v8a) is the
counterexample, and `OFFICIAL_SCORES.md`'s 10-19 lesson — "the MSE-credit model FAILED ... the
within-group layout drives MSE far more than min(B, claim) says" — is the mechanism.
