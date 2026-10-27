You are deciding the next move for a team in the Virtual Cell Challenge 2026.
Be brutal. We want YOUR decision, not a menu.

## NEW HARD MEASUREMENT (official pipeline, same 79 targets, same renderer, same anchors)

Two arms have finished official scoring. This is not a proxy or a response-space ratio.

|                      | B0 shipped | B1 manifold F1 |
|----------------------|-----------:|---------------:|
| official score_avg   |  +0.021555 |      -0.055694 |
| raw expression MSE   |     1.4107 |     **0.9519** |
| raw PDS cosine       | **0.8500** |         0.5547 |

B1 BEATS the expression-MSE term outright (0.9519, better even than the anchor baseline 0.9613) and
still LOSES by 0.077, because PDS collapsed 0.8500 -> 0.5547.

Point accounting (each scaled component is 1/6 of the composite):
  pds_cosine          -0.10282
  de fidelity         +0.01424
  de lfc nmae         +0.00807
  de reach            +0.00195
  expr MSE            +0.00164
  de jaccard          -0.00032
  TOTAL               -0.07725
PDS loss is 62.6x the MSE gain. B1 improved THREE of four DE sub-metrics and still lost badly.

## The full component map of our actual entry (B0, official)
| component | scaled | points | headroom to anchor baseline |
|---|---:|---:|---:|
| pds_cosine | +0.7312 | +0.1219 | -0.1219 (we are ABOVE baseline) |
| de fidelity | -0.4795 | -0.0799 | +0.0799 |
| de jaccard | -0.1803 | -0.0301 | +0.0301 |
| de reach | +0.0456 | +0.0076 | -0.0076 |
| de lfc nmae | +0.0124 | +0.0021 | -0.0021 |
| expr MSE | 0.0000 | 0.0000 | +0.0016 to baseline |

Our 0.0216 = +0.1592 of good components (almost all PDS) minus -0.1376 of bad ones.
PDS alone is 5.7x our entire composite. If every component merely reached the anchor baseline the
score would be exactly 0.000.

## The tension you are being asked to resolve
Both Fable and Kimi instructed: "contest expression MSE; stop all DE work." Today's numbers say:
- expression MSE is worth at most +0.0016 at baseline, and a careless expression arm costs 0.1028 via PDS;
- the two largest recoverable blocks are DE fidelity (-0.0799) and DE jaccard (-0.0301);
- those are exactly the two metrics our DE work failed to move (effect-size gating: fidelity 0.4487 ->
  0.1511; recipient conditioning: a 4-arm loss).

So "contest expression MSE" appears to be chasing the smallest lever while two 40-80x larger ones sit
unmoved. Either the advisors were wrong about the objective, or there is a reason the DE blocks are
unreachable that dominates their arithmetic.

## Still running (finishing within the hour)
B2 (manifold F3, keeps the incumbent's projected pooled direction - the only arm that might hold PDS
while gaining MSE), then B3 (pooled), B4 (training-mean), B5 (shipped shrunk 0.5 toward training-mean).
Also already measured in response space: B4 training-mean beats B0's MSE by MORE than B2 does.

## Answer decisively - pick ONE
1. Given PDS is 5.7x our score and MSE is worth at most +0.0016 at baseline, is "contest expression MSE"
   still the right objective? If not, what replaces it, given that our last three DE experiments all
   failed? If yes, justify it against this arithmetic.
2. If B2 holds PDS (say >=0.80) while gaining MSE, is that the answer, and what is the gate for
   packaging it? If B2 also collapses PDS, what runs Monday?
3. There is a specific hypothesis we have NOT tested: PDS and MSE are anti-correlated, so the right move
   may be to keep B0's direction and only rescale its amplitude (fixed per-family scalar, fitted on the
   168 by leave-one-out, never on the 79). Response-space evidence says B0's optimal scalar is ~0.10.
   Would amplitude-only rescaling of the shipped model gain MSE points while preserving PDS? Is that
   worth running Monday ahead of reliability weighting and Sidechain?
4. Is the realistic objective still (a) rank, (b) 0.25, or (c) a component? Pick one, and say what to STOP.
5. GPU: still no?
