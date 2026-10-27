# DECISION 2026-10-25 03:32 ET — `family_quality` was rejected by the wrong gate, and it is the next arm

## The claim
The strongest uneliminated lead in the project is not a new idea: it is an **already measured H1 result that
was thrown away by a gate that priced the wrong component**. `mse_direction_1024` screened `family_quality`
(per-family fusion weighted by `sqrt(min(n_cells,1000))`), then rejected it because its **scaled MSE** rose
only +0.0162 against a pre-registered +0.030 hurdle.

## What the same run actually measured (H1 SEARCH, 79 targets, official scorer, paired)
| component (scaled, official anchors) | incumbent | source_scorer_best | Δ |
|---|---:|---:|---:|
| **PDS** | 0.5616 | **0.6326** | **+0.0710** |
| MSE | 0.0064 | 0.0226 | +0.0162 |
| NMAE | 0.0689 | 0.0745 | +0.0056 |
| fidelity | −0.2002 | −0.1811 | +0.0191 |
| reach | 0.0439 | 0.0344 | −0.0095 |
| Jaccard | 0.0334 | 0.0312 | −0.0022 |
| **sum** | 0.5140 | 0.5932 | **+0.0792** → **+0.0132 overall** |

Directional rho at the scorer coordinate rose **0.121141 → 0.143296** (+0.0222, paired 95% CI +0.0137..+0.0308).

The gate that killed it valued one component. The official composite is the **unweighted mean of six**, so the
correct decision quantity was the **sum**, which was positive on five of six components and positive overall.
Two independent attempts to kill this gain share one defect: they judged a multi-component trade on MSE alone.

## Why the "CD4 does not replicate" verdict is weak evidence, not a null
CD4's source panel contains **one family** (CD4 `Ctrl`/`Stim8h`/`Stim48h` are the same family by construction;
the repo also lists the family as `CD4` for all three). A family-weighted average across a single family is
**identically the unweighted sum** — the intervention cannot change anything on CD4 by construction. So CD4 is
the **wrong control** for this intervention, not disconfirming evidence. Astra independently flagged the same
point and added a second: the recorded CD4 check built an `Ensemble('CD4')` that excludes CD4-family banks
altogether, so it did not even exercise the arm it claimed to test.

I record CD4 as **NON-DISCRIMINATING for this lever**, not as a failed replication. That is a correction of a
recorded verdict, and it is why the lever stays open.

## What is still adverse and must be carried forward
`source_scorer_best` had requested **B/claim 0.8673**, below the 0.90 requested-profile lower bound (emitted
capped check passed). That is a real defect of that specific run and it is why this arm must be rebuilt with
**budget-matched** calibration, not reused from the archive.

## The next arm (frozen, two arms only, H1 first)
- **A0 incumbent**: current shipped direction, re-run through the frozen label-free B=claim solver.
- **A1 family_quality**: the same recipe with per-family `sqrt(min(n_cells,1000))` weighting, then the **same**
  solver, same probes, same renderer, same call policy, same seeds.
- Everything else frozen: source files, recipient-family exclusion, genes, targets, controls, TS=50,000,
  layout, scorer. No new exponent, no amplitude grid, no combination with an unmask winner.

### Gate (multi-component, anchored, pre-registered)
Using recovered **official** anchors, orient every delta so positive is better, and require **both**:
- `C = mean(ΔP, ΔM, ΔN, ΔF, ΔR, ΔJ) > 0`
- `T = mean(ΔP, ΔM, ΔN, 0.1·ΔF, 1.5·ΔR, ΔJ/6) > 0`   (Astra's transfer-sensitivity sensitivity score)

plus budget/emission checks (G2, G3e upper ≤1.05, G6, G7, G8) and **no standalone MSE hurdle**. Report all six
deltas with paired target-bootstrap intervals. Promotion = "worth retaining", never "official gain proven".

## Second recipient
HCT116 is structurally legitimate (multi-family donor pool, HCT116 excluded) but **not ready**: local HCT116
SEARCH truth matrices are absent and no matched incumbent checkpoint exists. Preparation is a separate costed
step, not a free replication. Do not substitute HEK293T without the matched checkpoint.

## Explicit non-actions
No package change, no rebuild, no upload, no submission; safety entry `5NeeoWn96lfm23hgExcm` unchanged.
`expr_unmask_1025` continues uninterrupted; this arm starts only after it finishes.

## Sources
`mse_direction_1024/REPORT.md` (table + gate), `mse_direction_1024/artifacts/comparison.json`,
`artifacts/cross_context.json`, `consult_1025/astra_2230.md`, `rsi/artifacts/source_inventory` family list.

---

## CORRECTION (03:40 ET) — I was wrong about CD4, and the correction changes the ranking

I wrote above that CD4's "fails to replicate" verdict was a *structural no-op* because its source pool
contains one family. **That is false, and Astra flagged it before I checked.** Direct inspection:

```
model_billm/data/*.npz families = [CD4, H1, HCT116, HEK293T, HepG2, Jurkat, K562, RPE1]
Ensemble('CD4') -> Sources(['CD4']) -> donor pool = 7 families (only CD4 excluded)
```

`Sources.__init__` skips banks whose **family** is in `excluded`; excluding `CD4` removes the three CD4
banks and leaves seven other families. So `family_quality` **is** a non-trivial intervention on CD4, and
`cross_context.json` **is** a real test of the lever on a second recipient. It came back
**−0.0019 rho (median), −0.0022 (count-weighted)**, both bootstrap intervals including zero.

**Consequences, stated plainly:**
1. The lever is **H1-positive, CD4-null-to-slightly-negative**. That is the pattern this project has been
   burned by twice: a single-recipient local gain that does not reproduce. It is *not* a clean win.
2. The MSE-gate mispricing is still a real defect in how that result was judged, and the +0.0132 H1 sum is
   still the honest H1 number. But **a mispriced test is not a passing test**: correcting the decision rule
   does not supply the missing second-recipient evidence.
3. Astra's rank-1 was conditioned on the CD4 evidence being weak. With that condition withdrawn, the two
   leads are much closer, and the priority should be decided by the arm that repairs the *most* uncertainty
   per hour.

**Revised next-arm order** (to be finalised when `expr_unmask_1025` completes, since its 2x2 answers the
upstream question of whether unmasking raises rho on *both* recipients):
- If the unmask factorial produces a both-recipient direction gain, that lever dominates and this
  source-fusion arm drops to a follow-up.
- Otherwise run the queued two-arm `family_quality` test **anyway**, because it is ~40 minutes and it
  repairs a genuine budget defect (B/claim 0.8673 → in-band) while producing the first *properly priced,
  budget-matched* six-component measurement of the lever. Treat its CD4 adverse evidence as binding: no
  promotion to "candidate for an official slot" without a second recipient showing positive `C` and `T`.

**Process lesson recorded:** I asserted a structural argument from a family list I had not inspected, and
it would have discarded real contrary evidence had Astra not challenged it. Claims about what an
intervention *cannot* affect must be checked against the code path that runs, not the data layout.
