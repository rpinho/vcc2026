# DECISION 2026-10-25 03:40 ET — H1 factorial: the expression mask is PROTECTIVE; only the PC removal is (slightly) positive

## Complete 2x2 on H1_rich_search (paired, same run, official anchors, all required gates pass)
| arm | PC deflate | expr mask | raw PDS | raw MSE norm | alpha | B/claim |
|---|---|---:|---:|---:|---:|---:|
| `incumbent` | yes | yes | 0.753651 | 0.983537 | 0.4443 | 0.449 |
| `no_mask` | yes | no | 0.722650 | 0.991230 | 0.2611 | 0.428 |
| `no_pc` | no | yes | 0.756092 | 0.980166 | 0.4282 | 0.450 |
| `unmasked` | no | no | 0.734184 | 0.986038 | 0.2573 | 0.429 |

Official-anchored scaled deltas (positive = better):

| arm | ΔPDS | ΔMSE | ΔNMAE | Δfid | Δreach | ΔJaccard | **sum** | **Δoverall** | priced NAE rule |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `no_mask` | −0.0688 | −0.0081 | −0.0069 | +0.0047 | −0.0258 | −0.0012 | **−0.1060** | **−0.0177** | fail |
| `no_pc` | +0.0054 | +0.0035 | −0.0008 | +0.0067 | −0.0051 | −0.0001 | **+0.0096** | **+0.0016** | pass |
| `unmasked` | −0.0432 | −0.0026 | −0.0034 | +0.0019 | −0.0164 | −0.0016 | −0.0654 | −0.0109 | fail |

Values are from the final `score.json` files, repriced on the recovered **official** anchors by
`expr_unmask_1025/price_official.py` (the job's own scaled columns use local anchors, which clip MSE to 0
at our operating point and cannot resolve the endpoint).

**Verification status: PASS (arithmetic).** All four H1 arms have written `score.json`,
but the independent verifier crashed (`AttributeError: module 'verify' has no attribute 'verify_task'`)
after H1 finished; the driver was relaunched with `--resume-after-h1` and the verifier then ran to
completion. `artifacts/H1_rich_search/verification.json` reports **PASS independent reproduction** for all
four arms, with `archive_equal: true` for the incumbent (prediction SHA `2cbf35e2…` matches the archived
value), B recomputed independently, and fidelity/Jaccard recomputed from emitted significance states.
All required gates (`G2, G3_requested, G3e_actual_capped, G3e_input_control_cap, G6, G7, G8`) pass in all
four arms. The one earlier doubt — that a crashed verifier left `UNVERIFIED` scores — is resolved.

**Paired target bootstrap (2000 resamples, seed 20261025)** — this is the decisive column:

| arm | paired Δrho 95% | paired ΔMSE(norm) 95% | task_candidate |
|---|---|---|---|
| `no_mask` | [−0.0842, 0.0000] | [+0.0036, +0.0126] (worse) | **false** |
| `no_pc` | **[+0.0037, +0.0279]** | **[−0.0066, −0.0011]** (better) | **true** |
| `unmasked` | [−0.0369, +0.0016] | [−0.0003, +0.0059] | **false** |

Only `no_pc` has both intervals excluding zero in the beneficial direction. `no_mask` and `unmasked` are
rejected with resolved intervals, not merely point estimates.

**Caveat on the priced-NAE rule as implemented:** the bootstrap's NMAE convention is "increase in raw NMAE is
a loss", and under that convention `no_pc` shows a slightly *worse* raw NMAE (+0.00046). The official
anchor scaling reverses the sign because the recovered official NMAE anchor runs the other way from the local
one. Both readings are recorded; neither changes the verdict, because the NMAE move is ~0.0005 in either
direction and the PDS/MSE/rho intervals are resolved.

## Findings
1. **The mask is not suppressing information — it is removing noise.** Removing it (`no_mask`) costs
   **0.069 of PSD** — by far the largest single-metric loss we have measured on a direction arm — and also
   makes raw MSE *worse*. Cross-source sign agreement is doing real work.
2. **The 50-control-PC deflation is mildly harmful.** Removing it (`no_pc`) is the only positive arm:
   small gains in PDS, MSE, NMAE and fidelity, small reach loss, net **+0.0022 official**. It passes the
   pre-registered NAE break-even.
3. **There is a real interaction:** removing both (`unmasked`) is worse than removing only the PC
   (−0.0075). Without the mask, the PC removal no longer helps — the mask is what makes the deflated basis
   safe to use.
4. **The amplitude solver is the hidden variable.** `alpha` collapses from 0.4443 to ~0.26 whenever the mask
   is removed, and B/claim stays ~0.43 throughout. Unmasking therefore does not add usable amplitude: it
   changes the *shape*, and the B=constraint rescales it back down. Any future "more signal" experiment must
   report amplitude and B/claim alongside the shape, or it is measuring the solver.

## What this decides
- **Rejected for promotion:** the expression-mask removal, and the fully-unmasked arm. That line is closed on
  H1; CD4 is running for completeness under the frozen both-recipient rule.
- **Held for CD4 confirmation:** the PC-deflation removal (`no_pc`). +0.0016 official is **far too small to
  spend an official slot on** (it is inside the transfer haze), but it is the only direction arm in this
  family that is positive, and CD4 decides whether it is real or H1-specific.
- Expected decision on completion: if `no_pc` is positive on CD4 with an interval excluding zero, it becomes a
  *component* of a future package rather than a package on its own; if not, the whole expression-leg
  programme is closed and the queued `family_quality` arm is next.

## Contrast with the earlier MSE-gate failure
This factorial is the counter-example that justifies the multi-component gate adopted for
`family_quality`: `no_mask` looked like an MSE story (+0.0077 raw worse) but was really a **PDS** story
(−0.069). A single-component gate would have mis-described it in either direction.

No package change, no submission; safety entry `5NeeoWn96lfm23hgExcm` unchanged.
