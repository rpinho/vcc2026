# DECISION 2026-10-25 04:40 ET — no new score-seeking lever tonight; finish the registered family test

## Consulted (independently, no shared context)
- **Kimi** (`orq best`, answer `consult_0430/kimi_0430.md`).
- **Astra** (`codex exec -m gpt-6-astra`, answer `consult_0430/astra_0430.md`).
- Brief: `consult_0430/BRIEF_0430.md`, written after CD4 completed and with all closed levers listed.

Both answered **"no"** to "is there any remaining local experiment tonight with positive expected official
value", and independently gave the same reasoning on all three candidate directions:

| candidate | both said |
|---|---|
| (a) shrink submission-only calls to help NMAE | **No.** NMAE scores only the *reference's* significant genes; genes only the submission calls contribute nothing to its numerator. Also needs hidden reference info. Astra additionally corrects the 01:05 brief: NMAE is a sum of absolute errors over a **sum of absolute reference effects**, not a mean of per-gene relative errors, so small-effect genes do not automatically dominate. |
| (b) conditional dense layout where the reference is dense | **No.** The +0.009 is an *oracle* number keyed on hidden reference call count. Available proxies are weak (response magnitude ρ=0.347, v9 call count ρ=0.039), so any end-to-end gain would be privileged selection. |
| (c) rebuild-from-source reproducibility | **Secondary.** Real but modest: the published bytes already match the digest, 900/900 receipts re-audited, gates pass. Establishes a reproducibility *claim*, not submission correctness. |

## Decision
1. Complete the running `vcc_family_quality_matched_1025` H1 three-arm test under its frozen gate.
2. If `family_quality` passes, run the identical incumbent-vs-family_quality comparison on CD4 — the
   continuation **already registered** in `PREREGISTRATION_AMENDMENT.md`. Do not invent a replacement arm,
   add a variant, or relax the gate after seeing H1.
3. If H1 fails, the lever closed; stop and write the decision summary. **Do not** reopen a closed mechanism or start
   an unregistered sweep to keep the CPU busy — both consultants explicitly rejected CPU occupancy as a
   justification.
4. Leave a morning decision summary, containing: the
   family verdict against every frozen criterion, links to results/verification, the exact retained package
   path and full digest, the distinction between **verified artifact identity** and **untested source
   reproducibility**, any pending results marked pending, and an explicit recommendation.

## Why this is not idling
The instruction is to keep *useful* work going. The registered CD4 continuation is useful work with a
pre-committed gate; a fresh unregistered arm would be a search over levers that measured evidence has
already rejected. Kimi: "Once the bounded family programme and its decision receipt are complete, stop and
hand off. This is an evidence-based stopping decision after completed tests."

No submission, no package change. Safety entry `5NeeoWn96lfm23hgExcm` and package
`ce78c81c453f74676ea0e6d2e1298cccf58aa0c2d6649f930b428b62a4939dc4` unchanged.
