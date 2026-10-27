# Decision — record rebuild identity; gate bridge on AUDIT11

Date: 2026-10-25 07:45 UTC / 03:45 EDT
Scope: verification only; no score search, no package change, no submission.

## Decision
1. Treat `model_billmB/out/final_billmA_repro.vcc` as a **measured byte-identical reproduction**
   of the published safety package: 2,595,584,000 bytes, SHA-256
   `ce78c81c453f74676ea0e6d2e1298cccf58aa0c2d6649f930b428b62a4939dc4`.
   The 900-group build and the independent group monitor both completed PASS. This is a source
   reproducibility result, not a new score result and not a promotion.
2. Do **not** start the published-package bridge from an older audit. A supervisor waits for
   `package_bridge_1025/AUDIT11.md`, requires `check_audit_gate.py` to accept a clean PASS and all nine
   live source hashes, and only then execs `run_after_rebuild.sh`. The runner independently repeats the
   audit check and the receipt gate; neither layer can widen the other.
3. Keep the score search closed. The family lever failed its frozen `budget_match` gate
   (0.051027 > 0.020), and no new qualified candidate exists. The safety entry stays last.

## Why
The rebuild's container identity is now measured, but the verification chain still needs an independent
audit that can find false-pass paths. The bridge's purpose is evidence integrity; running it before that
audit would turn an unverified verifier into a PASS claim.

## Safety
Nothing submitted. Published package unchanged. No credentials, pods, public-repository writes, or
`OFFICIAL_SCORES.md` edits. The decode and audit run concurrently on local CPU only.

## Follow-up
- Wait for `panel_prep` decode, then the watcher's compare/preserved/sample/evidence chain and the finalizer.
- If AUDIT11 is PASS, the armed supervisor starts the bridge. If FAIL, freeze revision 12 and audit it;
  do not start the bridge on a partial or contradictory report.

---

## Correction (07:52Z, same night)

Two statements above are now outdated and must not be relied on:

- **AUDIT11 did not authorize the bridge.** It returned **FAIL** on a real false acceptance. The
  supervisor that was waiting on it exited with `SUPERVISOR_EXIT=10` without starting
  any bridge run. Point 2's *mechanism* is still correct (no bridge without a clean PASS at the exact
  running revision); only the expectation that AUDIT11 would supply it was wrong.
- **The rebuild's evidence preservation was not PASS at 07:45.** The 900-group monitor and the build
  identity were PASS, but the retained NPZ set was 85/900 corrupt (51 zero-byte, 34 truncated). It has
  since been reconstructed, verified against archived hashes, and repaired: see
  `decisions/20261025_0752_rebuild_pass_and_bridge_rev13_14.md` for the measured PASS and the finish-chain
  bugs that were fixed to reach it.

No package was changed and nothing was submitted.
