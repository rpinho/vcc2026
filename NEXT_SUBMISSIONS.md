# Next submissions — 2026-10-25

**Recommendation: do not submit.** The current last entry is already the safety package, and no new
candidate has passed the qualification gates. This file is guidance for Ricardo's manual decision; no
command here has been executed.

## Verified state as of 2026-10-25 09:55 EDT

Everything below was re-checked today, after the closeout verifier was hardened and independently
audited. **Precision matters here:** the gate predicates were re-evaluated against retained receipts,
not by rerunning `vcc prep`, redoing the historical 900 dense-matrix computations, or performing a
fresh decode. Those underlying artifacts are dated Oct 22 and are being *re-checked*, not regenerated.

- **Closeout verifier: `CLOSEOUT_AUDIT31.md` PASS.** Verifier `1f1df2b7…61d2`, suite `5efd424e…b35c`,
  mutation suite 54/54, both genuine archives PASS.
- **Full gate chain re-evaluated:** bridge audit gate PASS, receipt gate PASS
  (`complete, current-run 900/900`), decode qualification PASS, and package gates
  **G1-G8 all PASS** (`check_gates.py final_billmA`, `ALL GATES PASS`, exit 0). All eight were
  independently reproduced from primary artifacts on 2026-10-25
  (`closeout_repro_check_1025/REPRO_CHECK.md`), with the caveat above.
- **No certified exact local/official pair exists in the retained inventory.** A bounded read-only
  inventory found **0 EXACT-MATCHED** pairs out of 67 records (0 exact, 0 approximate, 23 unverified,
  44 no-official) under its four-field definition of recipe identity. See
  `matched_inventory_1025/INVENTORY.md`.
  This is a certification result over a bounded scope, not a proof that no such pair could exist
  anywhere. The operational decision — no promotion rule is justified on this evidence — stands, but
  it is a judgement, not a number. A looser "procedural" reading would make the v5a/v5b LOCO ports
  arguably recipe-matched, though still 79-target ports with `external_official_label: null` and
  therefore still not usable as official transfer observations.
- **Submission state (limited):** `vcc status 5NeeoWn96lfm23hgExcm` → `published`, `mode: FINAL`,
  panel `vcc2026-test-1`. `vcc whoami` → `can_submit: true`, deadline **Nov 05 23:59 UTC**.
  These two commands do **not** establish that no newer entry exists (status queries one ID) or that
  daily allowance is unused (that uses a separate `get_limits` endpoint, only exercised at submit
  time). They are consistent with a free slot; they do not prove one.

This does not change the recommendation: the safety entry is already last, is fully qualified, and
there is no evidence-backed reason to replace it.

## Current last entry (leave in place)

- Entry: `5NeeoWn96lfm23hgExcm` — "final billmA v9 layout r2"
- Package: `model_billmB/out/final_billmA.vcc`
- Size: `2,595,584,000` bytes
- SHA-256: `ce78c81c453f74676ea0e6d2e1298cccf58aa0c2d6649f930b428b62a4939dc4`
- Status: published; current last entry. No resubmission is needed.

## New candidates

None qualified. The score search is closed: the family lever failed its frozen budget gate
(`Δ 0.051027 > 0.020`), and prior Kimi/Astra consultations recommended stopping the family/amplitude/pool/
source search. The root `rsi/NEXT_SUBMISSIONS.md` records no qualified RSI candidates.

## Fallback command — NOT EXECUTED

Only if the last-entry state is lost and Ricardo explicitly decides to restore the safety package, this is
the exact form. Re-running it when the safety entry is already last is unnecessary and consumes the daily
allowance. Do **not** run it automatically.

```sh
vcc submit \
  model_billmB/out/final_billmA.vcc \
  -m "final billmA v9 layout r2" \
  -d "Bill M response-MLP blend on public source transfer, v5b amplitude, v9 call layout, panel-agnostic donor banks. No improvement claim." \
  --wait
```

## Evidence and open verification

- Source-reproducibility rebuild is **measured PASS**: 900/900 group receipts, 900/900 independent
  monitor, 900/900 preserved-NPZ rehash, 12/12 dense samples, identical structure, and a rebuilt package
  byte-identical to the published one (`ce78c81c…9dc4`). See `rebuild_repro_1025/final_receipt.json`.
- The rebuild evidence chain is **PASS** at its newest independent audit, `rebuild_repro_1025/AUDIT13.md`
  (revision 13, manifest SHA-256 `b99d5b2d…2735fc`, 27/27 live and frozen sources matched). The
  earlier AUDIT8 defects were repaired by revisions 9-13; AUDIT8 is superseded, not an authorization.
  An older note in this file said PARTIAL; that was stale as of 2026-10-25 05:46 EDT.
- Bridge revisions 1–24 are disqualified or superseded by independent audits. **Revision 25 is
  authorized by AUDIT27 PASS** (manifest SHA-256
  `72dbe540bc88a03e8f6f2dc16a8063a1d8b6cac99f0e7fb7d9d007248e5ae710`, 17/17 live and frozen sources
  matched), and revision-25 production passed 900/900 with a fresh independent 900/900 decode under
  two clean run IDs. `CLOSEOUT_AUDIT31.md` then returned **PASS** on the frozen closeout verifier
  (`1f1df2b7…61d2`, suite `5efd424e…b35c`, 54/54). The earlier "revision 16 / AUDIT16" line was stale.
- Nothing was submitted, and the published package was not modified.
