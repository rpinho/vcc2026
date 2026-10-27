# 2026-10-25 13:35Z — suite relocates the repo but `PACKAGE` stays bound to import-time path

## Finding

`verify_production_closeout.py` computes `PACKAGE` once at import time:

```python
SOURCE_ROOT = ROOT.parent
PACKAGE = str(SOURCE_ROOT / 'model_billmB/out/final_billmA.vcc')
```

`test_verify_production_closeout.py` relocates the verifier into a temporary tree and rebinds
`module.ROOT`, `LIVE_ROOT`, `SOURCE_ROOT`, `AUDIT`, `SAFETY`, `SAFETY_BYTES`, and `SAFETY_SHA`,
but it never rebinds `PACKAGE`. In a relocated fixture the expected package path therefore stays
pointed at the fixture copy while the receipt still records the original absolute path, so
`receipt package path` and `state package path` fail.

This was reproduced directly on the Kimi fixture:

```
import-time PACKAGE: .../closeout_audit31_scratch/suite_repo/model_billmB/out/final_billmA.vcc
receipt package   : .../vcc2026/model_billmB/out/final_billmA.vcc
equal: False
```

## Impact

- **Not a false PASS.** It is a false *rejection*, and only under repo relocation.
- The frozen suite in the real production layout still passes **54/54** — confirmed
  independently by Astra (`suite.stdout`, exit 0) and by Kimi re-running the exact frozen suite
  in the production layout (`suite_production.stdout`).
- Production verification is unaffected: in place, `SOURCE_ROOT` is the real repo, so `PACKAGE`
  matches the receipt exactly and both genuine archives PASS.

## Decision

Record the defect. Do **not** change frozen verifier/test bytes during an live audit; the audit
must complete on the hashes it was given. The clean fix, for a future revision, is either to
compute `PACKAGE` lazily inside `main()` from the current `SOURCE_ROOT`, or to have the suite
rebind `module.PACKAGE` alongside the other module constants.

Any auditor that reports this as a verifier FAIL is misattributing a fixture-environment artifact,
and any auditor that reports a 51/54 suite result without checking the production layout has not
finished the job.

No submission, upload, or package replacement. Verifier `1f1df2b7…61d2`, suite `5efd424e…b35c`,
safety package 2,595,584,000 bytes `ce78c81c…9dc4`.

---

## ERRATUM (2026-10-25 14:05Z, after independent reproducibility check)

The independent check (`closeout_repro_check_1025/REPRO_CHECK.md`, C11-C15) found that my
**51/54 figure in this record is wrong for the frozen suite**:

- Executing the **frozen** suite bytes (`5efd424e…b35c`) with a relocated import path gives
  **52/54**, not 51/54. Only `old_baseline` and `new_baseline` fail.
- The **51/54** number I originally quoted came from Kimi's scratch fixture, whose suite differs from
  the frozen bytes: it changed `module.LIVE_ROOT = HERE` to `module.LIVE_ROOT = root`
  (fixture suite SHA-256 `7e799008…330a`). That additional change also makes `audit_hash_changed`
  fail, because the temporary root lacks the live source files needed for the expected AUDIT27
  hash-row diagnostic.

So this record **conflated two distinct effects**: the genuine import-time `PACKAGE` binding defect
(52/54 under relocation) and a separate fixture modification (the extra 51st failure).

**What still holds:** the `PACKAGE` defect is real, it is a **false rejection** under relocation
rather than a false PASS, and rebinding `PACKAGE` restores exit 0 on the same genuine archive
(verified as C13 in the reproducibility check). Production verification is unaffected.

**What is withdrawn:** the claim that a relocated 51/54 result is "solely attributable to `PACKAGE`",
and the 51/54 count as a property of the frozen suite. The correct frozen-suite relocation figure is
**52/54**, failing exactly the two baseline controls.
