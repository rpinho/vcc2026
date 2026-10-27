# 2026-10-25 12:50Z — audit30 boolean type-confusion repair

## Trigger

Independent audit 30 returned **FAIL** on the closeout verifier. The archive
evidence remained sound: both genuine archives passed, all 14 audit28 probes and
all 24 audit29 probes rejected, and the then-current mutation suite passed
44/44. Audit30 nevertheless demonstrated six false PASSes caused by Python
equality treating `True == 1` and `False == 0` as equal:

- `summary_boolean_as_integer`
- `summary_zero_as_boolean`
- `receipt_meta_schema_boolean`
- `decode_meta_schema_boolean`
- `receipt_mismatch_boolean`
- `decode_mismatch_boolean`

A seventh probe, `log_invalid_utf8_rehashed`, failed closed with exit 1 but
emitted a traceback instead of a structured `FAIL`.

The frozen audit30 report is preserved unchanged at
`package_bridge_1025/CLOSEOUT_AUDIT30.md` (81,993 bytes).

## Repair

Only `package_bridge_1025/verify_production_closeout.py` and
`package_bridge_1025/test_verify_production_closeout.py` changed.

- Added recursive `same_value()` with exact type equality before value equality.
- Routed exact summary comparison, fixed metadata, input hashes, structure,
  archived structure, and state timing through `same_value()`.
- Required `type(value) is int` for receipt/decode mismatch counts and
  receipt/archived receipt counts.
- Caught `UnicodeDecodeError` while reading the production log and emitted
  `FAIL production log unreadable` without a traceback.

No production receipt, archive, log, manifest, AUDIT27 content, prior report,
or safety package was modified.

## Verification

- Repaired verifier SHA-256:
  `3ee73b3d301e05bbdaab9c345cd89f4774fac452a46064e7a3b30a230b673ed6`.
- Repaired suite SHA-256:
  `62958351cbc758da75c665adad43f93f376549de297aeca7b870f9fba2aefa51`.
- Mutation suite: **51/51**, exit 0, including both genuine archives and the six
  boolean probes.
- Both genuine archives: **PASS** in place.
- Audit30 regeneration replay: **34/34** probes reject for intended reasons.
- Zero unexpected PASSes, zero unstructured rejections, and frozen inputs
  unchanged before and after the replay.

Replay evidence is under
`package_bridge_1025/closeout_audit30_replay/`; the runner regenerates each
payload from genuine archives, preserving FIFO, hardlink, symlink, and
trust-root wrapper behavior instead of copying malformed trees.

## Decision

Freeze these repaired bytes and launch a fresh independent audit,
`vcc-closeout31`, rather than treating audit30's replay as authorization.
Audit31 must replay all audit28, audit29, and audit30 probes, independently
re-check both genuine archives, inspect the expanded suite, and add new
adversarial cases. A demonstrated false PASS is FAIL.

No submission, upload, package replacement, or new scientific arm. Safety entry
`5NeeoWn96lfm23hgExcm`; package `2,595,584,000` bytes, SHA-256
`ce78c81c453f74676ea0e6d2e1298cccf58aa0c2d6649f930b428b62a4939dc4`.
