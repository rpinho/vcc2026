# 2026-10-25 13:10Z — audit31 refusal; numeric-coercion repair

## Audit31 outcome

`vcc-closeout31` did not produce a scientific verdict. Its inventory and
independent checks completed, and its replay reached C076 with 76 invocations.
The provider then ended the run with a cybersecurity content-policy flag
(`exit=1`) before writing `CLOSEOUT_AUDIT31.md`. The frozen audit log is not exported. Its recorded
SHA-256 `3b7ff2b9361b25b6177bcc7dba3bdf44834092e4318cb5789d92b0f15be0f7cd`.
A policy refusal is not a PASS or FAIL on the verifier.

The audit's own suite wrapper also crashed while parsing a NaN test fixture
because its independent loader rejected nonstandard JSON constants outside a
per-case error handler. That is an audit-harness defect, not verifier evidence.

## Additional defects found by local review

Two genuine false PASSes remained in the repaired verifier:

- a manifest file record with `bytes: 377.0` was accepted for a 377-byte file;
- `decoded_shape: [360000.0, 18533]` was accepted for the required integer shape.

Both were ordinary loose numeric equality. Direct CLI probes against scratch
archives reproduced both on the audit31 target bytes.

## Repair

Only `package_bridge_1025/verify_production_closeout.py` and
`package_bridge_1025/test_verify_production_closeout.py` changed.

- Manifest byte counts must be strict `int`.
- `decoded_shape` is compared with type-exact `same_value`.
- Decode/receipt byte-count agreement uses `same_value`.
- JSON parsing rejects nonstandard constants (`NaN`, `Infinity`) through
  `parse_constant`.
- The canonical 900-name set is now constructed explicitly from
  `group_{0..2}_{000..299}.json`, independent of glob contents.
- The suite adds `manifest_byte_float`, `manifest_bytes_nan`, and
  `decode_shape_float` regressions.

## Verification

- Repaired verifier:
  `1f1df2b704c3d7d113b31d9d7cecb1d7d914723452b0ab80d44dcecabe6661d2`.
- Repaired suite:
  `5efd424ed632a75863ad9cf489e430d298b659482336587a534cc4cfd207b35c`.
- Mutation suite: **54/54**, exit 0.
- Both genuine archives: **PASS**.
- Direct probes reject float, NaN, boolean, and string manifest byte counts.
- `decode_shape_float` rejects with `FAIL decode shape`.

## Decision

Freeze these bytes and run the independent audit with Kimi,
because OpenAI refused the audit prompt. The audit must not claim a verdict from
the refused run. Preserve both the refusal and this repair rationale.

No submission, upload, package replacement, or new scientific arm. Safety entry
`5NeeoWn96lfm23hgExcm`; package `2,595,584,000` bytes, SHA-256
`ce78c81c453f74676ea0e6d2e1298cccf58aa0c2d6649f930b428b62a4939dc4`.
