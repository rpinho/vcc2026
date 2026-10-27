# 2026-10-25 04:43 EDT — revision 11 scope: wrapper exit, strict decode schema, source coverage

Consulted Kimi (`orq best`) on two AUDIT10 probe results and the audit-authorization scope. Full answer
was returned in-session; the operative points are below.

## Evidence

- `rebuild_repro_1025/audit10_scratch/supplemental_results.json`:
  - `chain_nonzero_wrapper_exit`: appending `EXIT=7` after `done final_billmA_repro` still yields
    aggregator `PASS` and a PASS headline.
  - `float_decode_shape`: a decode receipt with float-valued `structure.shape` still yields
    `decode_ok: true` and comparison `status: PASS`.
  - `missing_final_receipt_run_started_log_utc`: aggregator `PASS` (renderer still PARTIAL).
- The live rebuild log ends `done final_billmA_repro` then `EXIT=0`, so the legitimate production path
  already carries the wrapper status. The current parser strips it without checking its value.

## Kimi's assessment

1. **Nonzero wrapper exit is a genuine false-PASS defect.** Require the current attempt's wrapper exit to
   be exactly `0`. Missing exit evidence is pending, not success; malformed, duplicate, conflicting, or
   nonzero exits block PASS. Bind the attempt id, exit status, and finalized log length/hash into the
   completion receipt. The startup prefix cannot establish completion.
2. **Float decode fields are a schema false-PASS**, even when numerically equal. Require
   `type(bytes) is int` and positive; validate the decode structure with the same integer/range/hash
   schema as build structure; reject booleans and strings; do not coerce.
3. **Require an independently defined source set.** Enforce `required_paths ⊆ manifest_paths`, then
   verify every listed file against the frozen snapshot, live bytes, and audit report. The required set
   is a reviewed, versioned policy independent of the candidate manifest; a singleton unrelated manifest
   must not authorize.

## Revision 11 changes

- `run_state.py` / `attempt.py`: one shared terminal parser; `DONE` requires trailing `EXIT=0`;
  `EXIT=<nonzero>` is `FAILED`; missing/duplicate/conflicting exit is not `DONE`.
- `watch_and_finish.sh`: receipt records `wrapper_exit_code`, `build_log_bytes`, and
  `build_log_sha256`; `build_completed` requires exit 0.
- `verify_final_evidence.py`: require receipt `run_started_log_utc`, `build_log_bytes`,
  `build_log_sha256`, and `wrapper_exit_code` to match the live finalized log; strict decode types.
- `compare.py`: strict integer `bytes`, strict `structure` schema, exact digest string, integer exit.
- `render_handoff.py`: require the full 26-file verification-chain source set in the revision manifest
  (subset check, not candidate-supplied), then verify every manifest source live and hash-listed.
- Regressions for `EXIT=7`, missing exit, float decode bytes/shape, bool structure, missing run-start
  stamp, singleton manifest, and omitted required source.

Nothing submitted; no package change.
