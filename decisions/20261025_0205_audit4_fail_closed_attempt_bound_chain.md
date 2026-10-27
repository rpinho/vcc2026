# Decision — close AUDIT4 FAIL with one attempt-bound verification chain

- **When:** 2026-10-25 ~02:05 EDT (06:05Z)
- **Trigger:** `rebuild_repro_1025/AUDIT4.md` returned `Verdict: **FAIL**` on the verification chain, with
  five open findings (D2, D3, D5, D7, D8). The rebuild itself was healthy; the FAIL concerned whether a
  false PASS could reach the decision summary.
- **Decision:** patch all five findings at the root before allowing a final decision summary, then launch a fresh
  independent audit (revision 7) against the frozen patched source. Do not claim PASS from AUDIT4, and do not
  let the stale revision-6 audit queue run.
- **Why not defer:** the headline in the decision summary is a scientific claim of source reproducibility. A chain
  that can report PASS from a cached receipt after a later failure would make that claim unsupported.

## What changed
1. `attempt.py` — one immutable record per watcher run: UUID token, exact log-prefix bytes/SHA, start time.
   `check()` requires the prefix to still hash and the latest segment to be `done`.
2. All four evidence producers stamp `attempt_id`; the aggregator requires record + receipt + evidence to
   agree, plus the log-prefix and attempt-record digests.
3. `compare.py` — attainable-structure bounds closing D5.
4. `verify_samples.py` — canonical CSR and pre-densification range checks closing D8.
5. `watch_and_finish.sh` — bounded child waits closing D3.
6. Finalizer and renderer bind to the attempt record, not the mutable token file.

## Verification
120 regression cases across 11 suites, zero failures. The exact AUDIT4 reproducers (failed restart after a
cached success; old writer adopting a newer attempt; four impossible structures; duplicate int64 overflow)
now all fail closed.

## Constraints honored
No submission, no package modification, no credential or pod use.
The published package remains `2,595,584,000` bytes / `ce78c81c…`; safety entry remains `5NeeoWn96lfm23hgExcm`.

## Residual
D3 is now bounded in the watcher, but this is a liveness improvement, not a proof against every possible
hang inside a Python child. The revision-7 auditor is asked to state any remaining bound honestly. Evidence
files are still shared pathnames; the token makes cross-run mixing detectable rather than impossible.
