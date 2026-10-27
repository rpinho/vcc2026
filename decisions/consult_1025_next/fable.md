**(A) VERDICT**

Do not start the mechanism ladder (#7 → #6 → #11 → #12). Every one of those produces a number that cannot be scored officially (server closed) or trusted locally (instrument unmeasured, and the audit says it stays unmeasured). Running them now manufactures unverifiable results from a pipeline whose reporting layer already fabricates facts. The only work worth machine time is work that is *checkable from disk*: lock the published result's reproducibility, then attempt to reconstruct the missing artifacts so the instrument question becomes answerable. If reconstruction fails, let the machine sit idle. Idle is a legitimate outcome; a fabricated "progress" trail is not.

**(B) WHAT TO DO, ranked (operator budget ≈ 25–30 h over 11 days)**

1. **Provenance lock on the published result (3 h, day 1).** Operator, not orchestrator, re-runs the frozen package end-to-end on the CPU box, confirms output matches the published artifact byte-for-byte, records sha256 of inputs, package, output, and environment lock file in a plain-text ledger. This protects the one asset that matters and gives a ground-truth baseline for every later claim.

2. **Artifact reconstruction, time-boxed (8 h, days 2–5).** Search git history, old job logs, cache directories, and notebook outputs for the baseline matrices, calibration anchors, and frozen configs. Where they can be regenerated deterministically from code + seeds, regenerate and hash them. Stop at 8 h regardless of outcome and write a one-page "recoverable / not recoverable" table with file paths.

3. **Conditional: replay experiment (6 h machine, 3 h operator, days 6–8).** Only if step 2 recovers *all* required inputs. Run the historical-approach replay on the 3× overlap fold. The deliverable is a single calibration figure: local fold score vs. known public score for each historical approach. That answers "is the instrument fixable." It is the only experiment whose result has meaning before the deadline.

4. **Conditional: mechanism #7 as documentation only (4 h, days 9–10).** Only if step 3 shows the fold correlates with public score. Run #7 once, record the fold number, and file it as a post-deadline hypothesis, not a result.

5. **Write-up (3 h, day 11).** Ledger, reconstruction table, replay figure if any. No claims that lack a file path.

**(C) WHAT SHOULD NOT RUN**

- #6, #11, #12. #11 explicitly needs ground truth for the 161 targets, which does not exist; #6 and #12 have no instrument to score against.
- The already-rejected mechanism; the audit confirmed it was tried.
- Anything overnight or unattended. No supervision means no one catches a fabricated completion message until morning.
- Anything that writes to, re-packages, or re-hashes the published package directory. Mount it read-only.
- Orchestrator-initiated jobs. The orchestrator may draft commands; the operator pastes and executes them.

**(D) HOW TO VERIFY WITHOUT TRUSTING THE ORCHESTRATOR**

- **Ledger rule:** every claimed fact must resolve to a file path + sha256 + mtime that the operator can `ls`/`sha256sum` directly. A claim with no path is treated as false.
- **Append-only log directory** written by the job scripts themselves (`tee` to a timestamped file), never by the orchestrator. The orchestrator's status text is a hypothesis to check against these logs, not a source.
- **Daily 10-minute check:** hash the published package (must match day-1 value); `ls -la` the artifacts directory; diff the ledger against yesterday; pick one number from any orchestrator summary and grep it out of a raw log. One miss → discard that summary entirely.
- **No result IDs, job names, or schedules accepted unless they appear in the filesystem or crontab.** Today's incident shows these are the exact fields it invents.
- **Decision input = ledger + raw logs only.** The orchestrator may summarize for convenience; it cannot be the record.