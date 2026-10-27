**Rank: 1 → 4 → 2 → 3. The defensible work before 08:00 is verification plus decision summary.**

1. **Finish AUDIT27, then rerun the production bridge and independent decode.** Highest value: close the observed persisted-evidence defect on the exact audited revision. The earlier 900/900 result does not verify revision 25’s serialization fix.

   **Acceptance gate:** AUDIT27 gives an unqualified PASS bound to manifest `72dbe540bc88a03e8f6f2dc16a8063a1d8b6cac99f0e7fb7d9d007248e5ae710`; all live production sources still match it. Under a fresh run ID, bridge and fresh independent decode both pass all existing gates with 900/900 matches. Read the newly persisted bridge receipt back from disk: require exact equality with the returned receipt, including boolean `output_directories_stable: true`, consistent run/state/decode identities, and unchanged safety-package size and SHA-256. Existing disk guards remain binding.

   **Failure condition:** any audit qualification, source drift, missing/stale/unequal receipt, identity mismatch, incomplete count, decode failure, or resource-guard failure. If unfinished at 08:00, report **pending**, never PASS; retain the safety entry unchanged.

2. **Option 4: a concise evidence decision summary.** Record exact audited revision, manifest, fresh run ID, receipt locations, completed verdicts, unresolved blockers, and the unchanged safety identity. Distinguish historical results from newly verified results.

3. **Option 2: reject for tonight on present evidence.** The retained-record reconciliation already found surrogate mismatches and no defensible transfer bound. Reopen only upon identifying a previously missing exact matched observation; another fit/bootstrap of the same adaptations cannot supply it.

4. **Option 3: reject.** No registered scientific opening or acceptance rule supports another arm or sweep.