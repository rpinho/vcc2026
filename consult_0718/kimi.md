**The defensible priority before 08:00 ET is verification plus decision summary.** Rank: **1 → 4 → 2 → 3.**

1. **Finish AUDIT27, then rerun the audited production bridge and fresh independent decode.** This closes a demonstrated evidence defect. The previous 900/900 result cannot establish that revision 25 persists the corrected receipt.

   **Acceptance gate — all required:**
   - AUDIT27 explicitly passes revision 25, bound to manifest SHA-256 `72dbe540bc88a03e8f6f2dc16a8063a1d8b6cac99f0e7fb7d9d007248e5ae710`; execution sources match the audited hashes.
   - Production execution succeeds with a fresh supervisor-generated run ID; persisted receipt and run state pass their validators for that ID.
   - Persisted receipt equals the returned receipt and contains boolean `output_directories_stable: true`.
   - Bridge verifies exactly 900/900 groups, zero mismatches, exact identities/order/structure, and unchanged package bytes matching the stated safety size and SHA-256.
   - A newly executed independent decode starts no earlier than this bridge, passes its qualification checker, and confirms 900/900 with zero mismatches against that same package identity.

   **Failure condition:** any audit failure, source drift, missing/unequal/stale evidence, failed validator, mismatch, or incomplete execution blocks acceptance. Disk failure or reaching 08:00 unfinished means incomplete verification; retain the safety entry and report the gap explicitly.

2. **Option 4: a concise evidence decision summary.** Record audited revision/manifest, run ID, receipt locations, qualified result or exact outstanding failure, unchanged safety identity, and closed scientific levers.

3. **Option 2: defer before 08:00.** Retained-artifact reconciliation could support a later diagnostic, but no established transfer-aware acceptance rule currently makes it decision-relevant. Reject using it to justify promotion tonight.

4. **Option 3: reject.** A new arm or sweep lacks a registered decision basis and would reopen closed search.

No modifications made.