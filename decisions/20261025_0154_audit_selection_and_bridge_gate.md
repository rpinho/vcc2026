# Decision — audit provenance and bridge execution gate (2026-10-25 05:54Z)

**Decision.** Decision-summary audit selection is by highest numeric `AUDIT*.md` revision. The published-package
bridge may not run without a current bridge `AUDIT2.md` whose first line is exactly `Verdict: **PASS**`.

**Why.** A hard-coded AUDIT3 preference would have outlived its evidence after AUDIT4/AUDIT5. Separately,
a bridge audit that returned FAIL could previously arrive after the watcher had already begun a multi-GB
decode. Both are false-evidence risks, not cosmetic documentation issues.

**Evidence.** Renderer suite 26/26; bridge fixture suite 6/6; audit-gate suite 7/7; shell and Python syntax
pass. Live build unaffected: 259/900 at checkpoint, monitor 263/263 exact, 0 races/errors.

**Boundaries.** No score lever reopened, no package write, no submission. Published SHA remains
`ce78c81c453f74676ea0e6d2e1298cccf58aa0c2d6649f930b428b62a4939dc4`.
