# 2026-10-25 10:40Z — AUDIT24 FAIL is a brief path-scope error, not a revision-23 defect

## Verdict read

AUDIT24 reports **FAIL** on revision 23. The report's only blocking finding is
that the brief told the auditor to hash `audit23_scratch/sources/`, which
correctly still holds the revision-22 snapshot, while the revision-23 manifest
and its matching snapshot are under `audit24_scratch/sources/`.

The brief contains two stale references copied from the AUDIT23 template:

- "hash every frozen copy under `audit23_scratch/sources/`"
- "all writes belong under `package_bridge_1025/audit23_scratch/`"

The manifest path and digest in the same brief are revision 23's, and the
revision-23 freeze script created `audit24_scratch/sources/`. The stale lines
therefore create an impossible identity requirement: no directory can contain
both the old revision-22 timestamp suite and the new revision-23 suite at the
same pathname.

## What the auditor independently established

- **17/17 live** and **17/17 `audit24_scratch/sources/`** files matched the
  revision-23 manifest before and after the audit.
- The revision-23 behavioral timestamp suite rejects both required
  counterexamples: delayed capture immediately before result construction and
  an unconditional entry exception.
- The actual decoder emits entry time `09:59:00Z`, has a single wall-clock
  sample before output/input access, matches 900/900, leaves the package
  unchanged, and cleans up temporary outputs.
- Independent synthetic probes found no false qualified acceptance through the
  revision-23 launcher/orchestrator.
- All earlier fixes were preserved; exactly one file differs from revision 22,
  the timestamp suite; decoder bytes are unchanged.
- The auditor explicitly says the substantive result "does not authorize
  substituting a different frozen directory for the explicit identity
  requirement" and recommends correcting the brief.

## Decision

1. Preserve AUDIT24 as **FAIL**. Do not edit or relabel it.
2. Do not treat the revision-23 implementation as authorized yet. A PASS must
   name the exact frozen bytes in the correct directory.
3. Issue AUDIT25 with a corrected brief that names
   `audit24_scratch/sources/` for identity and `audit24_scratch/` for writes.
   Keep the manifest digest unchanged; no source changes are needed.
4. Preserve `audit23_scratch/` as the historical revision-22 snapshot. Do not
   copy revision-23 files over it.
5. Keep the safety package unchanged and do not submit.

## Why this is not test-shopping

The requested correction changes only the location that the auditor was told to
inspect, not the implementation, the manifest, the tests, or the pass criteria.
The corrected brief binds the same manifest digest and requires the same
independent adversarial probes. AUDIT24's report and its evidence remain in the
record.

## Safety

Nothing submitted, uploaded, or replaced. Safety entry
`5NeeoWn96lfm23hgExcm`; package SHA-256
`ce78c81c453f74676ea0e6d2e1298cccf58aa0c2d6649f930b428b62a4939dc4`
(2,595,584,000 bytes) unchanged; no production decode launched.
