# 2026-10-25 11:55Z — clean revision-25 production PASS; independent closeout audit

## Result

The clean end-to-end rerun authorized by AUDIT27 completed under fresh run ID
`20261025T114827Z-4b7be0bad0a0d1f3d7895805aefbffb1`.

- Bridge: **900/900**, zero mismatches, exact structure, safety hash match.
- Persisted receipt: `PASS`, boolean `output_directories_stable: true`.
- Returned-summary equality: every returned field agrees with the persisted
  receipt.
- Fresh independent decode: **900/900**, shape `[360000, 18533]`, nnz
  `2010253470`, package unchanged before and after.
- Bridge exit 0; decode qualification exit 0; clean production wrapper exit 0.

Evidence is preserved with a hash manifest under
`package_bridge_1025/production_runs/20261025T114827Z-4b7be0bad0a0d1f3d7895805aefbffb1/`.
The read-only closeout verifier passes on the archive, and its mutation suite
passes **11/11**.

## Why this is a real result and not an exit-code inference

The bridge process, receipt gate, fresh-decode qualification, wrapper equality
check, preserved-manifest rehash, and independent semantic checks all agree on
the same fresh run ID and the same safety package identity. The earlier
revision-24 production PASS is retained only as historical evidence and does
not authorize revision-25 bytes.

## Remaining verification gap

The verifier mutation suite hardcodes the earlier archive
`20261025T112331Z-2c3a6c7bfedc2642a44b9f3d3b06b2b2` as its fixture. The new
archive passes the verifier directly, but the mutant suite does not yet
exercise it. This is a coverage limit, not a demonstrated evidence defect.

The verifier also iterates the manifest's file mapping rather than requiring a
fixed set of six evidence files. An empty or reduced manifest could therefore
omit archive coverage. The verifier's semantic checks still require the named
receipt, state, decode, returned-summary, and run-ID files, but the manifest
completeness rule needs independent testing.

## Decision

Launch a bounded, read-only independent audit against the frozen sources.
It must:

1. Reconstruct the archive without calling the verifier.
2. Probe empty summaries, omitted/extra manifest entries, run-ID mismatch,
   log/run-ID mismatch, and extra unlisted files.
3. Run the supplied suite and state its fixture-coverage limit.
4. Verify the verifier does not write production evidence and that it binds
   AUDIT27's 17 sources, the safety package, and the archive.
5. Test cross-run isolation between the two preserved archives.

Do not launch a new scientific arm, sweep, package replacement, upload, or
submission. Do not edit production evidence to make an audit cleaner. A
demonstrated false PASS is `FAIL`; sound verification with incomplete tests is
`PASS WITH FIXES`.

## Safety

Nothing submitted, uploaded, or replaced. Safety entry
`5NeeoWn96lfm23hgExcm`; package `2,595,584,000` bytes, SHA-256
`ce78c81c453f74676ea0e6d2e1298cccf58aa0c2d6649f930b428b62a4939dc4`.
