# 2026-10-25 07:52Z — rebuild chain PASS; bridge revision 13 frozen, revision 14 adds archive pin

## Rebuild: PASS (this is now a measured result, not a claim)

The source rebuild of `final_billmA` completed and the whole finish chain exited 0:

| gate | result |
|---|---|
| build + `panel_prep` decode | PASS, exit 0 |
| group verification | PASS 900/900 |
| preserved NPZ rehash | PASS 900/900, zero-byte count 0 |
| retained-evidence repair | PASS, 0 objections (866 intact + 34 replaced) |
| dense sample arrays | PASS 12/12 |
| comparison vs published | PASS — identical 2,595,584,000 B / `ce78c81c…9dc4` |
| independent final-evidence aggregation | PASS, 0 failures |

Conclusion: the published artifact is **container-identical** to an independent rebuild, and its
per-group evidence is now intact and hash-matched to the archived receipts. The published package was
never written; it was rehashed and matches.

## Retention defect: what it actually was

The 03:30 checkpoint understated this as "51 zero-byte files". It was **85 corrupt of 900**: 51
zero-byte plus **34 truncated-but-nonempty** files that do not parse as NPZ at all (`BadZipFile`). 815
were intact. `preserve_state.json` reported COMPLETE because it tested inode independence, not content;
that record is not evidence of integrity and is not cited as such anywhere in the decision summary.

Repair was evidence-driven: each block was reconstructed from the independently decoded H5AD, then
required to match the archived `counts_sha256` exactly before any replacement. For the 34 truncated
files the retained copy *cannot* be evidence (it does not parse), so it could not veto a recovery that
already matched the archive. That rule is now explicit in `repair_retained_npz.py`.

## Finish-chain bugs found and fixed (fail-closed, no scientific gate relaxed)

1. `watch_and_finish.sh` and `finalize_after_receipt.sh` did `import attempt` without the script
   directory on `sys.path`. The finalizer therefore rejected a *valid* receipt with a traceback and
   would have waited forever. Fixed by explicit `sys.path` insertion.
2. `watch_and_finish.sh` required the build log's literal last line to be the completion marker, but the
   rebuild wrapper appends `EXIT=<code>` after it. `attempt.py` and `run_state.py` already ignore that
   wrapper line; the gate now applies the identical rule.
3. Three off-contract samples (`group_{0,1,2}_004.npz`, written 01:18 by an earlier sampler revision)
   made the 12-name sample gate report FAIL for the wrong reason. They are genuine evidence, so they were
   quarantined to `sample_npz_offcontract/` with a README rather than deleted.

## Bridge: still NOT authorized

AUDIT11 FAIL, and AUDIT12 (whose agent died on a content flag without writing its report) reproduced
real false acceptances from its completed probe logs. The most serious: **the receipt gate accepted a
receipt that merely asserted 900/900 about a 64-byte synthetic file.** It trusted the receipt's own
claims instead of recomputing anything.

Revision 13 closes that class:
- the gate recomputes the named package's size and SHA-256 itself;
- counts must be strict integers (`900.0`, `true`, `"900"` all rejected) and duplicate JSON keys fail;
- the 900 claimed group digests must cross-bind to `group_*.json` files the gate reads itself;
- receipt and run state may not share a pathname;
- group-index probes no longer leave a PASS state behind a FAIL receipt;
- malformed group evidence (array/null/scalar/invalid UTF-8) fails closed with a clean reason.

All 30 AUDIT12 package probes were re-run against the patched bridge: every failing case now fails
closed and **no failing case emits PASS evidence of either kind**.

Suites: bridge 7/7, findings 45/45, audit gate 43/43, receipt gate 20/20 = **115/115**.

## Revision 14: the hole I found while writing the audit spec

While drafting AUDIT13's spec I asked whether the *new* cross-binding could itself be spoofed, and it
could: a caller pointing `--archive-dir` at a directory it controls could supply 900 self-consistent
fabricated group receipts and satisfy the binding. That proves nothing about the published package.

Production was still protected by the safety SHA/size anchor, so this was not a live false-PASS path —
but it is the same "trust the input" mistake, so it is closed rather than documented: the archive is now
pinned to `model_billmB/artifacts/panels`, and an out-of-tree archive is refused unless a synthetic
fixture explicitly passes `allow_unpinned_archive=True`. Regression case added; the production control
(real package + real archive) still passes with zero failures.

Revision 14 is frozen at `audit14_scratch/revision14_manifest.json` and supersedes 13.

## Status and next gate

- AUDIT13 is auditing revision 13. Because revision 14 changes only `check_bridge_receipt.py`, AUDIT13's
  verdict on 13 remains informative but **cannot authorize the bridge**: the bridge must be audited at
  the exact revision that would run.
- Next: when AUDIT13 reports, read its verdict. If PASS, launch AUDIT14 against the revision-14 manifest.
  If FAIL, treat the findings as real, fix, and re-freeze. Do not start the bridge on any non-PASS.
- Nothing submitted. Safety entry `5NeeoWn96lfm23hgExcm`.
