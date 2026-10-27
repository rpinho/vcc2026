# 2026-10-25 13:30Z — Astra independent audit of frozen closeout verifier

## Why this record exists

The OpenAI `vcc-closeout31` run was refused by a cybersecurity content filter
before writing a report. The Kimi fallback was slow and repeatedly
lost time to its own harness bugs. A third, fully isolated audit was therefore
launched on `gpt-6-astra` against the same frozen bytes, writing only to
`closeout_audit31_astra_scratch/`.

## Outcome

Astra was **also stopped by the OpenAI cybersecurity filter** before it wrote
`CLOSEOUT_AUDIT31_ASTRA.md`. It did, however, complete and persist all of its
evidence before the refusal. That evidence is usable and is recorded here; the
refusal is **not** a verdict and is not treated as one.

## Completed evidence (all against the frozen hashes)

- Independent semantic verification: **1,989 checks, 0 failures**, covering
  exact six-file membership, byte counts, SHA-256, run-ID binding, summary/
  receipt agreement, structure, decode identity, timing, package paths, and
  **1,800** canonical anchor comparisons (900 per archive) plus 17 live source
  row bindings.
- Cross-run isolation: distinct content and inodes for both archives, link
  count 1, no symlinks on evidence members.
- Suite controls: **54/54** mutation cases PASS, exit 0, with both genuine
  archives PASS. Every listed mutation rejected for its intended reason;
  `audit_fail`, `audit_hash_changed`, and the duplicate/unsafe AUDIT27 row cases
  correctly show `audit_changed: true`.
- Replay matrices: **14/14** audit28, **24/24** audit29, **34/34** audit30 all
  rejected.
- New adversarial probes: 27 additional cases (exact-type variants, NaN/
  Infinity, numeric strings, nested containers, duplicate keys in every
  evidence file, path aliases, symlink and hardlink topology, cross-run reuse,
  summary-shape coercion). All reject.
- TOCTOU probes: receipt swapped after hashing but before parsing, as both a
  symlink and a regular file, are caught (`archive file changed` /
  `not a private regular file`). Read-only audit hooks blocked all write and
  process events during verification of both archives.

## The one reported "false PASS" — resolved as a harness bug

`replay30 symlink_archive_ancestor` initially reported PASS. The recorded
command shows Astra passed the **real** archive path, not the `alias/...` path,
so the verifier never saw a symlinked ancestor. Re-invoking with the alias path
gives `FAIL archive directory has a symlinked ancestor`, and Astra's own
corrected case (`replay30_symlink_ancestor_correct_invocation`) reproduces that
FAIL. This is a probe-construction defect, not a verifier defect.

## Current status

- Verifier `1f1df2b704c3d7d113b31d9d7cecb1d7d914723452b0ab80d44dcecabe6661d2`.
- Suite `5efd424ed632a75863ad9cf489e430d298b659482336587a534cc4cfd207b35c`.
- No demonstrated false PASS remains on these bytes across local review, the
  audit30 replay, and this Astra evidence.
- No report file exists, so there is still **no formal PASS verdict**. The Kimi
  audit is the run expected to produce `CLOSEOUT_AUDIT31.md`.

## Decision

Preserve this evidence and do not restate it as a verdict. Keep frozen bytes
unchanged while Kimi finishes. If Kimi also fails to produce a report, the next
step is a consultation on whether the adjudicated evidence (local review +
audit30 replay + Astra's independent 1,989 checks and 54/54 controls) is
sufficient to treat the verifier as qualified, or whether a different
independent auditor must be found first.

No submission, upload, package replacement, or new scientific arm. Safety entry
`5NeeoWn96lfm23hgExcm`; package `2,595,584,000` bytes, SHA-256
`ce78c81c453f74676ea0e6d2e1298cccf58aa0c2d6649f930b428b62a4939dc4`.
