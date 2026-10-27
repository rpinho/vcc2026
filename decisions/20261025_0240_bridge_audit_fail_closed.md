# Decision — close the bridge AUDIT FAIL before any decode

- **When:** 2026-10-25 ~02:40 EDT (06:40Z)
- **Trigger:** `package_bridge_1025/AUDIT.md` returned `Verdict: **FAIL**` with nine open findings (F1–F9).
- **Decision:** patch all nine, extend the suites, and queue a fresh independent audit (`AUDIT2.md`) against
  a frozen revision-2 manifest. The bridge stays blocked until that report passes with matching source hashes.
- **Why:** the bridge is the only independent check that the *published* artifact decodes to the claimed
  group content. A bridge that can report PASS while having decoded a different package (F2), or that can
  overwrite the package through its own output path (F1), produces evidence that cannot support a claim.

## Most important finding
F1 was not theoretical. The first patch attempt checked path safety but still wrote the FAIL receipt to the
requested output path, which the fixture had set to the package itself. `f1_input_bytes_preserved` caught the
package being rewritten. Path validation now runs before any write at all, and an unsafe request writes its
receipt under `tmp/` instead.

## Verification
Original 6/6, audit gate 8/8, receipt gate 11/11, audit findings 15/15 — 40 cases, zero failures. Fixtures
exercise the real verifier and the real gate scripts on tiny synthetic packages; the 2.6 GB artifact was
never opened, hashed, or decoded by any test.

## Constraints honored
No submission. Published package unchanged (`ce78c81c…`, 2,595,584,000 bytes). The rebuild process was not
touched. No credentials, pods, or public-repository writes.

## Residual
The descriptor binding protects against pathname substitution and in-place edits during decode, but not
against a filesystem that returns different bytes for the same descriptor (not a realistic local case).
Concurrent bridge runs are isolated by `run_id` but still share the `tmp/` directory; the auditor is asked
to state whether that is sufficient.
