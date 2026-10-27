# Decision — overnight re-audit of bridge revision 3

Date: 2026-10-26 07:00 UTC / 2026-10-25 03:00 EDT
Scope: rebuild bridge (bridge.py + gates), not the VCC model or submission.

## Context
The previous bridge audit (AUDIT2, `package_bridge_1025/AUDIT2.md`) failed. It found
two classes of false-pass risk: (1) in-place package mutation/restoration could leave
the verifier hashing A while decoding B, and (2) receipt/state replay or audit-gate
hash mixing could make stale evidence look current. Additional gaps existed in
post-check symlink handling, disk reserve arithmetic, and hash/parse races on group
receipts.

The user's instruction is explicit: keep progress going, never submit, and consult
Kimi/Astra when stuck. This was not a stuck decision; it was a mandatory
fail-closed correction of an audit finding, so no external consultation was needed.
The change is within the verification-only lane.

## Decision
Do not run the bridge until a fresh independent audit returns PASS against a frozen
revision. Revision 3 was built to make the verifier's claims structural:

- private immutable snapshot used for all decoding;
- four-point public package hash binding (pre-copy, post-copy, post-decode, pre-publish);
- atomic O_NOFOLLOW publication of receipt/state;
- launcher-supplied run ID, checked by the receipt gate;
- one-read evidence files with hash and parsed bytes from the same buffer;
- disk reserve subtracts the incoming block;
- strict two-column audit hash rows.

All 53 regression cases pass. The audit run itself was launched twice due to
environmental mistakes (missing executable bit, missing prompt path, unsupported
model name); none of those touched the package or rebuild.

## Rejected alternatives
- Treating AUDIT2 as advisory: rejected; its in-place mutation case is a genuine
  false PASS, and the bridge's purpose is evidence integrity.
- Adding more pathname checks without a snapshot: rejected as insufficient.
- Editing the audit report by hand: rejected; a new independent run is required.

## Safety
Nothing submitted. Published package was not opened or modified by any test. The
read-only audit is instructed to avoid the real package.

## Follow-up
1. Read audit 3 verdict when it lands.
2. If PASS, wire the exact audit hash into the bridge runner and let the bridge
   decode proceed.
3. If FAIL, fail closed again and repair only the reproduced findings.
4. Continue monitoring rebuild and revisions 7 audit; nothing in this lane changes
   the submission decision.
