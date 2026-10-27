# 2026-10-25 12:35Z — audit29 false-PASS repair and final verifier freeze

## Trigger

Independent audit 29 returned **FAIL** for the verifier produced by the first
repair. The archive evidence itself remained valid: both genuine archives passed,
and all 14 frozen audit28 probes were rejected. Ten new malformed cases were
accepted:

- duplicate `run_id`, `files`, or `status` JSON keys;
- an arbitrary 900-name receipt-hash mapping with zero digests;
- a symlink receipt and a hardlink receipt outside the archive;
- an unlisted subdirectory;
- whitespace or case changes in the run-ID text;
- `bridge exit 01`, which a substring check treated as success;
- a failed current run followed by a later successful run;
- a nonempty summary missing an ordinary field.

The audit also identified narrower unchecked-claim and trust-root gaps: state
timing/package identity, receipt notes and byte counts, the actual revision
manifest, AUDIT27's own digest, duplicate or unsafe source rows, and symlinked
archive paths.

## Repair

Only `package_bridge_1025/verify_production_closeout.py` and
`package_bridge_1025/test_verify_production_closeout.py` changed. No production
receipt, archive, log, manifest, AUDIT27 content, package, or prior audit report
was modified.

The verifier now:

- rejects duplicate JSON keys recursively before semantic evaluation;
- rejects symlinked archive directories or ancestors, and requires every evidence
  file and the manifest to be private regular files with link count 1;
- requires exact canonical run-ID text (`YYYYMMDDTHHMMSSZ-<32 lowercase hex>`);
- bounds the current log segment at its chain exit or next run start, requires one
  full decision summary marker and no duplicate start phrase, and parses exits exactly;
- requires all 900 canonical `group_{0..2}_{000..299}.json` names and verifies
  each against live bytes;
- binds receipt/state/decode package paths, digest, timestamps, ordering, notes,
  and the exact 2,608,164,839 decoded/temporary byte count;
- pins AUDIT27's SHA-256 and the actual revision-25 manifest SHA-256;
- rejects duplicate or unsafe AUDIT27 source-row paths.

## Verification

- Frozen verifier SHA-256:
  `3fb85315ab3cd8d18c30eba75ee2b86d0ebc551b8d5176b5d98966a2c445a759`.
- Frozen suite SHA-256:
  `1af0b07953f75bcc168e5b410b585362579a2ab66d58c61289d8003fafa3f7df`.
- Mutation suite: **44/44**, exit 0, including both genuine archives and every
  demonstrated false-PASS shape.
- Both genuine archives: **PASS**.
- Audit28 probes: **14/14** reject for intended reasons.
- Audit29 probes: **24/24** reject for intended reasons, with symlink and hardlink
  topology preserved during replay.

## Decision

Freeze these bytes and launch a fresh independent audit, `vcc-closeout30`, rather
than treating the stale audit29 run as authorization. Audit30 must use the frozen
hashes, replay both probe sets, add new adversarial cases, verify the suite and
false-rejection controls, and write `CLOSEOUT_AUDIT30.md`.

No submission, upload, package replacement, or new scientific arm. Safety entry
`5NeeoWn96lfm23hgExcm`; package `2,595,584,000` bytes, SHA-256
`ce78c81c453f74676ea0e6d2e1298cccf58aa0c2d6649f930b428b62a4939dc4`.
