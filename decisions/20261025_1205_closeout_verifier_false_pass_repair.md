# 2026-10-25 12:05Z — closeout verifier false-PASS repair

> **Correction, 12:35Z.** The figures below describe the first repair. Independent
> audit 29 then demonstrated ten additional false PASSes. The verifier and suite
> were repaired again and are frozen at verifier
> `3fb85315ab3cd8d18c30eba75ee2b86d0ebc551b8d5176b5d98966a2c445a759`
> and suite `1af0b07953f75bcc168e5b410b585362579a2ab66d58c61289d8003fafa3f7df`.
> The final suite is **44/44**, all 14 audit28 probes and all 24 audit29 probes
> reject for intended reasons, and both genuine archives pass. See
> `decisions/20261025_1235_closeout_audit29_repair.md`. The original text is kept
> as the historical record of the first repair, not as the current state.

## Trigger

The independent read-only audit `vcc-closeout28` returned
`Verdict: **FAIL**` in `package_bridge_1025/CLOSEOUT_AUDIT28.md`. The fresh
revision-25 archive itself passed all **155/155** independent reconstruction,
semantic, source-identity, and isolation checks. The failure was in
`verify_production_closeout.py`, which accepted malformed evidence:

- empty returned summary;
- manifest omitting one or all evidence files;
- extra manifest entry;
- wrong manifest run ID;
- wrong run ID in the production log;
- extra unlisted file;
- wrong decode snapshot digest;
- altered structure values;
- invalid decode timestamp;
- unknown/null returned-summary key;
- nonzero exits rewritten and rehashed.

The supplied 11/11 suite used only the earlier archive, so it did not cover the
new run or these malformed shapes.

## Repair

Only `verify_production_closeout.py` and its test were changed. Production
archives, receipts, logs, AUDIT27, the revision-25 manifest, and the safety
package were not modified.

The verifier now requires:

- exactly `{run_id, files}` in the manifest;
- `manifest.run_id` equal to the archive directory and to the run-ID file;
- exactly the six named evidence files, each with exact byte count and SHA-256;
- no extra or missing files in the archive directory;
- a nonempty returned summary with no unknown keys;
- receipt-minus-summary exactly the three intentionally filtered fields;
- exact receipt structure and archived-structure values;
- exact decode snapshot identity, valid UTC timestamps, and real datetime
  ordering;
- current-run log segment with `bridge exit 0`, `production chain exit 0`, and
  exactly one occurrence of the current run ID.

## Verification

- Repaired verifier: **PASS** on both preserved archives.
- Mutation suite at that time: **24/24**, with baseline coverage of both archives and all
  previously demonstrated false-PASS cases.
- Replay of all 12 original malformed audit probe directories: every one now
  exits 1 for its intended reason.
- `CLOSEOUT_AUDIT28.md` remains unchanged as the frozen pre-fix FAIL report.

## Decision

Launch a second independent audit (`vcc-closeout29`) against the patched
verifier and expanded suite. It must replay the original probe directories,
add its own independent malformed cases, and inspect whether the new checks
introduce false rejection or self-referential trust. Do not treat the first
audit as retroactively passing, and do not edit its report.

No submission, upload, package replacement, or new scientific arm. Safety
entry `5NeeoWn96lfm23hgExcm`; package `2,595,584,000` bytes, SHA-256
`ce78c81c453f74676ea0e6d2e1298cccf58aa0c2d6649f930b428b62a4939dc4`.
