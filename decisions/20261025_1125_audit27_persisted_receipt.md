# 2026-10-25 11:25Z — persisted-receipt defect; revision 25 audit and fresh production run

## Defect

AUDIT26 passed revision 24, and the audited production chain then passed:
900/900 groups, zero mismatches, identical structure, unchanged safety hash,
and an independently qualified fresh decode. That run exposed an evidence
defect in the producer.

`bridge.py` set `receipt['output_directories_stable'] = True` **after**
`_write_secure_json` serialized the PASS receipt. The returned Python object
therefore contained the field, but the persisted production receipt did not.
The production receipt is historical evidence and was not edited.

## Repair and revision 25

- `bridge.py` now sets the field before serialization, with the reason stated
  in the source.
- `test_bridge_gates.py` now requires every PASS case's persisted JSON to equal
  the returned receipt exactly, not merely to have status PASS.
- A mutant preserving the old post-write ordering was reproduced separately;
  the repaired suite fails it at **6/7** with `persisted receipt differs from
  returned receipt`.
- The repaired source passes **7/7**. Full local matrix: bridge 7/7, findings
  60/60, audit gate 43/43, receipt gate 35/35 x3, qualification 19/19,
  orchestration 9/9 x3, timestamp 3/3 x3.
- Revision 25 was frozen at
  `package_bridge_1025/audit27_scratch/revision25_manifest.json`, SHA-256
  `72dbe540bc88a03e8f6f2dc16a8063a1d8b6cac99f0e7fb7d9d007248e5ae710`.
  It contains 16 bridge files and one support file; all **34** live/frozen
  copies match their manifest byte counts and SHA-256 values.

Only `bridge.py` and `test_bridge_gates.py` differ from revision 24.

## Consultation

Question: before 08:00 ET, rank (1) finish AUDIT27 and rerun production,
(2) reconcile retained local/official artifacts, (3) a new scientific arm, or
(4) something else.

- **Kimi** (`orq best`, `consult_0718/kimi.md`): **1 -> 4 -> 2 -> 3**.
  Verification plus decision summary is the only defensible priority. AUDIT27 must
  name the exact manifest; production must use a fresh run ID; the persisted
  receipt must equal the returned receipt with a boolean
  `output_directories_stable: true`; bridge and fresh decode must each pass
  900/900; the package must remain byte-identical. Any failure, stale
  evidence, or incompleteness at 08:00 is a reported gap, not PASS.
- **Astra** (`codex exec -m gpt-6-astra`, `consult_0718/astra.md`): **1 -> 4 ->
  2 -> 3**. Same gate. Defer retained-artifact reconciliation unless a concrete
  unresolved question and a predeclared decision rule exist; reject a new arm
  or sweep because the registered levers are closed.
- **Opus 5.5** (`orq anthropic/claude-opus-5.5`, `consult_0718/opus.md`):
  **1 -> 4 -> 2 -> 3**. Same gate. Reject option 2 for tonight because another
  fit or bootstrap of the same adaptations cannot supply the missing exact
  matched observation.

All three were independent and agreed. Their acceptance gates match.

## Decision

1. Launch and read AUDIT27 against exactly revision 25. Do not infer PASS from
   the process exit code.
2. If AUDIT27 is PASS, rerun the audited production chain with a fresh run ID,
   read the **persisted** receipt from disk, require exact equality with the
   returned receipt and a boolean `output_directories_stable: true`, and
   require a fresh independent decode to qualify the same package identity.
3. Preserve the revision-24 production PASS as historical evidence; it does not
   authorize revision-25 bytes.
4. Do not launch a new scientific arm, sweep, package replacement, upload or
   submission. Prepare the morning decision summary with exact identities and any
   unfinished verification marked pending.
5. If AUDIT27 fails, preserve the report, repair only the named defect, freeze
   a new revision, and re-audit. Never edit the production receipt or weaken a
   test to obtain PASS.

## Safety

Nothing submitted, uploaded or replaced. Safety entry
`5NeeoWn96lfm23hgExcm`; package 2,595,584,000 bytes, SHA-256
`ce78c81c453f74676ea0e6d2e1298cccf58aa0c2d6649f930b428b62a4939dc4`.

## Outcome and wrapper correction

AUDIT27 returned **PASS** on revision 25, naming manifest
`72dbe540bc88a03e8f6f2dc16a8063a1d8b6cac99f0e7fb7d9d007248e5ae710`. Its
independent 900-group fixture produced a persisted receipt equal to the returned
receipt, and the old post-write ordering failed the unchanged suite.

The first fresh production run,
`20261025T112331Z-2c3a6c7bfedc2642a44b9f3d3b06b2b2`, completed the audited
chain: bridge 900/900, zero mismatches, exact structure, unchanged safety hash;
fresh independent decode 900/900, shape `[360000, 18533]`, nnz
`2010253470`; receipt validator and qualification both exit 0. The persisted
receipt contains boolean `output_directories_stable: true`, and all 17 live
sources still match AUDIT27.

My post-run wrapper nevertheless logged `FAIL bridge returned-receipt capture
is missing`. That was a wrapper false negative, not a bridge failure: the shell
had already loaded the old stdout path, and the bridge CLI intentionally omits
three large fields from its stdout summary. Direct comparison showed every
returned field agrees with persisted evidence; the three extra persisted fields
are `archived_receipt_sha256`, `archived_structure`, and `structure`. The
wrapper now compares returned fields as a required subset and is being rerun
under a fresh run ID for a clean end-to-end record. The failed wrapper log is
retained unchanged as historical evidence of the harness defect.

## Preserved evidence and morning checks

- Archive:
  `package_bridge_1025/production_runs/20261025T112331Z-2c3a6c7bfedc2642a44b9f3d3b06b2b2/`
- AUDIT27: **PASS** on the revision-25 manifest; all 17 live and frozen sources
  matched before and after testing.
- Closeout verifier: **PASS**; its six-case mutation suite is **6/6**.

```sh
../.venv/bin/python -B package_bridge_1025/verify_production_closeout.py \
  package_bridge_1025/production_runs/20261025T112331Z-2c3a6c7bfedc2642a44b9f3d3b06b2b2
../.venv/bin/python -B package_bridge_1025/test_verify_production_closeout.py
```
