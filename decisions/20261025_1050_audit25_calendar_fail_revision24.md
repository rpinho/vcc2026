# 2026-10-25 10:50Z — AUDIT25 FAIL closed the calendar gap; revision 24 frozen; AUDIT26 launched

## Verdict read

AUDIT25 is **FAIL** on revision 23. The behavioral timestamp finding from
AUDIT23 is closed: 17/17 live and 17/17 frozen files matched the manifest, the
unchanged suite rejected both the delayed-capture mutant and the entry
exception, the real decoder retained `09:59:00Z`, matched 900/900 and cleaned
up, and all 245 frozen-suite cases passed.

The blocker is a real contract defect inherited from revision 22. The real
fresh-decode orchestrator validated only the digit shape of `started_utc` and
compared timestamps lexicographically. It accepted three calendar-impossible
bridge starts and returned 0 instead of 15:

- `2026-10-00T10:00:00Z`
- `2026-02-30T10:00:00Z`
- `2026-10-25T00:00:99Z`

The earlier receipt gate already rejects these values. AUDIT25 did not
demonstrate a full-launcher bypass. The failure applies to the narrower
orchestration contract, and the report is preserved unedited as FAIL.

## Repair

Revision 24 makes the timestamp contract explicit at both layers:

- `run_decode_qualification.py` parses the bridge start with
  `datetime.strptime(..., '%Y-%m-%dT%H:%M:%SZ')` and refuses an invalid
  calendar instant before deleting a prior result or entering the decoder.
- `check_decode_qualification.py` applies the same real-calendar parse to both
  timestamps and compares parsed datetimes rather than strings.
- `test_decode_qualification.py` grows from 15/15 to **19/19** with zero-day,
  February 30, invalid-second and invalid decode-calendar cases.
- `test_decode_orchestration.py` grows from 6/6 to **9/9** with the same three
  bridge counterexamples and an entry marker that must remain absent.

The live reproduction `package_bridge_1025/revision24_repro.py` exercised both
direct and isolated-launcher-tail paths: **8/8** cases passed. The valid
control returned 0 and entered the decoder; all three invalid values returned
15 with no entry.

The full suite matrix completed with exit 0: bridge 7/7, audit findings 60/60,
audit gate 43/43, receipt gate 35/35 ×3, qualification 19/19, orchestration
9/9 ×3, behavioral timestamp 3/3 ×3. The decoder bytes remain unchanged
(`c407ff0b75050e9d6ecf0477498a361085b293f5fceb68a4e7aa1c8e24839cc3`).

## Freeze

- Manifest: `package_bridge_1025/audit26_scratch/revision24_manifest.json`
- SHA-256: `3d4c833dbf6b1487d4094f8b1f7abd437ed85e3517b54af10a8a6b488d68b4c0`
- Scope: 16 bridge files plus one rebuild support file
- Identity: **17/17 live and 17/17 frozen** copies match
- Supersedes: revision 23

## Audit

AUDIT26 is running with a brief that names only the
revision-24 tree and requires the three exact calendar counterexamples. It is
read-only. A PASS must name these bytes; an exit code is not a result.

## Why this is not test-shopping

The failure is an explicit contract requirement from AUDIT25's brief. The
repair changes only the validation needed to meet it and adds adversarial
coverage for the exact observed values. The audit criterion is unchanged, the
earlier receipt gate is untouched, and the new manifest is audited by a fresh
independent report. AUDIT24 and AUDIT25 remain FAIL in the record.

## Safety

Nothing submitted, uploaded, replaced or decoded from the production package.
Safety entry `5NeeoWn96lfm23hgExcm`; package remains 2,595,584,000 bytes with
SHA-256 `ce78c81c453f74676ea0e6d2e1298cccf58aa0c2d6649f930b428b62a4939dc4`.
