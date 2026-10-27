# 2026-10-25 08:04Z — AUDIT13 FAIL; all five findings fixed in revision 15

## AUDIT13 verdict: FAIL (real defects, 9/9 hashes matched, no drift)

AUDIT13 audited revision 13 and reproduced five issues. Four of them applied unchanged to revision 14,
so revision 14 was **not** a safe place to start the bridge either. Its suite-count critique was also
correct and worth recording: my "114/114" was a *live-tree* number, and it does not reproduce inside a
frozen snapshot, because `test_bridge_audit_findings.py` reaches for
`rebuild_repro_1025/finalize_after_receipt.sh` (outside the bridge directory) and
`test_bridge_receipt_gate.py` names the real production package. Both are scope artifacts, not gate
defects, but I should have stated the scope of the count.

| finding | what it was | fix in revision 15 |
|---|---|---|
| 1 | run-state counts were never checked, so `900.0`, `true`, `"900"` and duplicate keys passed in the state | state counts must be strict ints **and** equal the receipt's; duplicate keys rejected with an explicit `run state unreadable: duplicate key …` reason |
| 2 | `out.name == STATE_NAME` missed case aliases on this case-insensitive volume, so the receipt could overwrite the state | comparison is now `casefold()`-based **and** additionally uses `os.path.samefile` for existing aliases |
| 3 | early failures left a previous run's PASS state in place; a PASS receipt whose PASS-state write failed was still reported as success | rejections retire a leftover PASS state **only** in a provably safe output directory (never through the package inode or an ancestor), and the PASS-state write result is now checked |
| 4 | "900 self-consistent files" satisfied the archive binding, so `group_EVIL_000…899` passed | archive must match the exact 3×300 identity set with per-context `context`/`target`/`index` consistent with the filename |
| 5 | deeply nested JSON raised `RecursionError` outside the exception tuple → traceback | parsing bounded at 16 MiB and depth 64; `RecursionError`/`MemoryError` are clean rejections |

Two additional bugs were found *by the new tests while fixing these*, which is the point of writing them:
an unparsable archive filename (`group_EVIL_000.json`) raised `ValueError` from `int(parts[1])` instead of
failing closed, and a duplicate-key state degraded to `{}` so every field failed for the wrong reason
instead of reporting the real cause. Both fixed.

## What AUDIT13 correctly says is NOT closed

- A gate cannot **authenticate the provenance** of a reference archive, and cannot substitute for the
  bridge's actual decode. Pinning the archive to the real artifact tree is a trust-boundary statement,
  not a proof of authenticity. This is recorded as an explicit residual, not papered over.
- There is an unavoidable observation-to-publication window: bytes can change after the last hash. No
  number of pathname rehashes covers the future. The safety anchor plus the launcher's run identity are
  what make that window non-silent for production.

## AUDIT14: no verdict

AUDIT14 was launched against revision 14 and ended on a provider content-policy flag with no report —
the same failure mode as AUDIT12 and AUDIT6. A refusal is not a verdict, so it authorizes nothing.
Its frozen sources are revision 14, which is now superseded anyway.

## Current state

- Revision 15 frozen at `audit15_scratch/revision15_manifest.json`; AUDIT15 launched with a
  correctness-review brief (C1–C7) rather than an exploit-framed one, because the last two reviewers
  were cut off by policy rather than by evidence.
- Suites: bridge 7/7, findings 45/45, audit gate 43/43, receipt gate 27/27 = **122/122**, and the
  production control (real package + real archive + matching state) passes with zero failures.
- Rebuild chain remains closed PASS. Nothing submitted. Published package unchanged at
  `ce78c81c…9dc4`. Safety entry `5NeeoWn96lfm23hgExcm`.

## Rule I am keeping

Do not start the bridge on a verdict that does not name the exact bytes that would run. AUDIT13's PASS
would not have authorized revision 14 or 15, and AUDIT14's absence authorizes nothing.
