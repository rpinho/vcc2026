# 2026-10-25 14:10Z — independent reproducibility check: 6/6 headline claims reproduced, 3 overstatements corrected

## What the check did

`closeout_repro_check_1025/REPRO_CHECK.md` (job `vcc_closeout_reproducibility_1025`, exit 0) re-derived
tonight's headline claims from primary artifacts without trusting my reports. It found **all six
headline results REPRODUCED**:

| claim | measured | verdict |
|---|---|---|
| frozen suite identity | `5efd424e…b35c` | REPRODUCED |
| mutation suite | **54/54**, exit 0 | REPRODUCED |
| frozen verifier identity | `1f1df2b7…61d2` | REPRODUCED |
| both genuine archives | PASS, exit 0, both runs | REPRODUCED |
| bridge audit / receipt / decode gates | PASS, exit 0 (17 sources listed) | REPRODUCED |
| package G1-G8 | 298/300, alpha 1.092365/1.092389/1.134883, B/claim 1.0, mean L1 0.0020607, bulk L1 2.6e-08, `[360000,18533]`, 900/900 | REPRODUCED |
| inventory buckets | recounted 67 = 0/0/23/44 | REPRODUCED |
| safety package | 2,595,584,000 / `ce78c81c…9dc4`, same inode+mtime | REPRODUCED |

It also verified 109 referenced evidence-file hashes and 42 original-helper snapshot hashes, all
matching, and confirmed the `PACKAGE` defect is a false rejection (C13: unrebound exit 1, rebind → exit 0).

## Three overstatements it caught, all now corrected

**1. The 51/54 relocation figure was wrong for the frozen suite.** Executing the *frozen* suite
(`5efd424e…b35c`) with a relocated import path gives **52/54** — only `old_baseline` and
`new_baseline` fail. My 51/54 came from Kimi's scratch fixture, whose suite is a *different artifact*:
it changes `module.LIVE_ROOT = HERE` to `module.LIVE_ROOT = root` (fixture suite `7e799008…330a`),
which also makes `audit_hash_changed` fail for want of live source files. I independently confirmed
this diff. I had conflated the genuine `PACKAGE` binding defect with a separate fixture modification,
and overstated that a relocated 51/54 is "solely attributable to `PACKAGE`".
Erratum appended to `decisions/20261025_1335_suite_relocation_binding.md`.

**2. "No slot or daily allowance is blocking" did not follow from the two CLI commands I cited.**
`vcc status <id>` queries one submission; it cannot establish that no newer entry exists.
`vcc whoami` prints `can_submit` from `get_me`. Daily allowance comes from a separate
`GET /api/cli/submissions/limits` call that `submit.py` makes only during submission.
I verified this in the installed CLI source. Corrected in the decision summary and `NEXT_SUBMISSIONS.md`.

**3. "The evidence base does not exist and cannot be manufactured" overstated a bounded inventory.**
The inventory certifies 0 of 67 records under a four-field identity test and explicitly disclaims
full-repository absence. The operational conclusion — do not build a promotion rule — is a judgement
resting on that certification plus the deadline, not a number the inventory emits alone. Under a
looser procedural reading the v5a/v5b LOCO ports are *arguably* recipe-matched, though they remain
79-target ports with `external_official_label: null` and so still useless as official transfer pairs.
Corrected in `decisions/20261025_1430_matched_inventory_zero.md`, the decision summary, and `NEXT_SUBMISSIONS.md`.

## Also corrected

`package_bridge_1025/STATUS.md` mixed current and superseded states: its undated "Gate state" list and
per-revision sections still read as though AUDIT26 were running and the bridge were blocked. Added a
top banner stating that everything below it is chronology, that revision 25 is AUDIT27-authorized,
the verifier is AUDIT31-authorized, and nothing is blocked.

## Integrity note

The check's blanket "all 1,136 inputs unchanged" assertion failed (1,135/1,136) because **I appended my
own log entry mid-run**. It correctly preserved the failure as C14 exit 1 rather than erasing it, and
verified the original 139,330-byte log prefix still hashed to its recorded value. Lesson: do not write
to files a running job is fingerprinting.

## Outcome

The verification work stands — every substantive claim reproduced independently. What changed is
precision in three places where I had let a claim run ahead of its evidence. The recommendation is
unchanged: retain `5NeeoWn96lfm23hgExcm`; do not submit.

No submission, upload, or package replacement.
