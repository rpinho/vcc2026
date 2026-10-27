# DECISION 2026-10-25 05:05 ET — family lever closes on the budget gate

## Verdict before the third arm finished
`family_quality` fails the frozen H1 retention gate on **budget match**, independent of how
`family_equal` scores and independent of the bootstrap:

| quantity | A0 incumbent | A1 family_quality | frozen limit |
|---|---:|---:|---:|
| requested B/claim | 0.9090700116298247 | 0.9600968725138042 | — |
| absolute difference | — | **0.05102686088397945** | **0.020** |

The driver wrote this itself to `artifacts/H1_rich_search/budget_match_provisional.json` with
`budget_match: false` and the action "no CD4 because the H1 budget-match gate fails". This is exactly the
defect that invalidated the archived attempt (0.8673); re-calibration did not remove it.

## Point result, for the record only
The scored H1 arms are:

| arm | official overall | requested ratio | all required gates |
|---|---:|---:|---|
| `incumbent` | −0.01916405775558278 | 0.9090700116298247 | pass |
| `family_quality` | −0.0023418631033149256 | 0.9600968725138042 | pass |
| `family_equal` | pending at 31/79 when this was written | pending | pending |

`family_quality` is about **+0.0168 official overall** over the incumbent and raises PDS from 0.5298 to
0.6092 in raw terms, while MSE/NMAE/fidelity/reach/Jaccard move only slightly. That is a real local signal,
not a scorer artefact. It is **not promotable** as scored because the budget moved 2.5x past the
pre-registered tolerance, and the arm cannot be re-tuned after seeing the result without destroying the
test's validity.

## What this means scientifically
- At fixed `B = claim`, a direction replacement changes the jackknife denominator enough that the
  requested ratio shifts by 0.051. Any future direction candidate must be calibrated **against the
  final emitted layout**, not merely against the pre-blend profile, or it will fail this same check.
- The observed +0.0168 is direction quality, but the comparison is between two different *amplitudes* as
  well as two directions. This run does not isolate direction from the calibration response.
- `family_equal` remains a diagnostic. It cannot be promoted even if it wins; the amendment says so.

## Consequence
- **The family-direction lever is CLOSED for this campaign.** No CD4 continuation: the watcher requires all
  five gates, sees `budget_match: false`, and will exit without launching. That is correct.
- No fourth arm, no second anchor, no post-hoc rescaling to force the ratio into range. Re-solving amplitude
  after the blend was already measured at −0.041 scaled MSE elsewhere; repeating that here would be a
  different experiment chosen after seeing this result.
- Both independent consultants (`consult_0430/`, decision `20261025_0440`) said to stop after the registered
  programme and hand off. This result supplies the final piece of that programme.

## Pending, explicitly
- `family_equal` score and all independent verification/bootstrap/component artifacts were not present
  when this file was written. The report must append them rather than silently omitting them.
- No claim of official-score improvement. No package was built, uploaded, replaced or submitted.

Safety entry `5NeeoWn96lfm23hgExcm` and package
`ce78c81c453f74676ea0e6d2e1298cccf58aa0c2d6649f930b428b62a4939dc4` are unchanged.
