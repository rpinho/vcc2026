# Publication check — refreshed 2026-10-23

This source-only refresh was explicitly authorized for publication. It extends
our 2026-10-18 release with existing response-MLP, v7/v9 layout, final builder,
checked packaging and official-scorer ablation code. No new official improvement
or complete hidden-score reproduction is claimed.

## Rules and data boundary

Rechecked the live [rules](https://virtualcellchallenge.org/rules),
[FAQ](https://virtualcellchallenge.org/faq), [terms](https://virtualcellchallenge.org/terms)
and [challenge home](https://virtualcellchallenge.org/) during the local
2026-10-23 review (UTC 2026-10-24). The public Next.js rules and FAQ content was
read from the page's current served chunks:

- Rules: `2622-63c808ba4137a5e1.js`.
- FAQ: `334-142922380196709a.js`.
- Terms: `5584-2d1dd1b066f1342c.js`.

The FAQ encourages research publication subject to institutional permission.
The rules retain team eligibility and prohibit cross-team arrangements that
undermine challenge integrity. No blanket code-publication embargo was found;
this interpretation is not organizer approval of an arrangement. Public issues,
forks and general contributions do not change registered team membership.

The website terms do not establish blanket permission to redistribute source
screens, challenge controls or derived per-cell data. This repository publishes
**no datasets, weights, predictions, source banks or .vcc archives**. Readers
obtain authorized inputs themselves. MIT covers our code, not external data.
Original upstream notices are retained; see [THIRD_PARTY.md](THIRD_PARTY.md).

Final scores remain hidden. The recorded current last entry is
`5NeeoWn96lfm23hgExcm`, the billmA v9 layout safety recipe. This refresh does
not run a submission command or change that entry. The deadline remains
November 5, 2026, 23:59 UTC.

## Export and portability review

Files were selected individually from the private workbench. No private Git
history, notes, logs, launchd/agent supervisors, journals, raw leaderboard/member
payloads, credentials or account stores were copied. Our official aggregate score
history and a newly written explanation of lessons are the only result summaries.
Final entry IDs are public identifiers, not authentication tokens.

AtlasShift and Sidechain vendor source retains its notices. The independently
written numerical checks are in validation/. No new stage1-lab, numba-mwu,
State, Stack or PIE source/checkpoints are included. External installed packages
retain their licenses. Optional PIE data must be prepared separately.

Publication edits discover vcc on PATH, remove unsupported historical --force
from local prep commands, remove machine workload hooks/status-journal writes,
make code-module discovery independent of scratch-directory depth, and create
missing output/adapter directories on a fresh checkout. The numerical recipes
and fail-closed count/moment checks are retained. One truncated, unused tier
prototype was omitted. Historical receipts must be regenerated: source hashes
from the private checkout do not certify this export.

## Verification scope

Checks run in a disposable source copy, with synthetic data only. No production
weights/data were opened for testing; no prediction package was built or
submitted. Shared environment: Python 3.12 and cell-eval2 0.18.0.

| Check | Result |
|---|---|
| Exported Python syntax | PASS |
| Final safety shell syntax | PASS |
| Panel builder, panel packager, direction-knob builder CLI startup | PASS |
| Neural training and all three official-scorer ablation CLI startups | PASS |
| Synthetic neural recipient-family/sister-screen exclusion, measured masks, finite gradients and unit-direction blend | PASS |
| Feasible bulk projection against separate scalar solver | PASS |
| Regenerated bounded DE and bounded-scoring fixtures | PASS |
| Existing full adapter against installed official scorer | PASS |

The newly released v9 rank-test and independent count kernels are also checked
on synthetic tied count data: **PASS**, matching SciPy asymptotic MWU on
eligible genes and independently checking jackknife, row/gene totals and moments.
All-zero ineligible genes use the implementation's neutral p=1 convention.
Historical full-panel training, all-group package reconstruction and hidden
biological performance are not reproduced by these small checks.

## Required source/privacy checks

Run on the clean source tree before publication:

```sh
env -u GITLEAKS_CONFIG -u GITLEAKS_CONFIG_TOML gitleaks dir . \
  --no-banner --no-color --redact --ignore-gitleaks-allow
python tools/check_export.py
git diff --cached --check
```

The offline checker checks syntax, personal home paths, emails, credential
literals, private keys, binary/unexpected file types, symlinks and runtime/data
directories. It now also rejects files larger than 20 MB. The staged-file sizes
and Git tree are checked separately. SOURCE_MANIFEST.sha256 is refreshed for
source files only, excluding this audit and the manifest itself. No large input,
weight or package is checksummed as part of this publication job.

Final scan results are recorded after the source tree is complete. Any privacy
or secret finding blocks publication; none may be silently allowlisted.

## October 26 source refresh

This refresh selects five analysis/scoring modules from a fixed workbench
snapshot and publishes a newly written aggregate result summary. It preserves
the previous public source export and its license exclusions. Source normalization replaces a personal interpreter path with `python`.
Spaces around existing NumPy matrix products prevent email-detector false
positives without changing their calculations. No private
Git history, model/data arrays, receipts, raw responses or operational notes
are transferred. Official validation is unchanged; local SEARCH scores are
labelled separately. The public score ledger now includes the entries already
present in the README.

The 2026-10-23 publication-rules review above is the last substantive rules
review; the October 26 web reader returned only the site's loading shell. This
refresh does not assert a new organizer clearance or change participation.

Required publication checks: full source-tree secret and private-term scans,
outgoing commit/diff scans, source/privacy checker, source manifest, diff
whitespace, and synthetic aggregation/layout-bound checks. No full historical
campaign, private dataset, package build or challenge submission is run. The
selected source release is explicitly partial; its absent inputs must be
regenerated rather than bypassing provenance checks.
