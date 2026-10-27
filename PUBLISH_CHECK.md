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
directories. It now also rejects files larger than 5 MB. The staged-file sizes
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

## October 26 afternoon/evening refresh

This job curates existing calibration, sign, decoupled-amplitude, R9, MSE
instrumentation and biology studies, plus AI advisor briefs/answers and the
dated competitor-ideas document. It adds a retrospective, a second research
update and an unpublished thread draft. Exact scope and normalization are in
[EXPORT_20261026b.md](EXPORT_20261026b.md). The prior rules review above is
historical; no fresh rules clearance or new official score is claimed.

Before copying, free disk was 34.11 GiB, above the required 30.5 GiB. Local
commands use nice 10. No datasets, model/prediction arrays, fitted coefficient
arrays, raw competitor payloads, process listings, logs or error files are
published. Explicitly curated CSV/JSON diagnostics are below 1 MB each and
enumerated in tools/published_result_files.txt; all files are capped at 5 MB.
No real credentials, other people's private messages or personal knowledge
files were used as publication inputs. Workbench paths were normalized.

The new source is a partial historical release. Hashes in old result receipts
refer to original private-run inputs/sources; they do not certify edited export
files. Historical numerical PASS results are described as such. This job does
not rerun data-dependent fits, rendering, official scoring or challenge packaging.
All public Python sources are parsed with the reference Python 3.12 environment.
An initial system-Python 3.9 syntax finding was a version mismatch; the existing
source was retained without modification.

Final local verification is recorded below after the tree is frozen. The owner
also receives an external review inventory. No files are staged or committed;
no push, remote creation, submission or posting occurs.

### Final verification for this refresh

| Check | Result |
|---|---|
| Offline export/privacy checker, Python 3.12 | PASS: 716 files, 540 Python sources, zero findings |
| Gitleaks directory scan, redacted | PASS: zero leaks; repeated with inherited configuration removed and allow comments ignored |
| Synthetic export-checker policy probes | PASS: 17 cases, including size caps, explicit result list, paths, credentials, datasets, logs and symlinks |
| Curated JSON/CSV structure and Git visibility | PASS: 59 enumerated results, including 16 rectangular CSVs; no export file hidden by ignore rules |
| Relative Markdown file links in added/changed documents | PASS |
| Retrospective length | 1,976 whitespace-delimited words |
| Thread length | 9 posts; longest literal post 249 characters, including numbering and URL |
| Manifest | Regenerated with the original sorted SHA-256 procedure; 714 entries, independently checked for complete coverage and matching hashes |
| Whitespace and index | PASS: diff whitespace check; no staged changes |

The requested broad text scans find zero personal home-path prefixes, personal
mail-domain matches or forbidden provider names. Broad key terms and the short
credential-prefix fragment occur only in ordinary code/prose: dictionary/cache
keys, CLI mask options, disk/task/risk hyphenations and the checker/ignore rules.
Credential-shaped contexts were inspected; none contains a credential or private
key path. These benign text matches are not secret findings and no secret was
allowlisted. All true path/privacy findings in the initial copied material were
removed or normalized before the passing scans.

These are source publication checks, not new numerical campaign results or a
reproduction of official validation. The source workbench was read only; no
data/model artifacts were opened for testing. Final review is left to the owner.

## October 27 trim

The bulk decision record and competitor digests added on October 26 (commit c609890) were removed from `main`
to keep the repository to measured results and code; they remain in Git history. DECISIONS_INDEX.md now lists
the main rulings in one short table.
